*! hsmixture_certification_bivariate.do
*! Parameter-recovery certification for hsmixture_bivariate.
*!
*! The bivariate command admits a 2x2 grid of (v_T, v_Y) corners with a
*! free joint probability matrix. As the off-diagonal probabilities go to
*! zero the model approaches factor(separate) (two free diagonal shifts);
*! it matches factor(common) only if v_T2 = v_Y2. The softmax keeps every
*! cell strictly positive, so this is a limiting case. This script
*! exercises the off-diagonals.
*!
*! DGP rationale:
*!   - True structure has substantial mass on all four corners. The joint
*!     factor(separate) estimator (which constrains types to a 1-D locus)
*!     cannot fit this exactly; the bivariate can.
*!   - v_T2 = 1.5  (treatment-prone shift)
*!   - v_Y2 = 1.0  (outcome-prone shift)
*!   - pi_11 = 0.40 (low-T, low-Y)
*!   - pi_12 = 0.20 (low-T, high-Y)   <-- off-diagonal
*!   - pi_21 = 0.20 (high-T, low-Y)   <-- off-diagonal
*!   - pi_22 = 0.20 (high-T, high-Y)
*!   - delta = 0.5  (treatment raises outcome hazard ~65%)
*!   - N_persons = 6000 (4 corners need more events per cell than 2 mass points)
*!
*! Pass criteria:
*!   - e(converged) == 1 (BFGS converged AND |grad|/(1+|LL|) < 1e-5 AND V is PD)
*!   - |delta_hat - 0.5| < 3 * SE(delta)
*!   - |v_T2_hat| > 0.5  (treatment heterogeneity recovered, sign may flip
*!                       under labeling swap)
*!   - |v_Y2_hat| > 0.5  (outcome heterogeneity recovered)
*!   - Off-diagonal mass present: min(pi_12, pi_21) > 0.05
*!     (the bivariate is doing more than the 1-D-locus joint fit)
*!   - Implied rho recovered up to swap symmetry. Truth pi structure has:
*!       E[v_T] = (pi_21 + pi_22) * v_T2 = 0.4 * 1.5 = 0.60
*!       E[v_Y] = (pi_12 + pi_22) * v_Y2 = 0.4 * 1.0 = 0.40
*!       Cov   = pi_22 * v_T2 * v_Y2 - E[v_T]*E[v_Y] = 0.30 - 0.24 = 0.06
*!       Var_T = (pi_21 + pi_22) * v_T2^2 - E[v_T]^2 = 0.90 - 0.36 = 0.54
*!       Var_Y = (pi_12 + pi_22) * v_Y2^2 - E[v_Y]^2 = 0.40 - 0.16 = 0.24
*!       rho = 0.06 / sqrt(0.54 * 0.24) = 0.06 / 0.360 = 0.167
*!     With labeling-swap symmetry the absolute value is what matters.
*!     Pass if |rho_hat| within [0.05, 0.45].
*!
*! Authors: Jonghoon Park and R. Alan Seals (Auburn University)

clear all
set seed 20260507
set more off

capture log close _hsmix_biv
display as txt "Working directory: `c(pwd)'"
log using "hsmixture_certification_bivariate_log.txt", text replace name(_hsmix_biv)

* ============================================================================
* PART 1: Generate synthetic person-period data with bivariate-grid DGP
* ============================================================================

display _n as txt "{hline 70}"
display as txt "hsmixture_bivariate certification: 2x2 grid recovery test"
display as txt "{hline 70}"
display as txt "Generating synthetic data with bivariate-grid heterogeneity..."

local N_persons = 6000
local T_max     = 20

local delta_true = 0.5
local v_T2_true  = 1.5
local v_Y2_true  = 1.0
local pi_11_true = 0.40
local pi_12_true = 0.20
local pi_21_true = 0.20
local pi_22_true = 0.20
local beta_x1_true = 0.3

quietly {
    set obs `N_persons'
    gen long id = _n

    gen double x1 = rnormal()

    * Assign joint type from the four-cell grid
    gen double u_type = runiform()
    gen byte type_T = 0
    gen byte type_Y = 0
    replace type_T = 0
    replace type_Y = 1 if u_type < `pi_12_true'
    replace type_T = 1 if u_type >= `pi_12_true' & u_type < `pi_12_true' + `pi_21_true'
    replace type_Y = 0 if u_type >= `pi_12_true' & u_type < `pi_12_true' + `pi_21_true'
    replace type_T = 1 if u_type >= `pi_12_true' + `pi_21_true' & u_type < `pi_12_true' + `pi_21_true' + `pi_22_true'
    replace type_Y = 1 if u_type >= `pi_12_true' + `pi_21_true' & u_type < `pi_12_true' + `pi_21_true' + `pi_22_true'

    gen double v_T = type_T * `v_T2_true'
    gen double v_Y = type_Y * `v_Y2_true'

    expand `T_max'
    bysort id: gen int period = _n

    * Treatment hazard
    gen double eta_T = -3.0 + 0.02 * period + `beta_x1_true' * x1 + v_T
    gen double h_T = 1 - exp(-exp(eta_T))

    gen double u_T = runiform()
    gen byte treat_event = 0
    gen treat_time = .

    forvalues t = 1/`T_max' {
        bysort id (period): replace treat_event = (u_T < h_T) ///
            if period == `t' & treat_time == .
        bysort id (period): replace treat_time = `t' ///
            if treat_event == 1 & period == `t' & treat_time == .
    }

    bysort id: egen _tt = min(treat_time)
    replace treat_time = _tt
    drop _tt

    gen byte treated = (period > treat_time) if treat_time != .
    replace treated = 0 if treat_time == .

    replace treat_event = 0
    replace treat_event = 1 if period == treat_time

    * Outcome hazard
    gen double eta_Y = -3.0 + 0.02 * period + `beta_x1_true' * x1 ///
        + `delta_true' * treated + v_Y
    gen double h_Y = 1 - exp(-exp(eta_Y))

    gen double u_Y = runiform()
    gen byte outcome_event = 0
    gen outcome_time = .

    forvalues t = 1/`T_max' {
        bysort id (period): replace outcome_event = (u_Y < h_Y) ///
            if period == `t' & outcome_time == .
        bysort id (period): replace outcome_time = `t' ///
            if outcome_event == 1 & period == `t' & outcome_time == .
    }

    bysort id: egen _ot = min(outcome_time)
    replace outcome_time = _ot
    drop _ot

    replace outcome_event = 0
    replace outcome_event = 1 if period == outcome_time

    gen int censor_time = cond(outcome_time != ., outcome_time, `T_max')
    keep if period <= censor_time

    gen byte treat_at_risk = (treated == 0)

    drop eta_T h_T u_T eta_Y h_Y u_Y treat_time outcome_time censor_time ///
        type_T type_Y v_T v_Y u_type

    tab period, gen(pd_)
    drop pd_1
}

quietly count if treat_event == 1
local n_treat_events = r(N)
quietly count if outcome_event == 1
local n_outcome_events = r(N)
quietly tab id
local n_persons_observed = r(r)

display _n as txt "Realized event counts:"
display as txt "  Persons observed:" _col(35) as res `n_persons_observed'
display as txt "  Treatment events:" _col(35) as res `n_treat_events'
display as txt "  Outcome events:"   _col(35) as res `n_outcome_events'
display as txt "True DGP: delta=`delta_true', v_T2=`v_T2_true', v_Y2=`v_Y2_true'"
display as txt "  pi joint = (.40, .20 / .20, .20)"
display as txt "  Implied rho = " %5.3f 0.06 / sqrt(0.54 * 0.24)

* ============================================================================
* PART 2: Fit hsmixture_bivariate and assert recovery
* ============================================================================

display _n as txt "{hline 70}"
display as txt "Fitting hsmixture_bivariate (2x2 grid)..."
display as txt "{hline 70}"

hsmixture_bivariate (treat_event = x1 pd_*) ///
    (outcome_event = x1 pd_*, treat(treated)) ///
    , id(id) iterate(500) riskset(treat_at_risk) nstarts(6)

local conv      = e(converged)
local gn        = e(grad_norm)
local rg        = e(rel_grad)
local pd        = e(v_pd)
local delta_hat = e(delta)
local v_T2_hat  = e(v_T2)
local v_Y2_hat  = e(v_Y2)
local rho_hat   = e(rho)
local se_delta  = e(se_delta)

* Off-diagonal mass from joint probability matrix
matrix pimat = e(pi_joint)
local pi_12_hat = pimat[1, 2]
local pi_21_hat = pimat[2, 1]
local off_diag_min = min(`pi_12_hat', `pi_21_hat')

* Pass criteria
local pass_conv  = (`conv' == 1)
local pass_delta = abs(`delta_hat' - `delta_true') < 3 * `se_delta'
* Mass points magnitudes recovered (sign may flip under swap symmetry)
local pass_v_T   = abs(`v_T2_hat') > 0.5
local pass_v_Y   = abs(`v_Y2_hat') > 0.5
* Off-diagonal probability mass: at least one of (pi_12, pi_21) substantively
* nonzero. Under strict label-swap, true off-diag is min(0.20, 0.20) = 0.20.
local pass_offd  = (`off_diag_min' > 0.05)
* Implied correlation magnitude in plausible range. Truth |rho| ~ 0.17.
local pass_rho   = (abs(`rho_hat') > 0.05) & (abs(`rho_hat') < 0.45)

display _n as txt "{hline 70}"
display as txt "CERTIFICATION RESULTS (hsmixture_bivariate, 2x2 grid)"
display as txt "{hline 70}"

local status_delta = cond(`pass_delta', "{res}PASS", "{err}FAIL")
display as txt "  delta (treatment effect)" _col(35) %8.4f `delta_hat' ///
    "  SE " %6.4f `se_delta' "  truth " %5.2f `delta_true' "  `status_delta'"

local status_vT = cond(`pass_v_T', "{res}PASS", "{err}FAIL")
display as txt "  |v_T2| (truth 1.50)" _col(35) %8.4f abs(`v_T2_hat') ///
    "  threshold > 0.50  `status_vT'"

local status_vY = cond(`pass_v_Y', "{res}PASS", "{err}FAIL")
display as txt "  |v_Y2| (truth 1.00)" _col(35) %8.4f abs(`v_Y2_hat') ///
    "  threshold > 0.50  `status_vY'"

local status_offd = cond(`pass_offd', "{res}PASS", "{err}FAIL")
display as txt "  min(pi_12, pi_21)" _col(35) %8.4f `off_diag_min' ///
    "  threshold > 0.05  `status_offd'"
display as txt "    pi_12 hat = " %6.4f `pi_12_hat' ///
    "  pi_21 hat = " %6.4f `pi_21_hat'

local status_rho = cond(`pass_rho', "{res}PASS", "{err}FAIL")
display as txt "  |rho| (truth ~ 0.17)" _col(35) %8.4f abs(`rho_hat') ///
    "  band [0.05, 0.45]  `status_rho'"

local status_conv = cond(`pass_conv', "{res}PASS", "{err}FAIL")
display _n as txt "  Strict convergence (e(converged)==1):  `status_conv'"
display as txt "    |grad| =" %9.2e `gn'
display as txt "    |grad|/(1+|LL|) =" %9.2e `rg'
display as txt "    V positive definite = `pd'"

* ============================================================================
* PART 3: Final verdict
* ============================================================================

local n_pass = `pass_conv' + `pass_delta' + `pass_v_T' + `pass_v_Y' + ///
    `pass_offd' + `pass_rho'

display _n as txt "{hline 70}"
if `n_pass' == 6 {
    display as res "BIVARIATE CERTIFICATION PASSED (6/6 criteria met)"
    display as txt "  hsmixture_bivariate recovers the bivariate-grid DGP"
    display as txt "  within tolerance, including off-diagonal probability mass."
    display as txt "{hline 70}"
    capture log close _hsmix_biv
}
else {
    display as err "BIVARIATE CERTIFICATION FAILED (`n_pass'/6 criteria met)"
    display as err "  Inspect the model output above for details."
    display as err "  Possible causes: spike-and-slab corner solution, BFGS"
    display as err "  hit a local mode (re-run with from()), or genuine"
    display as err "  identification issue at this DGP scale."
    display as txt "{hline 70}"
    capture log close _hsmix_biv
    error 9
}
