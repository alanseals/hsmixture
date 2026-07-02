*! hsmixture_certification.do
*! Parameter-recovery certification for hsmixture_joint (joint timing-of-events model).
*!
*! Unlike hsmixture_example.do (which is a smoke test on a deliberately
*! small DGP), this script targets a sample size and effect size where
*! the joint K=2 Heckman-Singer model is well identified, and asserts
*! that the package recovers the true (delta, lambda_T, lambda_Y, pi_2)
*! within a stated tolerance. The script exits non-zero if the
*! certification criteria fail.
*!
*! DGP rationale:
*!   - N_persons = 4000  (vs 1000 in the smoke test): roughly 4x more
*!     events, so the K=2 mixture is identified rather than collapsing.
*!   - Treatment baseline hazard ~5% per period (vs ~1% in smoke test):
*!     ensures enough treated person-periods that delta is identified
*!     in both types.
*!   - Outcome baseline hazard ~3% per period (vs ~2%): ensures enough
*!     outcome events per type.
*!   - sigma_T = 0.7, sigma_Y = 1.0 (vs 0.5, 0.8): high-type hazard
*!     ratio is now exp(0.7)=2.0x for treatment and exp(1.0)=2.7x for
*!     outcome. Type separation is unambiguous.
*!   - pi_1 = 0.55, pi_2 = 0.45 (vs 0.6, 0.4): more balanced mixture.
*!   - delta = 0.5 (unchanged): treatment raises outcome hazard ~65%.
*!   - beta_x1 = 0.3 (unchanged).
*!
*! Pass criteria (for the joint K=2 model):
*!   - e(converged) == 1 (BFGS converged AND |grad|/(1+|LL|) < 1e-5 AND V is PD)
*!   - |delta_hat - 0.5| < 3 * SE(delta)   (95%-band recovery)
*!   - |lambda_T_hat - 0.7| < 3 * SE(lambda_T)
*!   - |lambda_Y_hat - 1.0| < 3 * SE(lambda_Y)
*!   - 0.30 < pi_2_hat < 0.60  (true 0.45)
*!
*! Authors: Jonghoon Park and R. Alan Seals (Auburn University)

clear all
set seed 20260501
set more off

capture log close _hsmixture_cert
display as txt "Working directory: `c(pwd)'"
log using "hsmixture_certification_log.txt", text replace name(_hsmixture_cert)

* ============================================================================
* PART 1: Generate synthetic person-period data (strong DGP)
* ============================================================================

display _n as txt "{hline 70}"
display as txt "hsmixture certification: parameter-recovery test"
display as txt "{hline 70}"
display as txt "Generating synthetic data with strong-identification DGP..."

local N_persons = 4000
local T_max     = 20

local delta_true   = 0.5
local sigma_T_true = 0.7
local sigma_Y_true = 1.0
local pi_1_true    = 0.55
local v_2_true     = 1.0
local beta_x1_true = 0.3

quietly {
    set obs `N_persons'
    gen long id = _n

    gen double x1 = rnormal()
    gen byte type = runiform() > `pi_1_true'
    gen double v = type * `v_2_true'

    expand `T_max'
    bysort id: gen int period = _n

    * Treatment hazard: baseline -3.0 (~4.7% per period for type 1)
    gen double eta_T = -3.0 + 0.02 * period + `beta_x1_true' * x1 ///
        + `sigma_T_true' * v
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

    * Outcome hazard: baseline -3.5 (~3.0% per period for untreated type 1)
    gen double eta_Y = -3.5 + 0.02 * period + `beta_x1_true' * x1 ///
        + `delta_true' * treated + `sigma_Y_true' * v
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

    drop eta_T h_T u_T eta_Y h_Y u_Y treat_time outcome_time censor_time type v

    label variable id "Person identifier"
    label variable period "Time period"
    label variable x1 "Covariate"
    label variable treat_event "Treatment event indicator"
    label variable outcome_event "Outcome event indicator"
    label variable treated "Post-treatment indicator"
    label variable treat_at_risk "At risk for treatment (1 = not yet treated)"

    tab period, gen(pd_)
    drop pd_1
}

* Report event counts. The certification DGP should produce far more events
* per type than the smoke-test DGP.
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
display as txt "True DGP: delta=`delta_true', sigma_T=`sigma_T_true', " ///
    "sigma_Y=`sigma_Y_true', pi_2=" %4.2f (1 - `pi_1_true')

* ============================================================================
* PART 2: Fit joint K=2 model and assert recovery
* ============================================================================

display _n as txt "{hline 70}"
display as txt "Fitting joint K=2 Heckman-Singer model..."
display as txt "{hline 70}"

hsmixture_joint (treat_event = x1 pd_*) ///
    (outcome_event = x1 pd_*, treat(treated)) ///
    , id(id) k(2) iterate(500) riskset(treat_at_risk)

* Read estimation results and the strict-convergence flag.
local conv      = e(converged)
local gn        = e(grad_norm)
local pd        = e(v_pd)
local delta_hat = e(delta)
local lT_hat    = e(lambda_T)
local lY_hat    = e(lambda_Y)
local pi_2_hat  = e(pi)[1, 2]
local se_delta  = _se[/delta]
local se_lT     = _se[/lambda_T]
local se_lY     = _se[/lambda_Y]

* Tolerance: 3 SEs for point estimates. pi_2 has a hard band [0.30, 0.60]
* because its SE under the softmax parameterization is awkward to compute
* on the same scale.
local pass_conv  = (`conv' == 1)
local pass_delta = abs(`delta_hat' - 0.5) < 3 * `se_delta'
local pass_lT    = abs(`lT_hat' - `sigma_T_true') < 3 * `se_lT'
local pass_lY    = abs(`lY_hat' - `sigma_Y_true') < 3 * `se_lY'
local pass_pi    = (`pi_2_hat' > 0.30) & (`pi_2_hat' < 0.60)

display _n as txt "{hline 70}"
display as txt "CERTIFICATION RESULTS (joint K=2)"
display as txt "{hline 70}"
display as txt ///
    "  Test                              Estimate    SE       Truth   Status"
display as txt "{hline 70}"

* delta
local status_delta = cond(`pass_delta', "{res}PASS", "{err}FAIL")
display as txt "  delta (treatment effect)" _col(35) %8.4f `delta_hat' ///
    "  " %6.4f `se_delta' "  " %6.2f 0.5 "    `status_delta'"

* lambda_T
local status_lT = cond(`pass_lT', "{res}PASS", "{err}FAIL")
display as txt "  lambda_T (treat loading)" _col(35) %8.4f `lT_hat' ///
    "  " %6.4f `se_lT' "  " %6.2f `sigma_T_true' "    `status_lT'"

* lambda_Y
local status_lY = cond(`pass_lY', "{res}PASS", "{err}FAIL")
display as txt "  lambda_Y (outcome loading)" _col(35) %8.4f `lY_hat' ///
    "  " %6.4f `se_lY' "  " %6.2f `sigma_Y_true' "    `status_lY'"

* pi_2
local status_pi = cond(`pass_pi', "{res}PASS", "{err}FAIL")
display as txt "  pi_2 (high-type share)" _col(35) %8.4f `pi_2_hat' ///
    "       --   " %6.2f (1 - `pi_1_true') "    `status_pi'"

* Strict convergence
local status_conv = cond(`pass_conv', "{res}PASS", "{err}FAIL")
display _n as txt "  Strict convergence (e(converged)==1):  `status_conv'"
display as txt "    Gradient norm =" %9.2e `gn'
display as txt "    V positive definite = `pd'"

* ============================================================================
* PART 3: Final verdict
* ============================================================================

local n_pass = `pass_conv' + `pass_delta' + `pass_lT' + `pass_lY' + `pass_pi'

display _n as txt "{hline 70}"
if `n_pass' == 5 {
    display as res "CERTIFICATION PASSED (5/5 criteria met)"
    display as txt "  hsmixture_joint recovers the true DGP parameters within tolerance."
    display as txt "{hline 70}"
    capture log close _hsmixture_cert
}
else {
    display as err "CERTIFICATION FAILED (`n_pass'/5 criteria met)"
    display as err "  Inspect the model output above for details."
    display as err "  Possible causes: BFGS hit a local mode (re-run with from()),"
    display as err "  rank-deficient Hessian (try smaller K), or genuine DGP issue."
    display as txt "{hline 70}"
    capture log close _hsmixture_cert
    error 9
}
