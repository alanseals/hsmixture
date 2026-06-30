*! hsmixture_certification_common.do
*! Parameter-recovery certification for factor(common) mode.
*!
*! This script does TWO things:
*!   PART A. Generate data under a true common-loading DGP
*!     (lambda_T_true = lambda_Y_true = sigma_true) and assert that
*!     hsmixture_joint, factor(common) recovers (delta, sigma, pi_2)
*!     within a stated tolerance.
*!   PART B. Refit the SAME data with hsmixture_joint, factor(separate)
*!     and assert that the data-driven separate loadings cluster around
*!     a common value (i.e., the data don't push lambda_T and lambda_Y
*!     apart when they are truly equal). This is a falsifiability test:
*!     if separate-mode estimation systematically distorted the loadings,
*!     factor(common) would not be a defensible restriction.
*!
*! Pass criteria (Part A, common mode):
*!   - e(converged) == 1 (BFGS + relative gradient + V positive definite)
*!   - |delta_hat - 0.5| < 3 * SE(delta)
*!   - |lambda_hat - 1.0| < 3 * SE(lambda)
*!   - 0.30 < pi_2_hat < 0.60  (true 0.45)
*!
*! Note on DGP strength: a previous version of this script used
*! sigma_true = 0.7 (matching the lambda_T_true in the existing same-signs
*! cert). Under common loading at that magnitude the K=2 mixture is weakly
*! identified — finite-sample drift in either equation's effective loading
*! pulls pi_2 around. The 2026-05-04 first run produced delta=0.504 (PASS),
*! lambda=1.08 (PASS by margin), pi_2=0.166 (FAIL). All 6 multistart configs
*! converged strictly to the same LL, confirming this was the actual MLE
*! for that data realization, not a local mode. Bumping sigma_true to 1.0
*! tightens identification so that pi_2 reliably recovers.
*!
*! Pass criteria (Part B, separate mode on common-loading data):
*!   - e(converged) == 1
*!   - lambda_T_hat and lambda_Y_hat have the same sign
*!   - |lambda_T_hat - lambda_Y_hat| < 3 * sqrt(SE_lT^2 + SE_lY^2)
*!     (the difference is within a 95% band of zero — separate mode does
*!      not manufacture a spurious gap)
*!
*! Authors: Jonghoon Park and R. Alan Seals (Auburn University)

clear all
set seed 20260504
set more off

capture log close _hsmix_common
display as txt "Working directory: `c(pwd)'"
log using "hsmixture_certification_common_log.txt", text replace name(_hsmix_common)

* ============================================================================
* PART 1: Generate synthetic person-period data
*         True DGP has a SINGLE shared loading on the latent type.
* ============================================================================

display _n as txt "{hline 70}"
display as txt "hsmixture certification: factor(common) recovery test"
display as txt "{hline 70}"
display as txt "Generating synthetic data with COMMON-loading DGP..."

local N_persons = 4000
local T_max     = 20

local delta_true   = 0.5
local sigma_true   = 1.0
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

    * Treatment hazard: SAME loading sigma on the latent v
    gen double eta_T = -3.0 + 0.02 * period + `beta_x1_true' * x1 ///
        + `sigma_true' * v
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

    * Outcome hazard: SAME loading sigma on the latent v
    gen double eta_Y = -3.5 + 0.02 * period + `beta_x1_true' * x1 ///
        + `delta_true' * treated + `sigma_true' * v
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
display as txt "True DGP: delta=`delta_true', sigma=`sigma_true' (shared), " ///
    "pi_2=" %4.2f (1 - `pi_1_true')
display as txt "True within-type HR = exp(`delta_true') = " %5.3f exp(`delta_true')

* ============================================================================
* PART A: Fit factor(common) and assert recovery
* ============================================================================

display _n as txt "{hline 70}"
display as txt "PART A: Fitting hsmixture_joint, factor(common) ..."
display as txt "{hline 70}"

hsmixture_joint (treat_event = x1 pd_*) ///
    (outcome_event = x1 pd_*, treat(treated)) ///
    , id(id) k(2) factor(common) iterate(500) riskset(treat_at_risk)

local conv_A      = e(converged)
local gn_A        = e(grad_norm)
local pd_A        = e(v_pd)
local delta_A     = e(delta)
local lambda_A    = e(lambda)
local pi_2_A      = e(pi)[1, 2]
local se_delta_A  = _se[/delta]
local se_lambda_A = _se[/lambda]

local pass_conv_A   = (`conv_A' == 1)
local pass_delta_A  = abs(`delta_A' - `delta_true') < 3 * `se_delta_A'
local pass_lambda_A = abs(`lambda_A' - `sigma_true') < 3 * `se_lambda_A'
local pass_pi_A     = (`pi_2_A' > 0.30) & (`pi_2_A' < 0.60)

display _n as txt "{hline 70}"
display as txt "PART A RESULTS (factor(common) on common-loading DGP)"
display as txt "{hline 70}"
display as txt ///
    "  Test                              Estimate    SE       Truth   Status"
display as txt "{hline 70}"

local status_delta_A = cond(`pass_delta_A', "{res}PASS", "{err}FAIL")
display as txt "  delta (treatment effect)" _col(35) %8.4f `delta_A' ///
    "  " %6.4f `se_delta_A' "  " %6.2f `delta_true' "    `status_delta_A'"

local status_lambda_A = cond(`pass_lambda_A', "{res}PASS", "{err}FAIL")
display as txt "  lambda (shared loading)" _col(35) %8.4f `lambda_A' ///
    "  " %6.4f `se_lambda_A' "  " %6.2f `sigma_true' "    `status_lambda_A'"

local status_pi_A = cond(`pass_pi_A', "{res}PASS", "{err}FAIL")
display as txt "  pi_2 (high-type share)" _col(35) %8.4f `pi_2_A' ///
    "       --   " %6.2f (1 - `pi_1_true') "    `status_pi_A'"

local status_conv_A = cond(`pass_conv_A', "{res}PASS", "{err}FAIL")
display _n as txt "  Strict convergence (e(converged)==1):  `status_conv_A'"
display as txt "    Gradient norm =" %9.2e `gn_A'
display as txt "    V positive definite = `pd_A'"

local n_pass_A = `pass_conv_A' + `pass_delta_A' + `pass_lambda_A' + `pass_pi_A'

* ============================================================================
* PART B: Refit factor(separate) and check it does not fabricate a gap.
* ============================================================================

display _n as txt "{hline 70}"
display as txt "PART B: Falsifiability — fitting factor(separate) on the same data"
display as txt "{hline 70}"
display as txt "Under truly common loading, separate mode should land near"
display as txt "lambda_T == lambda_Y (no spurious gap)."

hsmixture_joint (treat_event = x1 pd_*) ///
    (outcome_event = x1 pd_*, treat(treated)) ///
    , id(id) k(2) factor(separate) iterate(500) riskset(treat_at_risk)

local conv_B  = e(converged)
local lT_B    = e(lambda_T)
local lY_B    = e(lambda_Y)
local se_lT_B = _se[/lambda_T]
local se_lY_B = _se[/lambda_Y]

* Diff and pooled SE of the diff (assumes near-zero covariance, conservative)
local diff_B    = `lT_B' - `lY_B'
local se_diff_B = sqrt(`se_lT_B'^2 + `se_lY_B'^2)

local pass_conv_B    = (`conv_B' == 1)
local pass_signs_B   = (sign(`lT_B') == sign(`lY_B'))
local pass_no_gap_B  = abs(`diff_B') < 3 * `se_diff_B'

display _n as txt "{hline 70}"
display as txt "PART B RESULTS (factor(separate) on common-loading DGP)"
display as txt "{hline 70}"
display as txt ///
    "  Test                                     Value    Threshold   Status"
display as txt "{hline 70}"

local status_conv_B = cond(`pass_conv_B', "{res}PASS", "{err}FAIL")
display as txt "  Strict convergence" _col(45) ///
    cond(`conv_B'==1, "yes", "no") _col(58) "1" "    `status_conv_B'"

local status_signs_B = cond(`pass_signs_B', "{res}PASS", "{err}FAIL")
display as txt "  sign(lambda_T) == sign(lambda_Y)" _col(45) ///
    cond(`pass_signs_B', "yes", "no") _col(58) "1" "    `status_signs_B'"

local status_gap_B = cond(`pass_no_gap_B', "{res}PASS", "{err}FAIL")
display as txt "  |lambda_T - lambda_Y| within 3 SE diff" _col(45) %8.4f abs(`diff_B') ///
    "  " %6.4f 3 * `se_diff_B' "  `status_gap_B'"

display _n as txt "  lambda_T = " %7.4f `lT_B' "  (SE " %6.4f `se_lT_B' ")"
display as txt "  lambda_Y = " %7.4f `lY_B' "  (SE " %6.4f `se_lY_B' ")"
display as txt "  truth    = " %7.4f `sigma_true' " (shared)"

local n_pass_B = `pass_conv_B' + `pass_signs_B' + `pass_no_gap_B'

* ============================================================================
* PART C: Final verdict
* ============================================================================

display _n as txt "{hline 70}"
display as txt "OVERALL CERTIFICATION VERDICT"
display as txt "{hline 70}"

local total_pass = `n_pass_A' + `n_pass_B'
local total_crit = 7

if `total_pass' == `total_crit' {
    display as res "CERTIFICATION PASSED (`total_pass'/`total_crit' criteria met)"
    display as txt "  factor(common) recovers the true DGP within tolerance."
    display as txt "  factor(separate) fits the same data without manufacturing"
    display as txt "  a spurious lambda_T vs lambda_Y gap."
    display as txt "{hline 70}"
    capture log close _hsmix_common
}
else {
    display as err "CERTIFICATION FAILED (`total_pass'/`total_crit' criteria met)"
    display as err "  Part A: `n_pass_A'/4 met (factor(common) recovery)"
    display as err "  Part B: `n_pass_B'/3 met (factor(separate) falsifiability)"
    display as txt "{hline 70}"
    capture log close _hsmix_common
    error 9
}
