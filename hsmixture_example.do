*! hsmixture_example.do
*! SMOKE TEST and example for the hsmixture package.
*!
*! This script:
*!   1. Generates a synthetic person-period dataset with known DGP
*!   2. Estimates all three models (single, joint, bivariate)
*!   3. Runs postestimation diagnostics
*!   4. Demonstrates model comparison across K
*!
*! What this script PROVES: the three commands run end-to-end without
*! errors and produce ereturn outputs that the postestimation tools can
*! consume. It is NOT a parameter-recovery certification: the small-N /
*! low-event-rate DGP used here is at the edge of identifiability for a
*! K=2 mixture, and the joint and bivariate fits routinely terminate at
*! a spike-and-slab corner with non-PD V. That is a feature of the data
*! at this scale, not a bug in the estimator.
*!
*! For parameter-recovery certification with stricter tolerances and a
*! stronger DGP, see hsmixture_certification.do (separate file).
*!
*! Authors: Jonghoon Park and R. Alan Seals (Auburn University)

clear all
set seed 20260321
set more off

* Log output for debugging -- writes to Stata's current working directory
capture log close _hsmixture_log
display as txt "Working directory: `c(pwd)'"
log using "hsmixture_example_log.txt", text replace name(_hsmixture_log)

* ============================================================================
* PART 1: Generate synthetic person-period data
* ============================================================================

* DGP: Two latent types with correlated heterogeneity.
*   - Treatment is a time-varying event (e.g., program entry, job training)
*   - Outcome is a competing event (e.g., job exit, recidivism)
*   - Treatment increases outcome hazard (delta > 0)
*
* True parameters (must match the locals defined below):
*   delta     = 0.5  (treatment raises outcome hazard by ~65%)
*   sigma_T   = 0.5  (heterogeneity loads into treatment eq)
*   sigma_Y   = 0.8  (heterogeneity loads into outcome eq)
*   pi_1      = 0.6  (60% are type 1, v=0)
*   pi_2      = 0.4  (40% are type 2, v=1.0)
*   beta_x1   = 0.3  (covariate effect on both hazards)

display _n as txt "{hline 70}"
display as txt "hsmixture package: Example and smoke-test script"
display as txt "{hline 70}"
display as txt "Generating synthetic person-period data..."

local N_persons = 1000
local T_max     = 20

* True DGP parameters (moderate magnitudes for reliable convergence)
local delta_true   = 0.5
local sigma_T_true = 0.5
local sigma_Y_true = 0.8
local pi_1_true    = 0.6
local v_2_true     = 1.0
local beta_x1_true = 0.3

* Create person-level data first
quietly {
    set obs `N_persons'
    gen long id = _n

    * Covariate (standard normal)
    gen double x1 = rnormal()

    * Assign latent type
    gen byte type = runiform() > `pi_1_true'
    gen double v = type * `v_2_true'

    * Expand to person-period
    expand `T_max'
    bysort id: gen int period = _n

    * ----------------------------------------------------------------
    * Step 1: Generate treatment timing
    *   Treatment is NOT an absorbing state -- observation continues.
    *   Baseline hazard is low (~2% per period for type 1).
    * ----------------------------------------------------------------
    gen double eta_T = -4.5 + 0.04 * period + `beta_x1_true' * x1 ///
        + `sigma_T_true' * v
    gen double h_T = 1 - exp(-exp(eta_T))

    * Independent uniform draw per person-period
    gen double u_T = runiform()
    gen byte treat_event = 0
    gen treat_time = .

    * Determine treatment timing: first period where u < h, reading forward
    forvalues t = 1/`T_max' {
        bysort id (period): replace treat_event = (u_T < h_T) ///
            if period == `t' & treat_time == .
        bysort id (period): replace treat_time = `t' ///
            if treat_event == 1 & period == `t' & treat_time == .
    }

    * Fill treatment time across all periods for this person
    bysort id: egen _tt = min(treat_time)
    replace treat_time = _tt
    drop _tt

    * Treatment indicator: 1 in all periods AFTER treatment onset
    * (no-anticipation: treatment does not affect outcome in the same period)
    gen byte treated = (period > treat_time) if treat_time != .
    replace treated = 0 if treat_time == .

    * Clean treat_event to only fire at the treatment period
    replace treat_event = 0
    replace treat_event = 1 if period == treat_time

    * ----------------------------------------------------------------
    * Step 2: Generate outcome timing
    *   Outcome IS the absorbing state -- censor here.
    *   Hazard depends on treated (which switches on after treatment).
    * ----------------------------------------------------------------
    gen double eta_Y = -4.0 + 0.03 * period + `beta_x1_true' * x1 ///
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

    * ----------------------------------------------------------------
    * Step 3: Censor at OUTCOME only (outcome is absorbing)
    *   Treatment is NOT absorbing -- we keep post-treatment periods.
    *   Treatment risk set: at risk for treatment only until treated.
    * ----------------------------------------------------------------
    gen int censor_time = cond(outcome_time != ., outcome_time, `T_max')
    keep if period <= censor_time

    * Treatment risk set indicator: 1 if not yet treated (at risk for treatment)
    * After treatment, treat_event = 0 and individual is no longer at risk
    * for the treatment equation (but still observed for outcome equation)
    gen byte treat_at_risk = (treated == 0)

    * Clean up auxiliary variables
    drop eta_T h_T u_T eta_Y h_Y u_Y treat_time outcome_time censor_time type v

    * Label variables
    label variable id "Person identifier"
    label variable period "Time period"
    label variable x1 "Covariate"
    label variable treat_event "Treatment event indicator"
    label variable outcome_event "Outcome event indicator"
    label variable treated "Post-treatment indicator"
    label variable treat_at_risk "At risk for treatment (1 = not yet treated)"

    * Create period dummies manually (avoids fv version issues)
    tab period, gen(pd_)
    drop pd_1  // base category
}

display as txt "Data generated: " _continue
describe, short
display _n as txt "True DGP: delta = `delta_true', sigma_T = `sigma_T_true', " ///
    "sigma_Y = `sigma_Y_true'"
display as txt "          pi_1 = `pi_1_true', v_2 = `v_2_true', " ///
    "beta_x1 = `beta_x1_true'"

* Save for reuse (both tempfile and permanent for help file examples)
tempfile example_data
save `example_data'
save hsmixture_example_data, replace


* ============================================================================
* PART 2: Single-equation model (hsmixture)
* ============================================================================

display _n as txt "{hline 70}"
display as txt "EXAMPLE 1: Single-equation hazard with Heckman-Singer heterogeneity"
display as txt "{hline 70}"

* Diagnostic: which file is Stata loading?
which hsmixture
adopath
display "Stata version: `c(stata_version)'"

hsmixture outcome_event x1 pd_*, id(id) k(2)


* ============================================================================
* PART 3: Joint timing-of-events model (hsmixture_joint)
* ============================================================================

display _n as txt "{hline 70}"
display as txt "EXAMPLE 2: Joint Timing-of-Events (K=2)"
display as txt "{hline 70}"

hsmixture_joint (treat_event = x1 pd_*) ///
    (outcome_event = x1 pd_*, treat(treated)) ///
    , id(id) k(2) iterate(300) riskset(treat_at_risk)

estimates store m_joint_k2

* Postestimation diagnostics
hsmixture_joint_postestimation, all


* ============================================================================
* PART 4: Joint model with K=3
* ============================================================================

display _n as txt "{hline 70}"
display as txt "EXAMPLE 3: Joint Timing-of-Events (K=3)"
display as txt "{hline 70}"

hsmixture_joint (treat_event = x1 pd_*) ///
    (outcome_event = x1 pd_*, treat(treated)) ///
    , id(id) k(3) iterate(300) riskset(treat_at_risk)

estimates store m_joint_k3


* ============================================================================
* PART 5: Model comparison
* ============================================================================

display _n as txt "{hline 70}"
display as txt "EXAMPLE 4: Model comparison (K=2 vs K=3)"
display as txt "{hline 70}"

* NOTE: estimates stats computes BIC from the row count e(N) (person-periods).
* This package reports BIC on the person count e(N_persons), which can rank K
* differently because the larger row count inflates the per-parameter penalty.
* For the person-count BIC, run hsmixture_joint_postestimation, compare after
* each fit (it reports the person-count BIC for the active model).
estimates stats m_joint_k2 m_joint_k3

* LR test
hsmixture_joint_postestimation, lrtest(m_joint_k2)


* ============================================================================
* PART 6: Bivariate heterogeneity model
* ============================================================================

display _n as txt "{hline 70}"
display as txt "EXAMPLE 5: Bivariate Heterogeneity (2x2 grid)"
display as txt "{hline 70}"

use `example_data', clear

hsmixture_bivariate (treat_event = x1 pd_*) ///
    (outcome_event = x1 pd_*, treat(treated)) ///
    , id(id) nstarts(3) iterate(300) riskset(treat_at_risk)

estimates store m_bivariate

* Diagnostics
hsmixture_joint_postestimation, all


* ============================================================================
* PART 7: Summary
* ============================================================================

display _n as txt "{hline 70}"
display as txt "SUMMARY"
display as txt "{hline 70}"
display as txt "True DGP: delta=0.5, lambda_T=0.5, lambda_Y=0.8, pi_1=0.6, v_2=1.0"
display as txt ""

* NOTE: estimates stats uses row-count BIC (e(N)); the package reports
* person-count BIC (e(N_persons)). See the note in PART 5. For the person-count
* BIC per model, run hsmixture_joint_postestimation, compare after each fit.
estimates stats m_joint_k2 m_joint_k3 m_bivariate

* Honest convergence summary. Each estimator exposes a strict convergence
* flag (e(converged) requires BFGS converged AND |grad| small AND V positive
* definite). With this small synthetic dataset (1000 persons; the realized
* event counts are roughly 213 treatment events and 542 outcome events) the
* joint and bivariate fits typically hit the spike-and-slab corner and fail
* strict convergence. The single-equation hsmixture often does too. This is
* a property of the DGP at this sample size, not a defect in the estimator.
display _n as txt "Strict convergence by model:"
local n_strict 0
local n_total 0
foreach m in m_joint_k2 m_joint_k3 m_bivariate {
    quietly estimates restore `m'
    local conv = e(converged)
    local gn   = e(grad_norm)
    local pd   = e(v_pd)
    local ++n_total
    if `conv' == 1 {
        display as res "  `m': strict convergence" ///
            " (gradient norm = " %9.2e `gn' ", V positive definite)"
        local ++n_strict
    }
    else {
        display as err "  `m': DID NOT strictly converge" ///
            " (gradient norm = " %9.2e `gn' ", V positive definite = `pd')"
    }
}

display _n as txt "{hline 70}"
display as txt "SMOKE TEST RESULT"
display as txt "{hline 70}"
display as txt "  Commands executed without runtime errors:" _col(50) as res "PASS"
display as txt "  Strict convergence achieved on:" _col(50) ///
    as res "`n_strict' / `n_total' models"
if `n_strict' < `n_total' {
    display _n as err "  This is NOT a certification of parameter recovery."
    display as err "  The DGP in this script is at the edge of identifiability" ///
        " for K=2."
    display as err "  Run hsmixture_certification.do for the recovery test on a" ///
        " stronger DGP."
}
display as txt "  Always check e(converged), e(grad_norm), and the postestimation"
display as txt "  diagnostics before interpreting any HR or CI from a fitted model."
display as txt "{hline 70}"

capture log close _hsmixture_log
