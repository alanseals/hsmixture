*! hsmixture_certification_opposite_signs.do
*! Parameter-recovery certification for the OPPOSITE-SIGNS heterogeneity case.
*!
*! The standard hsmixture_certification.do tests recovery when lambda_T and
*! lambda_Y have the same sign (both positive). In applied joint timing-of-
*! events models the two loadings can take OPPOSITE signs -- a latent type
*! more prone to the treatment event can be less prone to the outcome event --
*! which can push the joint-ToE hazard ratio well above the simple cloglog HR.
*!
*! The standard cert never validates this region. This script does.
*!
*! DGP rationale:
*!   - lambda_T_true = +1.0  (type 2 more treatment-prone)
*!   - lambda_Y_true = -1.5  (type 2 LESS outcome-prone)
*!   - pi_2 = 0.45, delta = 1.0 (HR_within_type = exp(1) = 2.72)
*!
*! Predicted simple cloglog HR (from DGP arithmetic):
*!   - Treated avg outcome rate = 0.7 * exp(-1.5) * exp(1) + 0.3 * exp(1)
*!                              ≈ 0.43 + 0.82 = 1.25 * baseline
*!   - Untreated avg outcome rate = 0.3 * exp(-1.5) + 0.7
*!                                ≈ 0.07 + 0.70 = 0.77 * baseline
*!   - Cloglog HR ≈ 1.25 / 0.77 ≈ 1.62
*!   - True within-type HR = exp(1.0) = 2.72
*!   - So cloglog UNDERESTIMATES by ~40%, joint-ToE should correct upward.
*!
*! Pass criteria (joint K=2):
*!   - e(converged) == 1
*!   - |delta_hat - 1.0| < 3 * SE(delta)
*!   - sign(lambda_T) != sign(lambda_Y)  (opposite-signs structure preserved
*!     under labeling swap)
*!   - ||lambda_T| - 1.0| < 3 * SE(lambda_T)  (magnitude recovered)
*!   - ||lambda_Y| - 1.5| < 3 * SE(lambda_Y)  (magnitude recovered)
*!   - min(pi_1, pi_2) in [0.30, 0.60]  (mixture share recovered up to swap)
*!   - Cloglog HR (sanity) noticeably smaller than joint-ToE HR
*!
*! Note: lambda_T and lambda_Y signs are NOT individually identified -- the
*! likelihood is invariant under (pi_1 <-> pi_2, lambda_T -> -lambda_T,
*! lambda_Y -> -lambda_Y) plus baseline-constant adjustments. Both labelings
*! are valid MLEs. The cert tests what IS identified: opposite-signs
*! structure, magnitudes, delta, and the unordered mixture share.
*!
*! Authors: Jonghoon Park and R. Alan Seals (Auburn University)

clear all
set seed 20260503
set more off

capture log close _hsmix_opp
display as txt "Working directory: `c(pwd)'"
log using "hsmixture_certification_opposite_signs_log.txt", text replace name(_hsmix_opp)

* Provenance record: what code this certification actually ran on. Copy this
* block's output (plus the package git commit SHA) into CERTIFICATION.md
* when recording a release-gating run.
display _n as txt "{hline 70}"
display as txt "PROVENANCE"
display as txt "  Run date/time:  `c(current_date)' `c(current_time)'"
display as txt "  Stata version:  `c(stata_version)' (born `c(born_date)')"
display as txt "  OS / machine:   `c(os)' / `c(machine_type)'"
display as txt "  MP / cores:     `c(MP)' / `c(processors)'"
display as txt "  Command under test:"
which hsmixture_joint
display as txt "{hline 70}"

* ============================================================================
* PART 1: Generate synthetic person-period data with OPPOSITE-SIGNS DGP
* ============================================================================

display _n as txt "{hline 70}"
display as txt "hsmixture certification: OPPOSITE-SIGNS recovery test"
display as txt "{hline 70}"
display as txt "Generating synthetic data with negative-correlation heterogeneity..."

local N_persons = 6000
local T_max     = 20

local delta_true    = 1.0
local lambda_T_true = 1.0
local lambda_Y_true = -1.5
local pi_1_true     = 0.55
local v_2_true      = 1.0
local beta_x1_true  = 0.3

quietly {
    set obs `N_persons'
    gen long id = _n

    gen double x1 = rnormal()
    gen byte type = runiform() > `pi_1_true'
    gen double v = type * `v_2_true'

    expand `T_max'
    bysort id: gen int period = _n

    * Treatment hazard
    gen double eta_T = -3.0 + 0.02 * period + `beta_x1_true' * x1 ///
        + `lambda_T_true' * v
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

    * Outcome hazard with NEGATIVE lambda_Y -- type 2 has LOWER baseline outcome
    gen double eta_Y = -3.0 + 0.02 * period + `beta_x1_true' * x1 ///
        + `delta_true' * treated + `lambda_Y_true' * v
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
display as txt "True DGP: delta=`delta_true', lambda_T=`lambda_T_true', " ///
    "lambda_Y=`lambda_Y_true', pi_2=" %4.2f (1 - `pi_1_true')
display as txt "True within-type HR = exp(`delta_true') = " %5.3f exp(`delta_true')

* ============================================================================
* PART 2: Sanity check -- simple cloglog, expected to UNDERESTIMATE
* ============================================================================

display _n as txt "{hline 70}"
display as txt "Sanity check: simple cloglog (should be biased toward 1)..."
display as txt "{hline 70}"

quietly glm outcome_event treated x1 pd_*, ///
    family(binomial) link(cloglog) nolog
local hr_cloglog = exp(_b[treated])
local se_cloglog = _se[treated]
local z = invnormal(0.975)
local hr_lo_cl = exp(_b[treated] - `z' * `se_cloglog')
local hr_hi_cl = exp(_b[treated] + `z' * `se_cloglog')

display as txt "  Simple cloglog HR (biased)" _col(35) %8.4f `hr_cloglog' ///
    "  CI [" %5.3f `hr_lo_cl' ", " %5.3f `hr_hi_cl' "]"
display as txt "  True within-type HR" _col(35) %8.4f exp(`delta_true')

* ============================================================================
* PART 3: Fit joint K=2 model and assert recovery
* ============================================================================

display _n as txt "{hline 70}"
display as txt "Fitting joint K=2 Heckman-Singer model..."
display as txt "{hline 70}"

hsmixture_joint (treat_event = x1 pd_*) ///
    (outcome_event = x1 pd_*, treat(treated)) ///
    , id(id) k(2) iterate(500) riskset(treat_at_risk)

local conv      = e(converged)
local gn        = e(grad_norm)
local rg        = e(rel_grad)
local pd        = e(v_pd)
local delta_hat = e(delta)
local lT_hat    = e(lambda_T)
local lY_hat    = e(lambda_Y)
local pi_2_hat  = e(pi)[1, 2]
local hr_joint  = e(hr)
local se_delta  = _se[/delta]
local se_lT     = _se[/lambda_T]
local se_lY     = _se[/lambda_Y]

* Pass criteria. Signs of lambda individually are NOT identified (labeling
* swap symmetry); the cert tests what IS identified.
local pass_conv  = (`conv' == 1)
local pass_delta = abs(`delta_hat' - `delta_true') < 3 * `se_delta'
* Opposite-signs structure preserved under swap.
local pass_signs = (sign(`lT_hat') != sign(`lY_hat')) & ///
                   (sign(`lT_hat') != 0) & (sign(`lY_hat') != 0)
* Magnitudes (truth = |lambda_T|=1, |lambda_Y|=1.5).
local pass_lT_mag = abs(abs(`lT_hat') - abs(`lambda_T_true')) < 3 * `se_lT'
local pass_lY_mag = abs(abs(`lY_hat') - abs(`lambda_Y_true')) < 3 * `se_lY'
* Mixture share: under swap, the smaller of (pi_1, pi_2) corresponds to
* min(true pi_1, true pi_2) = 0.45.
local pi_1_hat = 1 - `pi_2_hat'
local pi_min_hat = min(`pi_1_hat', `pi_2_hat')
local pass_pi   = (`pi_min_hat' > 0.30) & (`pi_min_hat' < 0.60)

* Sanity: joint-ToE HR should be larger than cloglog HR (selection correction)
local pass_correction = (`hr_joint' > `hr_cloglog')

display _n as txt "{hline 70}"
display as txt "CERTIFICATION RESULTS (joint K=2, opposite-signs DGP)"
display as txt "{hline 70}"

local status_delta = cond(`pass_delta', "{res}PASS", "{err}FAIL")
display as txt "  delta (treatment effect)" _col(35) %8.4f `delta_hat' ///
    "  SE " %6.4f `se_delta' "  truth " %5.2f `delta_true' "  `status_delta'"

local status_signs = cond(`pass_signs', "{res}PASS", "{err}FAIL")
display as txt "  lambda_T, lambda_Y opposite signs" _col(40) ///
    cond(`pass_signs', "yes", "no") "  `status_signs'"

local status_lT = cond(`pass_lT_mag', "{res}PASS", "{err}FAIL")
display as txt "  |lambda_T| (truth 1.00)" _col(35) %8.4f abs(`lT_hat') ///
    "  SE " %6.4f `se_lT' "  truth " %5.2f abs(`lambda_T_true') "  `status_lT'"

local status_lY = cond(`pass_lY_mag', "{res}PASS", "{err}FAIL")
display as txt "  |lambda_Y| (truth 1.50)" _col(35) %8.4f abs(`lY_hat') ///
    "  SE " %6.4f `se_lY' "  truth " %5.2f abs(`lambda_Y_true') "  `status_lY'"

local status_pi = cond(`pass_pi', "{res}PASS", "{err}FAIL")
display as txt "  min(pi_1, pi_2) (truth 0.45)" _col(35) %8.4f `pi_min_hat' ///
    "  `status_pi'"

local status_corr = cond(`pass_correction', "{res}PASS", "{err}FAIL")
display _n as txt "  Selection correction (joint > cloglog):"
display as txt "    Cloglog HR:" _col(35) %8.4f `hr_cloglog'
display as txt "    Joint-ToE HR:" _col(35) %8.4f `hr_joint'
display as txt "    Joint > cloglog: `status_corr'"

local status_conv = cond(`pass_conv', "{res}PASS", "{err}FAIL")
display _n as txt "  Strict convergence (e(converged)==1):  `status_conv'"
display as txt "    |grad| =" %9.2e `gn'
display as txt "    |grad|/(1+|LL|) =" %9.2e `rg'
display as txt "    V positive definite = `pd'"
* See hsmixture_certification.do for why these two are printed. This fit
* produced the v_pd=0 verdict on 14jul2026 that prompted the PD
* investigation. The cause was not the 1e-8 floor. The minimum eigenvalue
* measured 1.85e-04, four orders of magnitude clear of it. The optimizer had
* aborted and the .ado posted the I*1e-20 placeholder V, so v_pd=0 was
* correct. The abort retry (v2.4.0) fixed that symptom. The scale-relative
* PD test is a separate and independently correct change.
display as txt "    V min eigenvalue =" %9.2e e(v_mineig)
display as txt "    V is placeholder/scaffold = " e(v_scaffold)

* ============================================================================
* PART 3b: Misspecification record — factor(common) on opposite-signs DGP.
* ============================================================================
* This is NOT a pass/fail criterion. We deliberately fit a model that the
* truth does not satisfy, to document the empirical signature. When the user
* sees this signature on real data (small lambda magnitude paired with
* convergence warnings and a delta close to the simple cloglog), they have
* informative evidence that opposite-signs heterogeneity is present in the
* data and that factor(common) is misspecified.

display _n as txt "{hline 70}"
display as txt "PART 3b (misspecification record): factor(common) on the same data"
display as txt "{hline 70}"
display as txt "factor(common) cannot represent opposite-signs heterogeneity by"
display as txt "construction. We fit it here to document what happens. The"
display as txt "expected signature is a small or unstable lambda and a delta"
display as txt "close to the simple cloglog HR, possibly with a convergence"
display as txt "warning. This is the diagnostic users should look for when"
display as txt "deciding which factor mode their real data wants."

* nstarts(1) iterate(30) caps this deliberately-misspecified fit, which is an
* informational record and NOT a scored criterion. factor(common) cannot
* represent opposite-signs heterogeneity, so it never strictly converges; it
* stalls on a flat region where each BFGS iteration degenerates into repeated
* failed line searches plus a full numeric gradient. Measured cost there on a
* 4-core Mac (15jul2026): ~48 seconds per iteration, five times the pace of
* the scored fit above. At the previous budget (6 default starts x 200
* iterations) this unscored demo could run 16 hours, dwarfing the entire rest
* of the certification suite. One start and 30 iterations is enough: the
* diagnostic signature is that the fit does NOT converge and delta collapses
* toward the naive cloglog, and that is visible immediately. Extra starts and
* iterations only buy a longer log of the same frozen likelihood. (The
* 02jul2026 Windows run stopped on its own near iteration 55 via a
* discontinuous-region abort; that early stop is a machine-dependent accident,
* not something to rely on -- hence an explicit cap.)
capture noisily {
    hsmixture_joint (treat_event = x1 pd_*) ///
        (outcome_event = x1 pd_*, treat(treated)) ///
        , id(id) k(2) factor(common) nstarts(1) iterate(30) riskset(treat_at_risk)
}
local rc_common = _rc

if `rc_common' == 0 {
    local conv_c   = e(converged)
    local rg_c     = e(rel_grad)
    local pd_c     = e(v_pd)
    local delta_c  = e(delta)
    local lambda_c = e(lambda)
    local hr_c     = e(hr)
    local pi_2_c   = e(pi)[1, 2]

    display _n as txt "  factor(common) results on opposite-signs DGP:"
    display as txt "    delta:" _col(28) %9.4f `delta_c' ///
        "  (truth " %5.2f `delta_true' "; cloglog " %5.2f ln(`hr_cloglog') ")"
    display as txt "    HR (exp delta):" _col(28) %9.4f `hr_c' ///
        "  (truth " %5.2f exp(`delta_true') "; cloglog " %5.2f `hr_cloglog' ")"
    display as txt "    lambda (shared):" _col(28) %9.4f `lambda_c'
    display as txt "    pi_2:" _col(28) %9.4f `pi_2_c'
    display as txt "    Strict convergence:" _col(28) %9.0f `conv_c'
    display as txt "    |grad|/(1+|LL|):" _col(28) %9.2e `rg_c'
    display as txt "    V positive definite:" _col(28) `pd_c'
    display _n as txt "  Interpretation: under opposite-signs DGP the common"
    display as txt "  factor cannot fit the structure. Compare against the"
    display as txt "  separate-mode results above (which DID recover the truth)."
}
else {
    display as err "  factor(common) failed to estimate on this data (rc=`rc_common')."
    display as txt "  Failure-to-converge is itself an informative misspecification"
    display as txt "  signal."
}

* ============================================================================
* PART 4: Final verdict (recovery test only; misspec record is informational)
* ============================================================================

local n_pass = `pass_conv' + `pass_delta' + `pass_signs' + ///
    `pass_lT_mag' + `pass_lY_mag' + `pass_pi' + `pass_correction'

display _n as txt "{hline 70}"
if `n_pass' == 7 {
    display as res "OPPOSITE-SIGNS CERTIFICATION PASSED (7/7 criteria met)"
    display as txt "  hsmixture_joint recovers opposite-signs heterogeneity correctly."
    display as txt "  Identified quantities (delta, |lambda|, sign structure,"
    display as txt "  unordered mixture share) all match truth within tolerance."
    display as txt "  This validates recovery in the opposite-signs regime."
    display as txt "{hline 70}"
    capture log close _hsmix_opp
}
else {
    display as err "OPPOSITE-SIGNS CERTIFICATION FAILED (`n_pass'/7 criteria met)"
    display as err "  Opposite-signs recovery is not validated until this is resolved."
    display as txt "{hline 70}"
    capture log close _hsmix_opp
    error 9
}
