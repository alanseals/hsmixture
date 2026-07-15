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
*! Pass criteria (for the joint K=2 model). The individual signs of the
*! loadings and the type labels are NOT identified -- the likelihood is
*! invariant under (pi_1 <-> pi_2, lambda_T -> -lambda_T, lambda_Y ->
*! -lambda_Y) with the baseline constants absorbing the shift, so BFGS can
*! terminate at either labeling and both are the MLE. The criteria score
*! only identified quantities (as hsmixture_certification_opposite_signs.do
*! already does for its DGP):
*!   - e(converged) == 1 (BFGS converged AND |grad|/(1+|LL|) < 1e-5 AND V is PD)
*!   - |delta_hat - 0.5| < 3 * SE(delta)   (95%-band recovery; delta IS identified)
*!   - ||lambda_T_hat| - 0.7| < 3 * SE(lambda_T)   (magnitude)
*!   - ||lambda_Y_hat| - 1.0| < 3 * SE(lambda_Y)   (magnitude)
*!   - sign(lambda_T) == sign(lambda_Y)  (same-sign structure; the RELATIVE
*!     sign survives the swap, which flips both loadings together)
*!   - 0.30 < min(pi_1_hat, pi_2_hat) < 0.60  (true unordered share 0.45)
*!
*! Data-contract criteria (v2.4.0; PART 3): the estimators must REJECT
*! malformed person-period panels --
*!   - C1: rows after the outcome event (joint)
*!   - C2: omitted riskset() when treated persons have post-event rows (joint)
*!   - C3: treated indicator reverting 1 -> 0 (joint)
*!   - C4: treated == 1 before the treatment event (joint)
*!   - C5: rows after the outcome event (single-equation hsmixture)
*!   - C6: rows after the outcome event (bivariate)
*!   - C7: omitted riskset() (bivariate)
*!   - C8: time() accepts a row-shuffled CONFORMING panel (no false rejection)
*!   - C9: time() still rejects a row-shuffled corrupted panel
*!
*! Authors: Jonghoon Park and R. Alan Seals (Auburn University)

clear all
set seed 20260501
set more off

capture log close _hsmixture_cert
display as txt "Working directory: `c(pwd)'"
log using "hsmixture_certification_log.txt", text replace name(_hsmixture_cert)

* Provenance record: what code this certification actually ran on. `which`
* prints the resolved .ado path and its version starbang. Copy this block's
* output (plus the git commit SHA of the package) into CERTIFICATION.md
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

* Tolerance: 3 SEs for point estimates. pi has a hard band [0.30, 0.60]
* because its SE under the softmax parameterization is awkward to compute
* on the same scale.
*
* Label-swap invariance. The individual SIGN of each loading and the type
* labels are not separately identified. The likelihood is invariant under
*     (pi_1 <-> pi_2,  lambda_T -> -lambda_T,  lambda_Y -> -lambda_Y,
*      cons_T -> cons_T + lambda_T,  cons_Y -> cons_Y + lambda_Y),
* because v_1=0 / v_2=1 fixes the mass-point SCALE but not which type
* carries the shift. Both labelings are the same mixture and the same MLE,
* and BFGS can terminate at either depending on the platform's
* floating-point path. The criteria below therefore test what IS
* identified: delta, the loading MAGNITUDES, the RELATIVE sign of the two
* loadings (invariant because the swap flips both), and the unordered
* mixture share. This mirrors hsmixture_certification_opposite_signs.do,
* which tests the same invariants for the opposite-signs DGP.
*
* The relative-sign test is not a weakening: this DGP is same-sign
* (sigma_T=0.7, sigma_Y=1.0 both positive), that structure survives the
* swap, and nothing previously tested it.
local pass_conv  = (`conv' == 1)
local pass_delta = abs(`delta_hat' - `delta_true') < 3 * `se_delta'
local pass_lT    = abs(abs(`lT_hat') - abs(`sigma_T_true')) < 3 * `se_lT'
local pass_lY    = abs(abs(`lY_hat') - abs(`sigma_Y_true')) < 3 * `se_lY'
local pass_signs = (sign(`lT_hat') == sign(`lY_hat')) & (sign(`lT_hat') != 0)
local pi_1_hat   = 1 - `pi_2_hat'
local pi_min_hat = min(`pi_1_hat', `pi_2_hat')
local pass_pi    = (`pi_min_hat' > 0.30) & (`pi_min_hat' < 0.60)

display _n as txt "{hline 70}"
display as txt "CERTIFICATION RESULTS (joint K=2)"
display as txt "{hline 70}"
display as txt ///
    "  Test                              Estimate    SE       Truth   Status"
display as txt "{hline 70}"

* Report the raw signed estimates first so the log records which labeling
* BFGS landed on. Either labeling is the MLE; the criteria below score the
* swap-invariant quantities.
display as txt "  Labeling reached: lambda_T =" %8.4f `lT_hat' ///
    ", lambda_Y =" %8.4f `lY_hat' ", pi_2 =" %6.4f `pi_2_hat'
display as txt "  (signs and type labels are not identified; magnitudes," ///
    " relative sign,"
display as txt "   delta and the unordered share are.)" _n

* delta
local status_delta = cond(`pass_delta', "{res}PASS", "{err}FAIL")
display as txt "  delta (treatment effect)" _col(35) %8.4f `delta_hat' ///
    "  " %6.4f `se_delta' "  " %6.2f `delta_true' "    `status_delta'"

* |lambda_T|
local status_lT = cond(`pass_lT', "{res}PASS", "{err}FAIL")
display as txt "  |lambda_T| (treat loading)" _col(35) %8.4f abs(`lT_hat') ///
    "  " %6.4f `se_lT' "  " %6.2f abs(`sigma_T_true') "    `status_lT'"

* |lambda_Y|
local status_lY = cond(`pass_lY', "{res}PASS", "{err}FAIL")
display as txt "  |lambda_Y| (outcome loading)" _col(35) %8.4f abs(`lY_hat') ///
    "  " %6.4f `se_lY' "  " %6.2f abs(`sigma_Y_true') "    `status_lY'"

* Relative sign: this DGP is same-sign, and that structure is identified.
local status_signs = cond(`pass_signs', "{res}PASS", "{err}FAIL")
display as txt "  lambda_T, lambda_Y same sign" _col(40) ///
    cond(`pass_signs', "yes", "no") "     (truth: yes)    `status_signs'"

* Unordered mixture share
local status_pi = cond(`pass_pi', "{res}PASS", "{err}FAIL")
display as txt "  min(pi_1, pi_2)" _col(35) %8.4f `pi_min_hat' ///
    "       --   " %6.2f min(`pi_1_true', 1 - `pi_1_true') "    `status_pi'"

* Strict convergence
local status_conv = cond(`pass_conv', "{res}PASS", "{err}FAIL")
display _n as txt "  Strict convergence (e(converged)==1):  `status_conv'"
display as txt "    Gradient norm =" %9.2e `gn'
display as txt "    |grad|/(1+|LL|) =" %9.2e e(rel_grad)
display as txt "    V positive definite = `pd'"
* Eigenvalue provenance for the v_pd verdict. e(v_pd) applies a
* scale-relative test (min eigenvalue > 0 and > 1e-12 * max eigenvalue;
* v2.4.0 -- the old absolute 1e-8 floor was scale-dependent and conflated
* an ill-conditioned V with the I*1e-20 placeholder). e(v_mineig)
* is the number tested and e(v_scaffold) says whether e(V) is the I*1e-20
* placeholder posted after Hessian-inversion failure (in which case the
* SEs are meaningless) rather than a real inverted Hessian. Printing both
* makes a v_pd=0 verdict diagnosable instead of opaque.
display as txt "    V min eigenvalue =" %9.2e e(v_mineig)
display as txt "    V is placeholder/scaffold = " e(v_scaffold)

* ============================================================================
* PART 3: Data-contract validation tests (v2.4.0 row contracts)
* ============================================================================
* Each test corrupts a copy of the certification panel in a specific way and
* asserts the estimator REFUSES it (exit code 198). These guard the row-level
* contracts the likelihood assumes: absorbing outcome (rows stop at the
* outcome event), one-time treatment (post-event rows excluded via riskset),
* and treated-indicator consistency (absorbing, never leads its own event).
* A contract test FAILS if the estimator silently accepts the bad panel.

display _n as txt "{hline 70}"
display as txt "Data-contract validation tests"
display as txt "{hline 70}"

* --- C1: rows after the outcome event must be rejected ---
* The duplicate row zeroes treat_event so a person whose treatment and
* outcome events coincide on the last row cannot trip the (pre-existing,
* earlier-running) multiple-treatment-event check; only the new
* absorbing-outcome row check can fire.
preserve
quietly {
    bysort id (period): gen byte _lastrow = (_n == _N)
    expand 2 if _lastrow & outcome_event == 1, gen(_dup)
    replace outcome_event = 0 if _dup == 1
    replace treat_event = 0 if _dup == 1
    replace period = period + 1 if _dup == 1
    sort id period
    count if _dup == 1
    local n_dup = r(N)
}
display as txt "    (post-outcome rows added: `n_dup')"
capture hsmixture_joint (treat_event = x1) ///
    (outcome_event = x1, treat(treated)) ///
    , id(id) k(2) riskset(treat_at_risk) nolog iterate(5)
local pass_c1 = (_rc == 198) & (`n_dup' > 0)
restore
local status_c1 = cond(`pass_c1', "{res}PASS", "{err}FAIL")
display as txt "  C1: post-outcome rows rejected" _col(55) "`status_c1'"

* --- C2: omitting riskset() on data with post-treatment-event rows must be
*     rejected (treated persons here are observed after their event, so an
*     all-ones treatment risk set would be mis-specified) ---
capture hsmixture_joint (treat_event = x1) ///
    (outcome_event = x1, treat(treated)) ///
    , id(id) k(2) nolog iterate(5)
local pass_c2 = (_rc == 198)
local status_c2 = cond(`pass_c2', "{res}PASS", "{err}FAIL")
display as txt "  C2: omitted riskset() rejected" _col(55) "`status_c2'"

* --- C3: treated reverting 1 -> 0 must be rejected ---
preserve
quietly {
    bysort id (period): gen byte _rev = (_n == _N) & (_N >= 2) & ///
        treated == 1 & treated[_N-1] == 1
    count if _rev
    local n_rev = r(N)
    replace treated = 0 if _rev
    * The reverted final rows leave the treatment risk set untouched, so
    * only the reversion check can fire.
}
display as txt "    (reverted final rows: `n_rev')"
capture hsmixture_joint (treat_event = x1) ///
    (outcome_event = x1, treat(treated)) ///
    , id(id) k(2) riskset(treat_at_risk) nolog iterate(5)
local pass_c3 = (_rc == 198) & (`n_rev' > 0)
restore
local status_c3 = cond(`pass_c3', "{res}PASS", "{err}FAIL")
display as txt "  C3: reverting treated rejected" _col(55) "`status_c3'"

* --- C4: treated == 1 before the treatment event must be rejected ---
* The corruption flips treated to 1 on BOTH the pre-event row and the event
* row, so the corrupted sequence (0,...,0,1,1,1,...) stays monotone and
* cannot trip the reversion check; only the before-event check can fire.
* (Flipping the pre-event row alone would create 0,...,1,0,1,... under the
* DGP's lagged convention and the reversion check would fire first,
* leaving the before-event check unexercised.)
preserve
quietly {
    bysort id (period): gen byte _preflip = (treat_event[_n+1] == 1)
    count if _preflip
    local n_pre = r(N)
    replace treated = 1 if _preflip | treat_event == 1
}
display as txt "    (pre-event rows flipped: `n_pre')"
capture hsmixture_joint (treat_event = x1) ///
    (outcome_event = x1, treat(treated)) ///
    , id(id) k(2) riskset(treat_at_risk) nolog iterate(5)
local pass_c4 = (_rc == 198) & (`n_pre' > 0)
restore
local status_c4 = cond(`pass_c4', "{res}PASS", "{err}FAIL")
display as txt "  C4: treated before its event rejected" _col(55) "`status_c4'"

* --- C5: single-equation hsmixture rejects post-outcome rows too ---
preserve
quietly {
    bysort id (period): gen byte _lastrow = (_n == _N)
    expand 2 if _lastrow & outcome_event == 1, gen(_dup)
    replace outcome_event = 0 if _dup == 1
    replace period = period + 1 if _dup == 1
    sort id period
    count if _dup == 1
    local n_dup5 = r(N)
}
capture hsmixture outcome_event x1, id(id) k(2) nolog iterate(5)
local pass_c5 = (_rc == 198) & (`n_dup5' > 0)
restore
local status_c5 = cond(`pass_c5', "{res}PASS", "{err}FAIL")
display as txt "  C5: post-outcome rows rejected (hsmixture)" _col(55) "`status_c5'"

* --- C6: hsmixture_bivariate rejects post-outcome rows too ---
* Same construction as C1 (treat_event zeroed on the duplicate row so only
* the absorbing-outcome row check can fire).
preserve
quietly {
    bysort id (period): gen byte _lastrow = (_n == _N)
    expand 2 if _lastrow & outcome_event == 1, gen(_dup)
    replace outcome_event = 0 if _dup == 1
    replace treat_event = 0 if _dup == 1
    replace period = period + 1 if _dup == 1
    sort id period
    count if _dup == 1
    local n_dup6 = r(N)
}
capture hsmixture_bivariate (treat_event = x1) ///
    (outcome_event = x1, treat(treated)) ///
    , id(id) riskset(treat_at_risk) nolog iterate(5) nstarts(1)
local pass_c6 = (_rc == 198) & (`n_dup6' > 0)
restore
local status_c6 = cond(`pass_c6', "{res}PASS", "{err}FAIL")
display as txt "  C6: post-outcome rows rejected (bivariate)" _col(55) "`status_c6'"

* --- C7: hsmixture_bivariate rejects omitted riskset() too ---
capture hsmixture_bivariate (treat_event = x1) ///
    (outcome_event = x1, treat(treated)) ///
    , id(id) nolog iterate(5) nstarts(1)
local pass_c7 = (_rc == 198)
local status_c7 = cond(`pass_c7', "{res}PASS", "{err}FAIL")
display as txt "  C7: omitted riskset() rejected (bivariate)" _col(55) "`status_c7'"

* --- C8/C9: time() supplies within-person order independent of row order ---
* The checks otherwise read order from physical row position. C8 shuffles a
* CONFORMING panel and asserts time(period) prevents a false rejection; C9
* shuffles a CORRUPTED panel and asserts time(period) still catches it.
* Together these prove the contract checks follow time(), not row order.
preserve
quietly {
    set seed 606
    gen double _shuf = runiform()
    sort _shuf
}
capture hsmixture_joint (treat_event = x1) ///
    (outcome_event = x1, treat(treated)) ///
    , id(id) k(2) riskset(treat_at_risk) time(period) nolog iterate(5) nstarts(1)
local rc_c8 = _rc
local pass_c8 = (`rc_c8' == 0)
restore
local status_c8 = cond(`pass_c8', "{res}PASS", "{err}FAIL")
display as txt "  C8: time() accepts shuffled conforming panel (rc=`rc_c8')" ///
    _col(55) "`status_c8'"

preserve
quietly {
    bysort id (period): gen byte _lastrow = (_n == _N)
    expand 2 if _lastrow & outcome_event == 1, gen(_dup)
    replace outcome_event = 0 if _dup == 1
    replace treat_event = 0 if _dup == 1
    replace period = period + 1 if _dup == 1
    set seed 707
    gen double _shuf = runiform()
    sort _shuf
}
capture hsmixture_joint (treat_event = x1) ///
    (outcome_event = x1, treat(treated)) ///
    , id(id) k(2) riskset(treat_at_risk) time(period) nolog iterate(5) nstarts(1)
local pass_c9 = (_rc == 198)
restore
local status_c9 = cond(`pass_c9', "{res}PASS", "{err}FAIL")
display as txt "  C9: time() still rejects shuffled bad panel" _col(55) "`status_c9'"

* ============================================================================
* PART 4: Final verdict
* ============================================================================

local n_pass = `pass_conv' + `pass_delta' + `pass_lT' + `pass_lY' + ///
    `pass_signs' + `pass_pi' ///
    + `pass_c1' + `pass_c2' + `pass_c3' + `pass_c4' + `pass_c5' ///
    + `pass_c6' + `pass_c7' + `pass_c8' + `pass_c9'

display _n as txt "{hline 70}"
if `n_pass' == 15 {
    display as res "CERTIFICATION PASSED (15/15 criteria met)"
    display as txt "  hsmixture_joint recovers the identified parameters of the true DGP"
    display as txt "  within tolerance (6/6) and the data-contract validations reject"
    display as txt "  malformed panels across all three estimators, without false"
    display as txt "  rejections (9/9)."
    display as txt "{hline 70}"
    capture log close _hsmixture_cert
}
else {
    display as err "CERTIFICATION FAILED (`n_pass'/15 criteria met)"
    display as err "  Recovery: `pass_conv'+`pass_delta'+`pass_lT'+`pass_lY'+`pass_signs'+`pass_pi' of 6;"
    display as err "  contracts: `pass_c1'+`pass_c2'+`pass_c3'+`pass_c4'+`pass_c5'+`pass_c6'+`pass_c7'+`pass_c8'+`pass_c9' of 9."
    display as err "  Inspect the model output above for details."
    display as err "  Possible causes: BFGS hit a local mode (re-run with from()),"
    display as err "  rank-deficient Hessian (try smaller K), a genuine DGP issue, or a"
    display as err "  validation regression (a contract test accepting a malformed panel)."
    display as txt "{hline 70}"
    capture log close _hsmixture_cert
    error 9
}
