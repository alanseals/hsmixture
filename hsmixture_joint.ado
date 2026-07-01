*! version 2.3.2  30jun2026
*! Heckman-Singer Joint Timing-of-Events Model
*! Two-equation discrete-time hazard with correlated unobserved heterogeneity.
*! Supports two factor structures via factor() option:
*!   common   - one shared loading lambda (manuscript-style one-factor MPH)
*!   separate - free loadings lambda_T, lambda_Y (default; admits opposite
*!              signs)
*!
*! Authors: Jonghoon Park and R. Alan Seals (Auburn University)
*!
*! References:
*!   Abbring, J.H. and G.J. van den Berg. 2003. The nonparametric
*!     identification of treatment effects in duration models.
*!     Econometrica 71(5): 1491-1517.
*!   Heckman, J.J. and B. Singer. 1984. Econometrica 52(2): 271-320.

program hsmixture_joint, eclass sortpreserve
    version 14

    if replay() {
        if "`e(cmd)'" != "hsmixture_joint" {
            error 301
        }
        Display `0'
        exit
    }

    * Parse the compound syntax
    * Syntax: (treat_dep = treat_indep) (outcome_dep = outcome_indep, treat(varname))
    *         , id() k() [options]

    gettoken eq1 0 : 0, parse("(") match(paren)
    gettoken eq2 0 : 0, parse("(") match(paren)

    * Parse equation 1 (treatment process)
    gettoken treat_dep eq1_rest : eq1, parse("=")
    gettoken eq_sign treat_indep : eq1_rest, parse("=")
    local treat_dep = strtrim("`treat_dep'")
    local treat_indep = strtrim("`treat_indep'")

    * Parse equation 2 (outcome process with treatment indicator)
    local treat_var ""

    * Check if treat() is specified
    if regexm("`eq2'", "treat\(([^\)]+)\)") {
        local treat_var = regexs(1)
        * Remove treat() from eq2
        local eq2 = regexr("`eq2'", ",[ ]*treat\([^\)]+\)", "")
    }

    gettoken outcome_dep eq2_rest : eq2, parse("=")
    gettoken eq_sign outcome_indep : eq2_rest, parse("=")
    local outcome_dep = strtrim("`outcome_dep'")
    local outcome_indep = strtrim("`outcome_indep'")

    * Parse remaining options. PRisk() is a deprecated alias for RISKset()
    * retained for back-compat with v2.0.0 production scripts.
    *
    * FACTor() selects the heterogeneity structure:
    *   separate (default): two free signed loadings lambda_T, lambda_Y
    *     for the treatment and outcome equations. Per-type shifts
    *     (lambda_T*v_k, lambda_Y*v_k) lie on a 1-D locus through the
    *     origin but the locus direction is free. This admits opposite
    *     signs (negative correlation between treatment-prone and
    *     outcome-prone latent types).
    *   common: one free signed loading lambda shared by both equations,
    *     so the type-k shift is (lambda*v_k, lambda*v_k). This is the
    *     classical Heckman-Singer one-factor mixed-proportional-hazards
    *     model: the same exp(eta_i) raises both hazards. Correlation
    *     between latent treatment and outcome propensity is mechanically
    *     positive.
    *
    * NSTarts(integer 7) selects the number of starting-value configurations
    * for multistart (best log-likelihood retained). Mixture surfaces are
    * multimodal, so single-start is unsafe at any K. Under factor(common)
    * the configuration grid is smaller (one fewer loading parameter);
    * the default 7 is capped to 6 in that mode.
    *
    * DIFFicult, TRace, GRADient, HESSian, TECHnique(), TOLerance(),
    * LTOLerance(), NRTOLerance() are legacy `ml model d0` options from
    * v2.0.0. They are silently accepted but ignored by the Mata optimize()
    * implementation, which uses BFGS with vtol=1e-10 + nrtol=1e-5.
    syntax [if] [in], ///
        ID(varname) ///
        [K(integer 2) ///
         FACTor(string) ///
         RISKset(varname) ///
         PRisk(varname) ///
         FROM(name) ///
         ITERate(integer 100) ///
         NSTarts(integer 7) ///
         noLOG ///
         Level(cilevel) ///
         DIFFicult ///
         TRace ///
         GRADient ///
         HESSian ///
         TECHnique(string) ///
         TOLerance(real 1e-6) ///
         LTOLerance(real 1e-7) ///
         NRTOLerance(real 1e-5)]

    * Resolve factor() mode. Default separate for back-compat with v2.x callers.
    if "`factor'" == "" {
        local factor "separate"
    }
    if !inlist("`factor'", "common", "separate") {
        display as error ///
            "factor() must be 'common' or 'separate' (got '`factor'')"
        exit 198
    }

    * Resolve prisk/riskset alias
    if "`prisk'" != "" & "`riskset'" != "" {
        display as error "cannot specify both prisk() and riskset(); use riskset()"
        exit 198
    }
    if "`prisk'" != "" & "`riskset'" == "" {
        local riskset "`prisk'"
        display as txt "(prisk() is a deprecated alias for riskset(); please update callers)"
    }

    * Mark sample. Include riskset so that observations with missing
    * riskset values are excluded from the estimation sample rather than
    * silently entering Mata as missings (which would multiply the
    * pregnancy log-likelihood contribution by missing).
    marksample touse
    markout `touse' `id' `treat_dep' `outcome_dep' `treat_var' `riskset'

    * Validate inputs
    if "`treat_var'" == "" {
        display as error "treat() required in second equation"
        exit 198
    }

    if `k' < 2 {
        display as error "k() must be at least 2"
        exit 198
    }

    * Expand wildcards and exclude obs with missing covariates from touse
    * BEFORE the data-contract assertions. If a row has missing covariates it
    * will not enter the estimation sample, so it must not be checked by the
    * 0/1 / single-event / riskset assertions either.
    unab treat_indep_exp : `treat_indep'
    unab outcome_indep_exp : `outcome_indep'
    markout `touse' `treat_indep_exp' `outcome_indep_exp'

    * Stage 3: data-contract validation. Errors (not warnings) on violations
    * because these break the likelihood. Risk-set rows are dropped via
    * markout above; here we check the values that survive.
    capture assert inlist(`treat_dep', 0, 1) if `touse'
    if _rc {
        display as error ///
            "treatment event variable `treat_dep' must be 0/1 in the estimation sample"
        exit 198
    }
    capture assert inlist(`outcome_dep', 0, 1) if `touse'
    if _rc {
        display as error ///
            "outcome event variable `outcome_dep' must be 0/1 in the estimation sample"
        exit 198
    }
    capture assert inlist(`treat_var', 0, 1) if `touse'
    if _rc {
        display as error ///
            "treatment indicator `treat_var' must be 0/1 in the estimation sample"
        exit 198
    }
    if "`riskset'" != "" {
        capture assert inlist(`riskset', 0, 1) if `touse'
        if _rc {
            display as error "riskset() variable `riskset' must be 0/1 in the estimation sample"
            exit 198
        }
    }

    * Stage 3 (extended): duration-data structural validation. The likelihood
    * assumes one-time treatment events (each person becomes treated at most
    * once) and an absorbing outcome (each person experiences the outcome at
    * most once). Violations of these assumptions will distort delta and the
    * mixture parameters in subtle ways. Errors here protect users from
    * silently feeding bad person-period panels into the estimator.
    tempvar __hsj_evcount
    quietly bysort `id': egen `__hsj_evcount' = total(`treat_dep') if `touse'
    capture assert `__hsj_evcount' <= 1 if `touse'
    if _rc {
        display as error ///
            "treatment event variable `treat_dep' has more than one event for some id"
        display as error ///
            "  hsmixture_joint requires one-time treatment timing (each person treated at most once)"
        exit 198
    }
    drop `__hsj_evcount'

    tempvar __hsj_evcount
    quietly bysort `id': egen `__hsj_evcount' = total(`outcome_dep') if `touse'
    capture assert `__hsj_evcount' <= 1 if `touse'
    if _rc {
        display as error ///
            "outcome event variable `outcome_dep' has more than one event for some id"
        display as error ///
            "  hsmixture_joint requires an absorbing outcome (each person experiences the outcome at most once)"
        exit 198
    }
    drop `__hsj_evcount'

    * If riskset is specified, require treat_event=1 only when riskset=1.
    if "`riskset'" != "" {
        capture assert `riskset' == 1 if `treat_dep' == 1 & `touse'
        if _rc {
            display as error ///
                "treatment events occur outside the treatment risk set (`treat_dep'=1 with `riskset'=0)"
            display as error ///
                "  treatment can only happen when the person is at risk (`riskset'=1)"
            exit 198
        }
    }

    * Compile Mata functions (silently skips if already loaded)
    capture findfile hsmixture_joint_mata.do
    if _rc {
        display as error "hsmixture_joint_mata.do not found on adopath"
        exit 601
    }
    capture quietly do "`r(fn)'"

    * ----------------------------------------------------------------
    * Set globals for likelihood evaluator (cleaned up on all exit paths)
    * ----------------------------------------------------------------
    global HSJ_K = `k'
    global HSJ_id "`id'"
    global HSJ_treat_event "`treat_dep'"
    global HSJ_outcome_event "`outcome_dep'"
    global HSJ_treat "`treat_var'"
    global HSJ_riskset "`riskset'"
    global HSJ_factor "`factor'"
    if "`log'" == "nolog" {
        global HSJ_nolog "1"
    }
    else {
        global HSJ_nolog ""
    }

    * Wrap estimation in capture for guaranteed cleanup on break/error
    capture noisily {

    * ----------------------------------------------------------------
    * Fit initial GLMs to get starting values and expand factor vars
    * ----------------------------------------------------------------
    if "`log'" != "nolog" {
        display _n as txt "Obtaining starting values from separate GLM models..."
    }

    * Treatment equation
    if "`riskset'" != "" {
        quietly glm `treat_dep' `treat_indep' if `touse' & `riskset' == 1, ///
            family(binomial) link(cloglog) nolog
    }
    else {
        quietly glm `treat_dep' `treat_indep' if `touse', ///
            family(binomial) link(cloglog) nolog
    }
    tempname b_treat
    matrix `b_treat' = e(b)

    * `treat_indep_exp' and `outcome_indep_exp' were expanded above (before
    * the data-contract asserts) so touse already excludes obs with missing
    * covariates.
    local treat_vars_mata "`treat_indep_exp'"

    * Outcome equation
    quietly glm `outcome_dep' `outcome_indep' if `touse', ///
        family(binomial) link(cloglog) nolog
    tempname b_outcome
    matrix `b_outcome' = e(b)

    local outcome_vars_mata "`outcome_indep_exp'"

    * Store variable names in globals for Mata
    global HSJ_treat_vars_mata "`treat_vars_mata'"
    global HSJ_outcome_vars_mata "`outcome_vars_mata'"

    * Get starting value for delta from simple cloglog with treatment
    quietly glm `outcome_dep' `treat_var' `outcome_indep' if `touse', ///
        family(binomial) link(cloglog) nolog
    local delta_init = _b[`treat_var']

    * ----------------------------------------------------------------
    * Build multistart starting-value configurations
    * ----------------------------------------------------------------
    * Layout under v_2=1 normalization:
    *   factor(separate): [beta_T | beta_Y | delta | lambda_T | lambda_Y |
    *                      v_3..v_K | eta_2..eta_K]
    *   factor(common):   [beta_T | beta_Y | delta | lambda |
    *                      v_3..v_K | eta_2..eta_K]
    * Loadings are signed real-line parameters. v_3..v_K (when K >= 3) are
    * free mass points. eta_k feeds the softmax mixture.
    *
    * Multistart probes regions of the loading space:
    *   separate: 7 configs. Configs 1-6 vary (lambda_T, lambda_Y) magnitudes
    *     and (for K>=3) v_3, with balanced initial mixture. Config 7 probes
    *     the small-loading + asymmetric-mixture region.
    *   common: 6 configs. The loading is one-dimensional, so we don't need
    *     the asymmetric (lambda_T != lambda_Y) probes from configs 3-4.
    *
    * Best log-likelihood across all starts wins.

    if "`factor'" == "common" {
        * Common-loading configs. lambda is the single shared loading.
        local sv_l_1   = 0.3
        local sv_v3_1  = 2.0
        local sv_eta2_1 = 0

        local sv_l_2   = 0.5
        local sv_v3_2  = 2.0
        local sv_eta2_2 = 0

        local sv_l_3   = 1.0
        local sv_v3_3  = 2.5
        local sv_eta2_3 = 0

        local sv_l_4   = 1.0
        local sv_v3_4  = 3.0
        local sv_eta2_4 = 0

        local sv_l_5   = 1.5
        local sv_v3_5  = 4.0
        local sv_eta2_5 = 0

        * Config 6: small-loading, asymmetric-mixture probe
        local sv_l_6   = 0.05
        local sv_v3_6  = 2.0
        local sv_eta2_6 = 1.9

        local n_configs = min(`nstarts', 6)
        if `n_configs' < 1 local n_configs = 1

        forvalues s = 1/`n_configs' {
            tempname b0_`s'
            matrix `b0_`s'' = `b_treat', `b_outcome', `delta_init'

            * Single shared loading
            matrix `b0_`s'' = `b0_`s'', `sv_l_`s''

            * v_3..v_K (empty for K=2)
            forvalues j = 3/`k' {
                local v_init = `sv_v3_`s'' + (`j' - 3)
                matrix `b0_`s'' = `b0_`s'', `v_init'
            }

            * eta_2..eta_K
            matrix `b0_`s'' = `b0_`s'', `sv_eta2_`s''
            forvalues j = 3/`k' {
                matrix `b0_`s'' = `b0_`s'', 0
            }
        }
    }
    else {
        * Separate-loading configs (current default). Per-start
        * (lambda_T, lambda_Y, v_3, eta_2) values; v_4..v_K built off v_3.
        local sv_lT_1   = 0.3
        local sv_lY_1   = 0.3
        local sv_v3_1   = 2.0
        local sv_eta2_1 = 0

        local sv_lT_2   = 0.5
        local sv_lY_2   = 0.5
        local sv_v3_2   = 2.0
        local sv_eta2_2 = 0

        local sv_lT_3   = 1.0
        local sv_lY_3   = 0.5
        local sv_v3_3   = 2.0
        local sv_eta2_3 = 0

        local sv_lT_4   = 0.5
        local sv_lY_4   = 1.0
        local sv_v3_4   = 2.5
        local sv_eta2_4 = 0

        local sv_lT_5   = 1.0
        local sv_lY_5   = 1.0
        local sv_v3_5   = 3.0
        local sv_eta2_5 = 0

        local sv_lT_6   = 1.5
        local sv_lY_6   = 1.5
        local sv_v3_6   = 4.0
        local sv_eta2_6 = 0

        * Config 7: small-loading, asymmetric-mixture probe
        local sv_lT_7   = 0.04
        local sv_lY_7   = 0.06
        local sv_v3_7   = 2.0
        local sv_eta2_7 = 1.9

        local n_configs = min(`nstarts', 7)
        if `n_configs' < 1 local n_configs = 1

        forvalues s = 1/`n_configs' {
            tempname b0_`s'
            matrix `b0_`s'' = `b_treat', `b_outcome', `delta_init'

            * lambda_T, lambda_Y (signed real, no exp transform)
            matrix `b0_`s'' = `b0_`s'', `sv_lT_`s'', `sv_lY_`s''

            * v_3..v_K (empty for K=2)
            forvalues j = 3/`k' {
                local v_init = `sv_v3_`s'' + (`j' - 3)
                matrix `b0_`s'' = `b0_`s'', `v_init'
            }

            * eta_2..eta_K
            matrix `b0_`s'' = `b0_`s'', `sv_eta2_`s''
            forvalues j = 3/`k' {
                matrix `b0_`s'' = `b0_`s'', 0
            }
        }
    }

    * Use user-supplied starting values if provided (replaces config 1 only)
    if "`from'" != "" {
        matrix `b0_1' = `from'
    }

    * Display header
    if "`log'" != "nolog" {
        display _n as txt "Joint Timing-of-Events Model with Heckman-Singer Heterogeneity"
        display as txt "K = `k' mass points (v_1 = 0, v_2 = 1 normalized)"
        if "`factor'" == "common" {
            display as txt "Factor structure: common (one shared loading lambda)"
        }
        else {
            display as txt "Factor structure: separate (lambda_T, lambda_Y free)"
        }
        display as txt "Multistart: `n_configs' starting configurations"
        if "`riskset'" != "" {
            display as txt "Treatment risk set: `riskset' == 1"
        }
        display as txt "{hline 70}"
        display as txt "Treatment equation: `treat_dep' = `treat_indep'"
        display as txt "Outcome equation:   `outcome_dep' = `outcome_indep' + delta*`treat_var'"
        display as txt "{hline 70}"
    }

    * Ensure data is sorted by id for Mata panelsetup
    sort `id'

    * Initialize Mata cache
    capture mata: _hsj_cleanup()

    * ----------------------------------------------------------------
    * Run optimization with Mata optimize()
    * ----------------------------------------------------------------

    local best_ll = .
    local best_converged = 0

    forvalues s = 1/`n_configs' {
        if "`log'" != "nolog" {
            if "`factor'" == "common" {
                if `k' >= 3 {
                    display _n as txt "--- Starting configuration `s'/`n_configs': " ///
                        "lambda=`sv_l_`s'', " ///
                        "v_3=`sv_v3_`s'', eta_2=`sv_eta2_`s'' ---"
                }
                else {
                    display _n as txt "--- Starting configuration `s'/`n_configs': " ///
                        "lambda=`sv_l_`s'', " ///
                        "eta_2=`sv_eta2_`s'' ---"
                }
            }
            else {
                if `k' >= 3 {
                    display _n as txt "--- Starting configuration `s'/`n_configs': " ///
                        "lambda_T=`sv_lT_`s'', lambda_Y=`sv_lY_`s'', " ///
                        "v_3=`sv_v3_`s'', eta_2=`sv_eta2_`s'' ---"
                }
                else {
                    display _n as txt "--- Starting configuration `s'/`n_configs': " ///
                        "lambda_T=`sv_lT_`s'', lambda_Y=`sv_lY_`s'', " ///
                        "eta_2=`sv_eta2_`s'' ---"
                }
            }
        }

        matrix __hsj_start = `b0_`s''
        scalar __hsj_has_result = 0
        scalar __hsj_converged = 0

        capture noisily mata: _hsj_run_optimize("__hsj_start", `iterate')

        * Check if Mata posted a result (finite LL)
        capture confirm scalar __hsj_has_result
        if _rc == 0 & scalar(__hsj_has_result) == 1 {
            local ll_s = scalar(__hsj_ll)
            if "`log'" != "nolog" {
                display as txt "  Log-likelihood: " %12.4f `ll_s'
            }

            if `ll_s' > `best_ll' | `best_ll' == . {
                local best_ll = `ll_s'
                local best_start = `s'
                local best_converged = scalar(__hsj_converged)
                local best_ic = scalar(__hsj_ic)
                matrix __hsj_best_b = __hsj_b
                capture confirm matrix __hsj_V
                if _rc == 0 {
                    matrix __hsj_best_V = __hsj_V
                }
                else {
                    * Hessian inversion failed. Post a finite scaffold V so
                    * `ereturn post` succeeds and callers can still read
                    * point estimates from e(b). A missing-value V (J(np,np,.))
                    * is rejected by Stata's `ereturn post`, which would
                    * leave the entire fit unreadable. The eigenvalue check
                    * below sets e(v_pd)=0 and the Display block warns
                    * prominently; callers should gate on e(converged)==1
                    * (which requires V positive definite) before trusting SEs.
                    local _np = colsof(__hsj_best_b)
                    matrix __hsj_best_V = I(`_np') * 1e-20
                }
                capture confirm matrix __hsj_g
                if _rc == 0 {
                    matrix __hsj_best_g = __hsj_g
                }
                else {
                    capture matrix drop __hsj_best_g
                }
            }
        }
        else {
            if "`log'" != "nolog" {
                display as txt "  Configuration `s' failed to converge"
            }
        }

        capture matrix drop __hsj_start __hsj_b __hsj_V __hsj_g
        capture scalar drop __hsj_ll __hsj_has_result __hsj_converged __hsj_ic
    }

    * ----------------------------------------------------------------
    * Check that at least one configuration produced results
    * ----------------------------------------------------------------

    if `best_ll' == . {
        display as error "All starting configurations failed to converge"
        exit 430
    }

    if "`log'" != "nolog" {
        display _n as txt "Best result from configuration `best_start'" ///
            " (log-lik: " %12.4f `best_ll' ")"
    }

    * ----------------------------------------------------------------
    * Post results to ereturn
    * ----------------------------------------------------------------

    * Equation labels use the actual dependent-variable names so that
    * post-estimation `_b[depvar:varname]` calls match what the user wrote.
    * v_2 is fixed at 1 (scale normalization) and is not a free parameter,
    * so v starts at v_3 here.
    local colnames ""
    foreach v of local treat_indep_exp {
        local colnames "`colnames' `treat_dep':`v'"
    }
    local colnames "`colnames' `treat_dep':_cons"
    foreach v of local outcome_indep_exp {
        local colnames "`colnames' `outcome_dep':`v'"
    }
    local colnames "`colnames' `outcome_dep':_cons"
    if "`factor'" == "common" {
        local colnames "`colnames' delta:_cons lambda:_cons"
    }
    else {
        local colnames "`colnames' delta:_cons lambda_T:_cons lambda_Y:_cons"
    }
    forvalues j = 3/`k' {
        local colnames "`colnames' v_`j':_cons"
    }
    forvalues j = 2/`k' {
        local colnames "`colnames' eta_`j':_cons"
    }

    matrix colnames __hsj_best_b = `colnames'
    matrix colnames __hsj_best_V = `colnames'
    matrix rownames __hsj_best_V = `colnames'

    quietly count if `touse'
    local N_obs = r(N)

    * Person-level count for AIC/BIC. The IID unit in this mixture model is
    * the person; the latent type is integrated out at the person level (see
    * Mata _hsj_compute_ll). Using person-period N would inflate the BIC
    * penalty. e(N) follows Stata's row-count convention; e(N_persons) is
    * the statistically appropriate denominator for cross-K comparisons.
    tempvar __hsj_pid_tag
    quietly egen byte `__hsj_pid_tag' = tag(`id') if `touse'
    quietly count if `__hsj_pid_tag' == 1
    local N_persons = r(N)
    drop `__hsj_pid_tag'

    ereturn post __hsj_best_b __hsj_best_V, obs(`N_obs')
    ereturn scalar ll = `best_ll'
    ereturn scalar N_persons = `N_persons'

    * Clean up Mata cache
    capture mata: _hsj_cleanup()

    * Strict-convergence diagnostics. e(converged)=1 only if (a) BFGS reported
    * convergence, (b) the relative gradient |grad|/(1+|LL|) < 1e-5, and
    * (c) the variance matrix is positive definite. The BFGS flag alone is
    * preserved as e(converged_bfgs).
    *
    * The relative metric is used because absolute |grad| does not scale with
    * sample size. For LL on the order of 10^4 the locally flat optimum has
    * |grad| ~ 5e-3 by Mata's BFGS at vtol=1e-10/nrtol=1e-5; tightening
    * further sends BFGS into a contraction-reset loop. The relative metric
    * matches the form of Mata's own nrtol stopping criterion.
    *
    * Production callers gating on e(converged) will therefore distinguish
    * a clean interior optimum from a boundary "convergence" that satisfies
    * BFGS's relative-LL stopping rule but is not a stationary point.
    local strict_converged = `best_converged'
    local grad_norm = .
    local rel_grad = .
    local v_pd = 0

    capture confirm matrix __hsj_best_g
    if _rc == 0 {
        capture mata: st_numscalar("__hsj_gn", norm(st_matrix("__hsj_best_g")))
        if _rc == 0 {
            local grad_norm = scalar(__hsj_gn)
            local rel_grad = `grad_norm' / (1 + abs(`best_ll'))
            if `rel_grad' > 1e-5 local strict_converged = 0
        }
        else {
            local strict_converged = 0
        }
        capture scalar drop __hsj_gn
    }
    else {
        local strict_converged = 0
    }

    capture mata: st_numscalar("__hsj_meig", min(Re(eigenvalues(st_matrix("e(V)")))))
    if _rc == 0 {
        if scalar(__hsj_meig) > 1e-8 local v_pd = 1
    }
    if !`v_pd' local strict_converged = 0
    capture scalar drop __hsj_meig

    * Store results
    ereturn local cmd "hsmixture_joint"
    ereturn local cmdline "hsmixture_joint `0'"
    ereturn local treat_depvar "`treat_dep'"
    ereturn local outcome_depvar "`outcome_dep'"
    ereturn local treat_var "`treat_var'"
    ereturn local idvar "`id'"
    ereturn local riskset_var "`riskset'"
    ereturn local factor "`factor'"
    ereturn scalar K = `k'
    ereturn scalar k = colsof(e(b))
    * Force estimates stats / IC machinery to use the *design* parameter
    * count rather than rank(V). When the Hessian is rank-deficient
    * (e.g., a mass-point parameter pinned to a boundary, missing SE),
    * Stata's default would compute df = rank(V) < colsof(e(b)) and
    * report an AIC/BIC that disagrees with the package's own display.
    ereturn scalar rank = colsof(e(b))
    ereturn scalar df_m = colsof(e(b))
    ereturn scalar converged = `strict_converged'
    ereturn scalar converged_bfgs = `best_converged'
    ereturn scalar grad_norm = `grad_norm'
    ereturn scalar rel_grad = `rel_grad'
    ereturn scalar v_pd = `v_pd'
    ereturn scalar ic = `best_ic'
    ereturn scalar n_starts = `n_configs'
    ereturn scalar best_start = `best_start'

    * Post gradient at best-start optimum (consumed by
    * hsmixture_joint_postestimation, convergence)
    capture confirm matrix __hsj_best_g
    if _rc == 0 {
        matrix colnames __hsj_best_g = `colnames'
        ereturn matrix gradient = __hsj_best_g
    }

    * Compute and store transformed parameters. Loadings are signed real-line
    * (Heckman-Singer factor loadings under v_2=1 scale normalization). The
    * v2.0.0 sigma_T/sigma_Y/sigma_P/sigma_D scalars are preserved as aliases
    * pointing to the same value, since under v_2=1 the "factor loading" and
    * the "type-2 risk shift" coincide.
    *
    * Under factor(common), the single shared loading lambda is exposed
    * directly and also as both _T and _Y aliases (and corresponding sigma
    * aliases). Existing v2.x callers reading e(lambda_T) and e(lambda_Y)
    * therefore work without modification — they just see equal values.
    tempname lambda_T lambda_Y delta_est hr
    scalar `delta_est' = _b[/delta]
    if "`factor'" == "common" {
        scalar `lambda_T' = _b[/lambda]
        scalar `lambda_Y' = `lambda_T'
        ereturn scalar lambda = `lambda_T'
    }
    else {
        scalar `lambda_T' = _b[/lambda_T]
        scalar `lambda_Y' = _b[/lambda_Y]
    }
    scalar `hr' = exp(`delta_est')

    ereturn scalar delta = `delta_est'
    ereturn scalar lambda_T = `lambda_T'
    ereturn scalar lambda_Y = `lambda_Y'
    ereturn scalar sigma_T = `lambda_T'
    ereturn scalar sigma_Y = `lambda_Y'
    ereturn scalar sigma_P = `lambda_T'
    ereturn scalar sigma_D = `lambda_Y'
    ereturn scalar hr = `hr'

    * Compute mixture probabilities
    tempname sum_exp_eta pi
    scalar `sum_exp_eta' = 1
    forvalues j = 2/`k' {
        scalar `sum_exp_eta' = `sum_exp_eta' + exp(_b[/eta_`j'])
    }

    matrix `pi' = J(1, `k', .)
    matrix `pi'[1, 1] = 1 / `sum_exp_eta'
    forvalues j = 2/`k' {
        matrix `pi'[1, `j'] = exp(_b[/eta_`j']) / `sum_exp_eta'
    }
    ereturn matrix pi = `pi'

    * Compute mass points. v_1 = 0 and v_2 = 1 are fixed by the
    * normalization; v_3..v_K are read from the coefficient vector.
    tempname v
    matrix `v' = J(1, `k', .)
    matrix `v'[1, 1] = 0
    if `k' >= 2 {
        matrix `v'[1, 2] = 1
    }
    forvalues j = 3/`k' {
        matrix `v'[1, `j'] = _b[/v_`j']
    }
    ereturn matrix v = `v'

    * Compute CI for hazard ratio
    local se_delta = _se[/delta]
    local z = invnormal(1 - (1 - `level'/100)/2)
    local hr_lo = exp(`delta_est' - `z' * `se_delta')
    local hr_hi = exp(`delta_est' + `z' * `se_delta')
    ereturn scalar hr_ci_lo = `hr_lo'
    ereturn scalar hr_ci_hi = `hr_hi'
    ereturn scalar level = `level'

    * Display results
    Display, level(`level')

    } // end capture noisily

    local rc = _rc

    * Guaranteed cleanup on all exit paths (normal, break, error)
    capture mata: _hsj_cleanup()
    capture macro drop HSJ_K HSJ_id HSJ_treat_event HSJ_outcome_event HSJ_treat
    capture macro drop HSJ_riskset HSJ_factor HSJ_nolog HSJ_treat_vars_mata HSJ_outcome_vars_mata
    capture matrix drop __hsj_start __hsj_b __hsj_V __hsj_g __hsj_best_b __hsj_best_V __hsj_best_g
    capture scalar drop __hsj_ll __hsj_has_result __hsj_converged __hsj_ic

    * Only propagate rc when the substantive estimation didn't complete.
    * If `e(cmd)` was set to "hsmixture_joint", `ereturn post` and the
    * identifying `ereturn local cmd` both succeeded, so the caller has
    * valid e(b), e(V), e(hr), e(delta) and friends. A non-zero rc from
    * a post-post bookkeeping step (gradient confirm, scalar drop, etc.)
    * should not be propagated as an estimation failure.
    if `rc' & "`e(cmd)'" != "hsmixture_joint" {
        exit `rc'
    }
    * Estimation succeeded — explicit clean exit so a non-zero rc lingering
    * from a `capture` cleanup above does not become the program's return rc.
    exit 0
end

program Display
    syntax [, Level(cilevel)]

    local K = e(K)
    local level_val = cond("`level'" != "", `level', e(level))

    display _n as txt "{hline 78}"
    display as txt "Joint Timing-of-Events Model with Heckman-Singer Heterogeneity"
    display as txt "{hline 78}"
    display as txt "Number of obs" _col(50) "=" _col(55) %12.0fc e(N)
    display as txt "K mass points" _col(50) "=" _col(55) %12.0f e(K)
    display as txt "Log likelihood" _col(50) "=" _col(55) %12.4f e(ll)
    if "`e(factor)'" == "common" {
        display as txt "Factor structure" _col(50) "=" _col(55) ///
            "       common (one shared loading)"
    }
    else {
        display as txt "Factor structure" _col(50) "=" _col(55) ///
            "     separate (lambda_T, lambda_Y)"
    }
    if "`e(riskset_var)'" != "" {
        display as txt "Treatment risk set" _col(50) "=" _col(55) "`e(riskset_var)' == 1"
    }
    display as txt "{hline 78}"

    * Strict-convergence warning. e(converged) = 1 only when BFGS converged AND
    * relative gradient |grad|/(1+|LL|) < 1e-5 AND V is positive definite. If
    * any of those failed, the reported point estimate may not be the MLE and
    * the SEs may not be trustworthy. Surfacing this in the main output, not
    * just postestimation.
    if e(converged) == 0 {
        display _n as err "{bf:WARNING: strict convergence not achieved.}"
        if e(rel_grad) != . & e(rel_grad) > 1e-5 {
            display as err "  Relative gradient |grad|/(1+|LL|) = " ///
                %10.2e e(rel_grad) " (threshold 1e-5; |grad| = " ///
                %10.2e e(grad_norm) ")."
        }
        if e(v_pd) == 0 {
            display as err "  Variance matrix is not positive definite."
        }
        if e(converged_bfgs) == 0 {
            display as err "  BFGS did not reach its own convergence criterion."
        }
        display as err "  The reported HR and CI may be unreliable. Run "
        display as err "  hsmixture_joint_postestimation, all  for details."
        display as err "{hline 78}"
    }

    * Display coefficient table
    ereturn display, level(`level_val')

    * Key result: Treatment effect
    local z = invnormal(1 - (1 - `level_val'/100)/2)
    local se_delta = _se[/delta]
    local hr_lo = exp(e(delta) - `z' * `se_delta')
    local hr_hi = exp(e(delta) + `z' * `se_delta')

    display _n as txt "{hline 78}"
    display as txt "{bf:Treatment Effect (Equation 1 -> Equation 2)}"
    display as txt "{hline 78}"
    display as txt "delta (log hazard ratio)" _col(30) "=" _col(35) %9.4f e(delta) ///
        _col(50) "SE = " %7.4f `se_delta'
    display as txt "Hazard Ratio" _col(30) "=" _col(35) %9.2f e(hr)
    display as txt "`level_val'% CI" _col(30) "=" _col(35) ///
        "[" %5.2f `hr_lo' ", " %5.2f `hr_hi' "]"

    * Factor loadings (signed under v_2=1 normalization)
    display _n as txt "Factor Loadings (signed; v_1=0, v_2=1 normalization):"
    if "`e(factor)'" == "common" {
        display as txt "  lambda (shared)" _col(30) "=" _col(35) %9.4f e(lambda_T)
    }
    else {
        display as txt "  lambda_T (treatment)" _col(30) "=" _col(35) %9.4f e(lambda_T)
        display as txt "  lambda_Y (outcome)" _col(30) "=" _col(35) %9.4f e(lambda_Y)
    }

    * Mixture distribution
    display _n as txt "Mixture Distribution:"
    display as txt "  Type" _col(15) "Prob (pi)" _col(30) "Mass Pt (v)" ///
        _col(50) "lambda_T*v" _col(65) "lambda_Y*v"
    display as txt "  {hline 70}"
    forvalues k = 1/`K' {
        local pi_k = e(pi)[1, `k']
        local v_k = e(v)[1, `k']
        local lT_v = e(lambda_T) * `v_k'
        local lY_v = e(lambda_Y) * `v_k'
        display as txt "  `k'" _col(15) %7.4f `pi_k' _col(30) %9.4f `v_k' ///
            _col(50) %9.4f `lT_v' _col(65) %9.4f `lY_v'
    }

    * Model fit statistics. BIC denominator is e(N_persons) because the IID
    * unit in this mixture model is the person; using person-period N would
    * inflate the penalty. e(N) is preserved for Stata convention.
    local ll = e(ll)
    local k_params = e(k)
    local N_p = e(N_persons)
    local aic = -2 * `ll' + 2 * `k_params'
    local bic = -2 * `ll' + `k_params' * ln(`N_p')

    display _n as txt "Model Fit:"
    display as txt "  AIC" _col(30) "=" _col(35) %12.2f `aic'
    display as txt "  BIC" _col(30) "=" _col(35) %12.2f `bic' ///
        "  (N = " %7.0fc `N_p' " persons)"
    display as txt "  Parameters" _col(30) "=" _col(35) %12.0f `k_params'

    display as txt "{hline 78}"
end
