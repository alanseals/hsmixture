*! version 2.3.1  07may2026
*! Bivariate Heterogeneity Joint Timing-of-Events Model
*! Two-dimensional discrete-time hazard with SEPARATE latent types for
*! treatment (v_T) and outcome (v_Y) on a 2x2 joint grid.
*!
*! Authors: Jonghoon Park and R. Alan Seals (Auburn University)
*!
*! Relation to hsmixture_joint:
*!   hsmixture_joint, factor(common):    one shared loading; per-type shifts
*!     (lambda*v_k, lambda*v_k) lie on the 45-degree line. Manuscript-style
*!     one-factor MPH.
*!   hsmixture_joint, factor(separate):  two free signed loadings; per-type
*!     shifts (lambda_T*v_k, lambda_Y*v_k) lie on a 1-D locus through the
*!     origin (sign of locus is data-determined, opposite signs admissible).
*!   hsmixture_bivariate (this command): 2x2 grid of corner shifts (0,0),
*!     (0,v_Y2), (v_T2,0), (v_T2,v_Y2) with full joint probability matrix.
*!     Strictly nests both joint variants when the off-diagonal probabilities
*!     are unconstrained.

program hsmixture_bivariate, eclass sortpreserve
    version 14

    if replay() {
        if "`e(cmd)'" != "hsmixture_bivariate" {
            error 301
        }
        Display `0'
        exit
    }

    * Parse the compound syntax
    gettoken eq1 0 : 0, parse("(") match(paren)
    gettoken eq2 0 : 0, parse("(") match(paren)

    * Parse equation 1 (treatment process)
    gettoken treat_dep eq1_rest : eq1, parse("=")
    gettoken eq_sign treat_indep : eq1_rest, parse("=")
    local treat_dep = strtrim("`treat_dep'")
    local treat_indep = strtrim("`treat_indep'")

    * Parse equation 2 (outcome process with treatment indicator)
    local treat_var ""
    if regexm("`eq2'", "treat\(([^\)]+)\)") {
        local treat_var = regexs(1)
        local eq2 = regexr("`eq2'", ",[ ]*treat\([^\)]+\)", "")
    }

    gettoken outcome_dep eq2_rest : eq2, parse("=")
    gettoken eq_sign outcome_indep : eq2_rest, parse("=")
    local outcome_dep = strtrim("`outcome_dep'")
    local outcome_indep = strtrim("`outcome_indep'")

    * Parse remaining options. PRisk() is a deprecated alias for RISKset()
    * retained for back-compat with v2.0.0 production scripts.
    *
    * DIFFicult, TRace, GRADient, HESSian, TECHnique(), TOLerance(),
    * LTOLerance(), NRTOLerance() are legacy `ml model d0` options from
    * v2.0.0. They are silently accepted but ignored by the Mata optimize()
    * implementation.
    syntax [if] [in], ///
        ID(varname) ///
        [RISKset(varname) ///
         PRisk(varname) ///
         FROM(name) ///
         ITERate(integer 200) ///
         noLOG ///
         NSTarts(integer 6) ///
         Level(cilevel) ///
         DIFFicult ///
         TRace ///
         GRADient ///
         HESSian ///
         TECHnique(string) ///
         TOLerance(real 1e-6) ///
         LTOLerance(real 1e-7) ///
         NRTOLerance(real 1e-5)]

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

    * Validate
    if "`treat_var'" == "" {
        display as error "treat() required in second equation"
        exit 198
    }

    * Expand wildcards and exclude obs with missing covariates from touse
    * BEFORE the data-contract assertions. If a row has missing covariates it
    * will not enter the estimation sample, so it must not be checked by the
    * 0/1 / single-event / riskset assertions either.
    unab treat_indep_exp : `treat_indep'
    unab outcome_indep_exp : `outcome_indep'
    markout `touse' `treat_indep_exp' `outcome_indep_exp'

    * Stage 3: data-contract validation. Errors on violations.
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

    * Stage 3 (extended): duration-data structural validation.
    tempvar __hsb_evcount
    quietly bysort `id': egen `__hsb_evcount' = total(`treat_dep') if `touse'
    capture assert `__hsb_evcount' <= 1 if `touse'
    if _rc {
        display as error ///
            "treatment event variable `treat_dep' has more than one event for some id"
        display as error ///
            "  hsmixture_bivariate requires one-time treatment timing"
        exit 198
    }
    drop `__hsb_evcount'

    tempvar __hsb_evcount
    quietly bysort `id': egen `__hsb_evcount' = total(`outcome_dep') if `touse'
    capture assert `__hsb_evcount' <= 1 if `touse'
    if _rc {
        display as error ///
            "outcome event variable `outcome_dep' has more than one event for some id"
        display as error ///
            "  hsmixture_bivariate requires an absorbing outcome"
        exit 198
    }
    drop `__hsb_evcount'

    if "`riskset'" != "" {
        capture assert `riskset' == 1 if `treat_dep' == 1 & `touse'
        if _rc {
            display as error ///
                "treatment events occur outside the treatment risk set"
            exit 198
        }
    }

    * Compile Mata functions (silently skips if already loaded)
    capture findfile hsmixture_bivariate_mata.do
    if _rc {
        display as error "hsmixture_bivariate_mata.do not found on adopath"
        exit 601
    }
    capture quietly do "`r(fn)'"

    * Set globals for likelihood evaluator
    global HSB_id "`id'"
    global HSB_treat_event "`treat_dep'"
    global HSB_outcome_event "`outcome_dep'"
    global HSB_treat "`treat_var'"
    global HSB_riskset "`riskset'"
    if "`log'" == "nolog" {
        global HSB_nolog "1"
    }
    else {
        global HSB_nolog ""
    }

    * Wrap estimation in capture for guaranteed cleanup on break/error
    capture noisily {

    * ----------------------------------------------------------------
    * Initial values from separate GLMs (also expands factor variables)
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
    global HSB_treat_vars_mata "`treat_vars_mata'"
    global HSB_outcome_vars_mata "`outcome_vars_mata'"

    * Delta from cloglog with treatment
    quietly glm `outcome_dep' `treat_var' `outcome_indep' if `touse', ///
        family(binomial) link(cloglog) nolog
    local delta_init = _b[`treat_var']

    * ----------------------------------------------------------------
    * Starting value configurations for v_T2 and v_Y2
    * ----------------------------------------------------------------

    local sv_vT2_1 = 0.5
    local sv_vY2_1 = 0.5
    local sv_vT2_2 = 1.0
    local sv_vY2_2 = 1.0
    local sv_vT2_3 = 1.0
    local sv_vY2_3 = 1.5
    local sv_vT2_4 = 1.5
    local sv_vY2_4 = 1.0
    local sv_vT2_5 = 1.5
    local sv_vY2_5 = 2.0
    local sv_vT2_6 = 2.0
    local sv_vY2_6 = 1.5

    * Display header
    if "`log'" != "nolog" {
        display _n as txt "Bivariate Heterogeneity Joint Timing-of-Events Model"
        display as txt "2x2 joint grid: K_T=2, K_Y=2 (4 joint types)"
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
    capture mata: _hsb_cleanup()

    * ----------------------------------------------------------------
    * Run with multiple starting values using Mata optimize()
    * ----------------------------------------------------------------

    local best_ll = .
    local best_converged = 0
    local nstarts_use = min(`nstarts', 6)

    forvalues s = 1/`nstarts_use' {
        * Build starting vector for this configuration
        tempname b0_`s'
        matrix `b0_`s'' = `b_treat', `b_outcome', `delta_init', ///
            `sv_vT2_`s'', `sv_vY2_`s'', 0, 0, 0

        * Use user-supplied starting values if provided (only first config)
        if "`from'" != "" & `s' == 1 {
            matrix `b0_`s'' = `from'
        }

        if "`log'" != "nolog" {
            display _n as txt ///
                "--- Starting configuration `s'/`nstarts_use': " ///
                "v_T2=`sv_vT2_`s'', v_Y2=`sv_vY2_`s'' ---"
        }

        * Post starting values and run Mata optimizer
        matrix __hsb_start = `b0_`s''
        scalar __hsb_has_result = 0
        scalar __hsb_converged = 0

        capture noisily mata: _hsb_run_optimize("__hsb_start", `iterate')

        * Check if Mata posted a result (finite LL)
        capture confirm scalar __hsb_has_result
        if _rc == 0 & scalar(__hsb_has_result) == 1 {
            local ll_s = scalar(__hsb_ll)
            if "`log'" != "nolog" {
                local _pos_d = colsof(`b_treat') + colsof(`b_outcome') + 1
                tempname _btemp
                matrix `_btemp' = __hsb_b
                display as txt "  Log-likelihood: " %10.2f `ll_s' ///
                    "  HR: " %5.3f exp(`_btemp'[1, `_pos_d'])
            }

            if `ll_s' > `best_ll' | `best_ll' == . {
                local best_ll = `ll_s'
                local best_start = `s'
                local best_converged = scalar(__hsb_converged)
                local best_ic = scalar(__hsb_ic)
                matrix __hsb_best_b = __hsb_b
                capture confirm matrix __hsb_V
                if _rc == 0 {
                    matrix __hsb_best_V = __hsb_V
                }
                else {
                    * Hessian inversion failed. Post a finite scaffold V so
                    * `ereturn post` succeeds. Stata rejects all-missing V
                    * with "matrix has missing values". Callers should gate
                    * on e(converged)==1 (V positive definite) before trusting
                    * SEs.
                    local _np = colsof(__hsb_best_b)
                    matrix __hsb_best_V = I(`_np') * 1e-20
                }
                capture confirm matrix __hsb_g
                if _rc == 0 {
                    matrix __hsb_best_g = __hsb_g
                }
                else {
                    capture matrix drop __hsb_best_g
                }
            }
        }
        else {
            if "`log'" != "nolog" {
                display as txt "  Configuration `s' failed to converge"
            }
        }

        * Clean up per-iteration temporaries
        capture matrix drop __hsb_start __hsb_b __hsb_V __hsb_g
        capture scalar drop __hsb_ll __hsb_has_result __hsb_converged __hsb_ic
    }

    * ----------------------------------------------------------------
    * Check that at least one configuration produced results
    * ----------------------------------------------------------------

    if `best_ll' == . {
        display as error "All starting configurations failed to converge"
        exit 430
    }

    if "`log'" != "nolog" {
        display _n as txt "Best result from starting configuration `best_start'" ///
            " (log-lik: " %10.2f `best_ll' ")"
    }

    * ----------------------------------------------------------------
    * Post results to ereturn
    * ----------------------------------------------------------------

    * Equation labels use the actual dependent-variable names so that
    * post-estimation `_b[depvar:varname]` calls match what the user wrote.
    local colnames ""
    foreach v of local treat_indep_exp {
        local colnames "`colnames' `treat_dep':`v'"
    }
    local colnames "`colnames' `treat_dep':_cons"
    foreach v of local outcome_indep_exp {
        local colnames "`colnames' `outcome_dep':`v'"
    }
    local colnames "`colnames' `outcome_dep':_cons"
    local colnames "`colnames' delta:_cons v_T2:_cons v_Y2:_cons"
    local colnames "`colnames' logit_pi_12:_cons logit_pi_21:_cons logit_pi_22:_cons"

    matrix colnames __hsb_best_b = `colnames'
    matrix colnames __hsb_best_V = `colnames'
    matrix rownames __hsb_best_V = `colnames'

    quietly count if `touse'
    local N_obs = r(N)

    * Person-level count for AIC/BIC. The IID unit in this mixture model is
    * the person; the latent type is integrated out at the person level (see
    * Mata _hsb_compute_ll). Using person-period N would inflate the BIC
    * penalty. e(N) follows Stata's row-count convention; e(N_persons) is
    * the statistically appropriate denominator for cross-K comparisons.
    tempvar __hsb_pid_tag
    quietly egen byte `__hsb_pid_tag' = tag(`id') if `touse'
    quietly count if `__hsb_pid_tag' == 1
    local N_persons = r(N)
    drop `__hsb_pid_tag'

    ereturn post __hsb_best_b __hsb_best_V, obs(`N_obs')
    ereturn scalar ll = `best_ll'
    ereturn scalar N_persons = `N_persons'

    * Clean up Mata cache
    capture mata: _hsb_cleanup()

    * Strict-convergence diagnostics. See hsmixture_joint.ado for rationale.
    * Uses relative gradient |grad|/(1+|LL|) < 1e-5.
    local strict_converged = `best_converged'
    local grad_norm = .
    local rel_grad = .
    local v_pd = 0

    capture confirm matrix __hsb_best_g
    if _rc == 0 {
        capture mata: st_numscalar("__hsb_gn", norm(st_matrix("__hsb_best_g")))
        if _rc == 0 {
            local grad_norm = scalar(__hsb_gn)
            local rel_grad = `grad_norm' / (1 + abs(`best_ll'))
            if `rel_grad' > 1e-5 local strict_converged = 0
        }
        else {
            local strict_converged = 0
        }
        capture scalar drop __hsb_gn
    }
    else {
        local strict_converged = 0
    }

    capture mata: st_numscalar("__hsb_meig", min(Re(eigenvalues(st_matrix("e(V)")))))
    if _rc == 0 {
        if scalar(__hsb_meig) > 1e-8 local v_pd = 1
    }
    if !`v_pd' local strict_converged = 0
    capture scalar drop __hsb_meig

    * Store results
    ereturn local cmd "hsmixture_bivariate"
    ereturn local cmdline "hsmixture_bivariate `0'"
    ereturn local treat_depvar "`treat_dep'"
    ereturn local outcome_depvar "`outcome_dep'"
    ereturn local treat_var "`treat_var'"
    ereturn local idvar "`id'"
    ereturn local riskset_var "`riskset'"
    ereturn scalar n_starts = `nstarts_use'
    ereturn scalar best_start = `best_start'
    ereturn scalar k = colsof(e(b))
    * Force estimates stats / IC machinery to use the *design* parameter
    * count rather than rank(V). When the Hessian is rank-deficient
    * (typical for the bivariate softmax when one off-diagonal mass cell
    * is empirically near zero), Stata's default reduces df below the
    * actual number of free parameters and produces an AIC/BIC that
    * disagrees with the package's own display.
    ereturn scalar rank = colsof(e(b))
    ereturn scalar df_m = colsof(e(b))
    ereturn scalar converged = `strict_converged'
    ereturn scalar converged_bfgs = `best_converged'
    ereturn scalar grad_norm = `grad_norm'
    ereturn scalar rel_grad = `rel_grad'
    ereturn scalar v_pd = `v_pd'
    ereturn scalar ic = `best_ic'

    * Post gradient at best-start optimum (consumed by
    * hsmixture_joint_postestimation, convergence)
    capture confirm matrix __hsb_best_g
    if _rc == 0 {
        matrix colnames __hsb_best_g = `colnames'
        ereturn matrix gradient = __hsb_best_g
    }

    * Treatment effect
    tempname delta_est hr
    scalar `delta_est' = _b[/delta]
    scalar `hr' = exp(`delta_est')
    ereturn scalar delta = `delta_est'
    ereturn scalar hr = `hr'

    * CI for hazard ratio
    local se_delta = _se[/delta]
    local z = invnormal(1 - (1 - `level'/100)/2)
    ereturn scalar hr_ci_lo = exp(`delta_est' - `z' * `se_delta')
    ereturn scalar hr_ci_hi = exp(`delta_est' + `z' * `se_delta')
    ereturn scalar se_delta = `se_delta'
    ereturn scalar level = `level'

    * Mass points
    ereturn scalar v_T2 = _b[/v_T2]
    ereturn scalar v_Y2 = _b[/v_Y2]

    * Joint probabilities via softmax
    tempname lp12 lp21 lp22 sum_e
    scalar `lp12' = _b[/logit_pi_12]
    scalar `lp21' = _b[/logit_pi_21]
    scalar `lp22' = _b[/logit_pi_22]
    scalar `sum_e' = 1 + exp(`lp12') + exp(`lp21') + exp(`lp22')

    tempname pi_mat
    matrix `pi_mat' = J(2, 2, .)
    matrix `pi_mat'[1, 1] = 1 / `sum_e'
    matrix `pi_mat'[1, 2] = exp(`lp12') / `sum_e'
    matrix `pi_mat'[2, 1] = exp(`lp21') / `sum_e'
    matrix `pi_mat'[2, 2] = exp(`lp22') / `sum_e'
    matrix rownames `pi_mat' = "v_T=0" "v_T=v_T2"
    matrix colnames `pi_mat' = "v_Y=0" "v_Y=v_Y2"
    ereturn matrix pi_joint = `pi_mat'

    * Implied correlation
    local pi11 = 1 / `sum_e'
    local pi12 = exp(`lp12') / `sum_e'
    local pi21 = exp(`lp21') / `sum_e'
    local pi22 = exp(`lp22') / `sum_e'

    local E_vT = (`pi21' + `pi22') * _b[/v_T2]
    local E_vY = (`pi12' + `pi22') * _b[/v_Y2]
    local Var_vT = (`pi21' + `pi22') * _b[/v_T2]^2 - `E_vT'^2
    local Var_vY = (`pi12' + `pi22') * _b[/v_Y2]^2 - `E_vY'^2
    local Cov_TY = `pi22' * _b[/v_T2] * _b[/v_Y2] - `E_vT' * `E_vY'

    if `Var_vT' > 0 & `Var_vY' > 0 {
        ereturn scalar rho = `Cov_TY' / sqrt(`Var_vT' * `Var_vY')
    }
    else {
        ereturn scalar rho = .
    }

    * Model fit -- compute from actual parameter count
    local n_treat_params = colsof(`b_treat')
    local n_outcome_params = colsof(`b_outcome')
    local n_total_params = `n_treat_params' + `n_outcome_params' + 6
    ereturn scalar n_params = `n_total_params'

    * Display results
    Display, level(`level')

    } // end capture noisily

    local rc = _rc

    * Guaranteed cleanup on all exit paths (normal, break, error)
    capture mata: _hsb_cleanup()
    capture macro drop HSB_id HSB_treat_event HSB_outcome_event HSB_treat HSB_riskset
    capture macro drop HSB_nolog HSB_treat_vars_mata HSB_outcome_vars_mata
    capture matrix drop __hsb_start __hsb_b __hsb_V __hsb_g __hsb_best_b __hsb_best_V __hsb_best_g
    capture scalar drop __hsb_ll __hsb_has_result __hsb_converged __hsb_ic

    * Only propagate rc when the substantive estimation didn't complete.
    * If `e(cmd)` was set to "hsmixture_bivariate", `ereturn post` and the
    * identifying `ereturn local cmd` both succeeded, so the caller has
    * valid e() results regardless of any post-post bookkeeping rc.
    if `rc' & "`e(cmd)'" != "hsmixture_bivariate" {
        exit `rc'
    }
    * Estimation succeeded — explicit clean exit so a non-zero rc lingering
    * from a `capture` cleanup above does not become the program's return rc.
    exit 0
end

program Display
    syntax [, Level(cilevel)]

    local level_val = cond("`level'" != "", `level', e(level))

    display _n as txt "{hline 78}"
    display as txt "Bivariate Heterogeneity Joint Timing-of-Events Model"
    display as txt "Two-dimensional unobserved heterogeneity (2x2 joint grid)"
    display as txt "{hline 78}"
    display as txt "Number of obs" _col(50) "=" _col(55) %12.0fc e(N)
    display as txt "Log likelihood" _col(50) "=" _col(55) %12.4f e(ll)
    display as txt "Starting configs tried" _col(50) "=" _col(55) %12.0f e(n_starts)
    display as txt "Best starting config" _col(50) "=" _col(55) %12.0f e(best_start)
    if "`e(riskset_var)'" != "" {
        display as txt "Treatment risk set" _col(50) "=" _col(55) "`e(riskset_var)' == 1"
    }
    display as txt "{hline 78}"

    * Strict-convergence warning. See hsmixture_joint.ado Display for rationale.
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
        display as err "  The reported HR and CI may be unreliable."
        display as err "{hline 78}"
    }

    * Display coefficient table
    ereturn display, level(`level_val')

    * Treatment effect
    local z = invnormal(1 - (1 - `level_val'/100)/2)
    local se_delta = e(se_delta)
    local hr_lo = exp(e(delta) - `z' * `se_delta')
    local hr_hi = exp(e(delta) + `z' * `se_delta')

    display _n as txt "{hline 78}"
    display as txt "{bf:Treatment Effect (Equation 1 -> Equation 2)}"
    display as txt "{hline 78}"
    display as txt "delta (log hazard ratio)" _col(30) "=" _col(35) %9.4f e(delta) ///
        _col(50) "SE = " %7.4f e(se_delta)
    display as txt "Hazard Ratio" _col(30) "=" _col(35) %9.2f e(hr)
    display as txt "`level_val'% CI" _col(30) "=" _col(35) ///
        "[" %5.2f `hr_lo' ", " %5.2f `hr_hi' "]"

    * Mass points
    display _n as txt "Bivariate Mass Points:"
    display as txt "  Treatment: v_T = {0, " %6.3f e(v_T2) "}"
    display as txt "  Outcome:   v_Y = {0, " %6.3f e(v_Y2) "}"

    * Joint probability matrix
    display _n as txt "Joint Probability Matrix (pi_{jk}):"
    matrix list e(pi_joint), format(%7.4f)

    * Implied correlation
    display _n as txt "Implied Correlation:"
    display as txt "  rho(v_T, v_Y)" _col(30) "=" _col(35) %7.4f e(rho)

    * Model fit. BIC denominator is e(N_persons) because the IID unit in
    * this mixture model is the person; using person-period N would inflate
    * the penalty. e(N) is preserved for Stata convention.
    local ll = e(ll)
    local np = e(n_params)
    local N_p = e(N_persons)
    local aic = -2 * `ll' + 2 * `np'
    local bic = -2 * `ll' + `np' * ln(`N_p')

    display _n as txt "Model Fit:"
    display as txt "  Parameters" _col(30) "=" _col(35) %12.0f `np'
    display as txt "  AIC" _col(30) "=" _col(35) %12.2f `aic'
    display as txt "  BIC" _col(30) "=" _col(35) %12.2f `bic' ///
        "  (N = " %7.0fc `N_p' " persons)"

    display as txt "{hline 78}"
end
