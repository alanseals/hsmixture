*! version 2.3.2  30jun2026
*! Heckman-Singer Mixture Model - Single Equation
*! Discrete-time hazard model with unobserved heterogeneity
*! Mata-accelerated likelihood evaluator
*!
*! Authors: Jonghoon Park and R. Alan Seals (Auburn University)
*!
*! Reference: Heckman, J.J. and B. Singer. 1984. A method for minimizing
*!   the impact of distributional assumptions in econometric models for
*!   duration data. Econometrica 52(2): 271-320.

program hsmixture, eclass sortpreserve
    version 14

    if replay() {
        if "`e(cmd)'" != "hsmixture" {
            error 301
        }
        Display `0'
        exit
    }

    syntax varlist(min=1) [if] [in], ///
        ID(varname) ///
        [K(integer 2) ///
         FROM(name) ///
         ITERate(integer 200) ///
         noLOG ///
         Level(cilevel)]

    * Parse dependent and independent variables
    gettoken depvar indepvars : varlist

    * Mark sample
    marksample touse
    markout `touse' `id' `depvar'

    * Expand wildcards in indepvars and exclude obs with missing covariates
    * BEFORE the data-contract assertions. If a row has missing covariates it
    * will not enter the estimation sample, so it should not be checked by the
    * 0/1 / single-event assertions either.
    if "`indepvars'" != "" {
        unab indep_exp : `indepvars'
        markout `touse' `indep_exp'
    }
    else {
        * Intercept-only model: no covariates supplied. The Mata design always
        * appends a constant column, so an empty covariate list yields a valid
        * constant-only hazard. indepvars is documented as optional; guarding
        * the unab here keeps that promise (unab errors on an empty varlist).
        local indep_exp ""
    }

    * Validate K
    if `k' < 2 {
        display as error "k() must be at least 2"
        exit 198
    }

    * Stage 3: data-contract validation. Errors on violations.
    capture assert inlist(`depvar', 0, 1) if `touse'
    if _rc {
        display as error ///
            "outcome variable `depvar' must be 0/1 in the estimation sample"
        exit 198
    }

    * At most one event per id (absorbing outcome assumption).
    tempvar __hs_evcount
    quietly bysort `id': egen `__hs_evcount' = total(`depvar') if `touse'
    capture assert `__hs_evcount' <= 1 if `touse'
    if _rc {
        display as error ///
            "outcome variable `depvar' has more than one event for some id"
        display as error ///
            "  hsmixture requires an absorbing outcome (each person experiences the event at most once)"
        exit 198
    }
    drop `__hs_evcount'

    * Compile Mata functions (silently skips if already loaded)
    capture findfile hsmixture_mata.do
    if _rc {
        display as error "hsmixture_mata.do not found on adopath"
        exit 601
    }
    capture quietly do "`r(fn)'"

    * Set globals for likelihood evaluator
    global HS_K = `k'
    global HS_id "`id'"
    global HS_depvar "`depvar'"
    if "`log'" == "nolog" {
        global HS_nolog "1"
    }
    else {
        global HS_nolog ""
    }

    * Wrap estimation in capture for guaranteed cleanup
    capture noisily {

    * Fit initial GLM to get starting values and expand factor variables
    quietly glm `depvar' `indepvars' if `touse', family(binomial) link(cloglog) nolog

    tempname b_glm
    matrix `b_glm' = e(b)

    * `indep_exp' was expanded above (before the data-contract asserts) so
    * touse already excludes obs with missing covariates.
    local xvars_mata "`indep_exp'"
    global HS_xvars_mata "`xvars_mata'"

    * Ensure data is sorted by id for Mata panelsetup
    sort `id'

    * Build starting vector from GLM coefficients.
    * Layout under v_2=1 normalization:
    *   [beta (k_x) | lambda | v_3..v_K | eta_2..eta_K]
    * Lambda is a signed real-line factor loading. v_2 is fixed at 1 in the
    * Mata likelihood and is not a free parameter.
    tempname b0
    matrix `b0' = `b_glm'

    * lambda: start at 0.5 (signed, no log transform)
    matrix `b0' = `b0', 0.5

    * v_3, ..., v_K: starting at 2.0 with unit increments. For K=2 this loop
    * is empty (v_2=1 is the second mass point and is fixed).
    forvalues j = 3/`k' {
        local v_init = `j' - 1
        matrix `b0' = `b0', `v_init'
    }

    * eta_2, ..., eta_K: start at 0 (equal probabilities)
    forvalues j = 2/`k' {
        matrix `b0' = `b0', 0
    }

    * Use user-supplied starting values if provided
    if "`from'" != "" {
        matrix `b0' = `from'
    }

    * Display header
    if "`log'" != "nolog" {
        display _n as txt "Heckman-Singer Mixture Model (K = `k' mass points)"
        display as txt "Single-equation discrete-time hazard"
        display as txt "Mass-point normalization: v_1 = 0, v_2 = 1"
        display as txt "{hline 60}"
    }

    * Initialize Mata cache
    capture mata: _hs_cleanup()

    * Run Mata optimizer
    matrix __hs_start = `b0'
    scalar __hs_has_result = 0
    scalar __hs_converged = 0

    capture noisily mata: _hs_run_optimize("__hs_start", `iterate')

    * Check results
    capture confirm scalar __hs_has_result
    if _rc != 0 | scalar(__hs_has_result) != 1 {
        display as error "Optimization failed to converge"
        exit 430
    }

    local best_ll = scalar(__hs_ll)
    local best_converged = scalar(__hs_converged)
    local best_ic = scalar(__hs_ic)
    matrix __hs_best_b = __hs_b
    capture confirm matrix __hs_V
    if _rc == 0 {
        matrix __hs_best_V = __hs_V
    }
    else {
        * Hessian inversion failed. Post a finite scaffold V so `ereturn post`
        * succeeds (Stata rejects all-missing V). Callers should gate on
        * e(converged)==1 before trusting SEs.
        local _np = colsof(__hs_best_b)
        matrix __hs_best_V = I(`_np') * 1e-20
    }
    capture confirm matrix __hs_g
    if _rc == 0 {
        matrix __hs_best_g = __hs_g
    }

    capture matrix drop __hs_start __hs_b __hs_V __hs_g
    capture scalar drop __hs_ll __hs_has_result __hs_converged __hs_ic

    * ----------------------------------------------------------------
    * Post results to ereturn
    * ----------------------------------------------------------------

    * Build column names. Under v_2=1 normalization, the free mass-point
    * parameters are v_3..v_K (v_1=0 and v_2=1 are fixed).
    local colnames ""
    foreach v of local indep_exp {
        local colnames "`colnames' xb:`v'"
    }
    local colnames "`colnames' xb:_cons"
    local colnames "`colnames' lambda:_cons"
    forvalues j = 3/`k' {
        local colnames "`colnames' v_`j':_cons"
    }
    forvalues j = 2/`k' {
        local colnames "`colnames' eta_`j':_cons"
    }

    matrix colnames __hs_best_b = `colnames'
    matrix colnames __hs_best_V = `colnames'
    matrix rownames __hs_best_V = `colnames'

    quietly count if `touse'
    local N_obs = r(N)

    * Person-level count for AIC/BIC. The IID unit in this mixture model is
    * the person; the latent type is integrated out at the person level (see
    * Mata _hs_compute_ll). Using person-period N would inflate the BIC
    * penalty. e(N) follows Stata's row-count convention; e(N_persons) is
    * the statistically appropriate denominator for cross-K comparisons.
    tempvar __hs_pid_tag
    quietly egen byte `__hs_pid_tag' = tag(`id') if `touse'
    quietly count if `__hs_pid_tag' == 1
    local N_persons = r(N)
    drop `__hs_pid_tag'

    ereturn post __hs_best_b __hs_best_V, obs(`N_obs')
    ereturn scalar ll = `best_ll'
    ereturn scalar N_persons = `N_persons'

    * Clean up Mata cache
    capture mata: _hs_cleanup()

    * Strict-convergence diagnostics (see hsmixture_joint.ado for rationale).
    * Uses relative gradient |grad|/(1+|LL|) < 1e-5.
    local strict_converged = `best_converged'
    local grad_norm = .
    local rel_grad = .
    local v_pd = 0

    capture confirm matrix __hs_best_g
    if _rc == 0 {
        capture mata: st_numscalar("__hs_gn", norm(st_matrix("__hs_best_g")))
        if _rc == 0 {
            local grad_norm = scalar(__hs_gn)
            local rel_grad = `grad_norm' / (1 + abs(`best_ll'))
            if `rel_grad' > 1e-5 local strict_converged = 0
        }
        else {
            local strict_converged = 0
        }
        capture scalar drop __hs_gn
    }
    else {
        local strict_converged = 0
    }

    capture mata: st_numscalar("__hs_meig", min(Re(eigenvalues(st_matrix("e(V)")))))
    if _rc == 0 {
        if !missing(scalar(__hs_meig)) & scalar(__hs_meig) > 1e-8 local v_pd = 1
    }
    if !`v_pd' local strict_converged = 0
    capture scalar drop __hs_meig

    * Store estimation results
    ereturn local cmd "hsmixture"
    ereturn local cmdline "hsmixture `0'"
    ereturn local depvar "`depvar'"
    ereturn local idvar "`id'"
    ereturn scalar K = `k'
    * Force estimates stats / IC machinery to use the *design* parameter
    * count rather than rank(V). The single-equation eta_2 row commonly
    * has a missing SE when pi_2 is near a boundary; without this override
    * `estimates stats` would report df = colsof(e(b)) - 1 and produce an
    * AIC/BIC that disagrees with the package's own display.
    ereturn scalar rank = colsof(e(b))
    ereturn scalar df_m = colsof(e(b))
    ereturn scalar converged = `strict_converged'
    ereturn scalar converged_bfgs = `best_converged'
    ereturn scalar grad_norm = `grad_norm'
    ereturn scalar rel_grad = `rel_grad'
    ereturn scalar v_pd = `v_pd'
    ereturn scalar ic = `best_ic'

    * Post gradient at optimum (consumed by postestimation diagnostics)
    capture confirm matrix __hs_best_g
    if _rc == 0 {
        matrix colnames __hs_best_g = `colnames'
        ereturn matrix gradient = __hs_best_g
    }

    * Compute transformed parameters. lambda is a signed real-line factor
    * loading. The legacy `e(sigma)` is preserved as an alias.
    tempname lambda
    scalar `lambda' = _b[/lambda]
    ereturn scalar lambda = `lambda'
    ereturn scalar sigma = `lambda'

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
    ereturn scalar level = `level'

    * Display results
    Display, level(`level')

    } // end capture noisily

    local rc = _rc

    * Guaranteed cleanup
    capture mata: _hs_cleanup()
    capture macro drop HS_K HS_id HS_depvar HS_nolog HS_xvars_mata
    capture matrix drop __hs_start __hs_b __hs_V __hs_g __hs_best_b __hs_best_V __hs_best_g
    capture scalar drop __hs_ll __hs_has_result __hs_converged __hs_ic

    if `rc' {
        exit `rc'
    }
end

program Display
    syntax [, Level(cilevel)]

    local level_val = cond("`level'" != "", `level', e(level))

    display _n as txt "Heckman-Singer Mixture Model" ///
        _col(50) "Number of obs" _col(67) "=" _col(69) %10.0fc e(N)
    display as txt "K = " e(K) " mass points" ///
        _col(50) "Log likelihood" _col(67) "=" _col(69) %10.4f e(ll)
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
        display as err "  Estimates may be unreliable."
        display as err "{hline 78}"
    }

    * Display coefficients
    ereturn display, level(`level_val')

    * Display transformed parameters
    display _n as txt "Transformed Parameters (signed; v_1=0, v_2=1 normalization):"
    display as txt "{hline 60}"
    display as txt "lambda" _col(20) "=" _col(25) %8.4f e(lambda)

    display _n as txt "Mixture Probabilities:"
    local K = e(K)
    forvalues k = 1/`K' {
        display as txt "  pi_`k'" _col(20) "=" _col(25) %8.4f e(pi)[1, `k']
    }

    display _n as txt "Mass Points:"
    forvalues k = 1/`K' {
        display as txt "  v_`k'" _col(20) "=" _col(25) %8.4f e(v)[1, `k']
    }

    * AIC and BIC. BIC denominator is e(N_persons) because the IID unit in
    * this mixture model is the person; using person-period N would inflate
    * the penalty. e(N) is preserved for Stata convention.
    local ll = e(ll)
    local k_params = colsof(e(b))
    local N_p = e(N_persons)
    local aic = -2 * `ll' + 2 * `k_params'
    local bic = -2 * `ll' + `k_params' * ln(`N_p')
    display _n as txt "AIC" _col(20) "=" _col(25) %12.2f `aic'
    display as txt "BIC" _col(20) "=" _col(25) %12.2f `bic' ///
        "  (N = " %7.0fc `N_p' " persons)"

end
