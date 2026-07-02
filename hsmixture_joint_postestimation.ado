*! version 2.3.3  02jul2026
*! Heckman-Singer Joint ToE Model - Postestimation Diagnostics
*! Convergence checks, mixture-distribution summary, and model comparison.
*! Adapts display under factor(common) when called after hsmixture_joint.
*!
*! Authors: Jonghoon Park and R. Alan Seals (Auburn University)

program hsmixture_joint_postestimation, rclass
    version 14

    syntax [, ///
        CONVergence ///      Display convergence diagnostics
        POSTerior ///        Display mixture distribution summary
        COMPare ///          Compare models (requires stored estimates)
        LRtest(name) ///     LR test vs specified model
        ALL ///              Display all diagnostics
        ]

    * Verify estimation context
    if "`e(cmd)'" != "hsmixture_joint" & "`e(cmd)'" != "hsmixture_bivariate" {
        display as error "hsmixture_joint_postestimation requires " ///
            "hsmixture_joint or hsmixture_bivariate estimates"
        exit 301
    }

    local is_bivariate = ("`e(cmd)'" == "hsmixture_bivariate")

    * Factor mode of the underlying joint fit. Empty for bivariate or for a
    * pre-factor()-option joint fit; treated as "separate" in either case.
    local factor_mode = "`e(factor)'"
    if "`factor_mode'" == "" local factor_mode "separate"

    if `is_bivariate' {
        local K = 2
    }
    else {
        local K = e(K)
    }
    local N = e(N)
    * Person count for IC denominators. Falls back to row count e(N) for
    * pre-v2.3.1 fits that did not store e(N_persons). New fits use the
    * person count (the IID unit for this mixture model).
    capture local N_p = e(N_persons)
    if "`N_p'" == "" | `N_p' == . local N_p = `N'
    local ll = e(ll)
    local k_params = e(k)

    * Default: show convergence if nothing specified
    if "`convergence'`posterior'`compare'`lrtest'`all'" == "" {
        local convergence "convergence"
    }

    if "`all'" != "" {
        local convergence "convergence"
        local posterior "posterior"
        local compare "compare"
    }

    *==========================================================================
    * CONVERGENCE DIAGNOSTICS
    *==========================================================================
    if "`convergence'" != "" {
        display _n as txt "{hline 78}"
        display as txt "{bf:CONVERGENCE DIAGNOSTICS}"
        display as txt "{hline 78}"

        * Check if converged
        capture local converged = e(converged)
        if _rc == 0 {
            if `converged' == 1 {
                display as txt "Convergence status:" _col(35) as res "Converged"
            }
            else {
                display as txt "Convergence status:" _col(35) as err "DID NOT CONVERGE"
                display as err "  WARNING: Results may be unreliable"
            }
        }
        else {
            display as txt "Convergence status:" _col(35) as txt "(not available)"
            local converged = .
        }

        capture display as txt "Iterations:" _col(35) as res e(ic)
        display as txt "Log likelihood:" _col(35) as res %12.4f `ll'

        * Gradient check. Reports both absolute and relative gradient norm;
        * the relative metric |grad|/(1+|LL|) is the strict-convergence
        * criterion (threshold 1e-5).
        local grad_norm = .
        local rel_grad = .
        tempname gradient
        capture matrix `gradient' = e(gradient)
        if _rc == 0 {
            mata: st_numscalar("r(grad_norm)", norm(st_matrix("`gradient'")))
            local grad_norm = r(grad_norm)
            local rel_grad = `grad_norm' / (1 + abs(`ll'))
            display as txt "Gradient norm:" _col(35) as res %12.6g `grad_norm'
            display as txt "  |grad|/(1+|LL|):" _col(35) as res %12.6g `rel_grad'
            if `rel_grad' != . & `rel_grad' > 1e-5 {
                display as err "  WARNING: relative gradient > 1e-5 suggests incomplete convergence"
            }
        }
        else {
            display as txt "Gradient norm:" _col(35) as txt "(not available)"
        }

        * Hessian check (negative definite at maximum)
        display _n as txt "Variance matrix check:"
        tempname V
        matrix `V' = e(V)

        * Check for missing/infinite values
        local vcov_ok = 1
        local nrows = rowsof(`V')
        forvalues i = 1/`nrows' {
            forvalues j = 1/`nrows' {
                if `V'[`i', `j'] == . | abs(`V'[`i', `j']) > 1e10 {
                    local vcov_ok = 0
                }
            }
        }

        if `vcov_ok' == 1 {
            * Match the estimation-time convergence gate. Each estimator posts a
            * scaffold V = I*1e-20 when Hessian inversion fails and rejects it
            * with a min-eigenvalue > 1e-8 test (setting e(v_pd)=0). A bare ">0"
            * test here would instead certify that placeholder as positive
            * definite. Prefer the estimation-time verdict e(v_pd) when present;
            * otherwise recompute with the same 1e-8 floor.
            capture local vpd = e(v_pd)
            if "`vpd'" != "" & "`vpd'" != "." {
                local pd_ok = (`vpd' == 1)
            }
            else {
                mata: st_numscalar("r(min_eigen)", min(Re(eigenvalues(st_matrix("`V'")))))
                local pd_ok = (!missing(r(min_eigen)) & r(min_eigen) > 1e-8)
            }
            if `pd_ok' {
                display as txt "  Variance matrix:" _col(35) as res "Positive definite"
            }
            else {
                display as txt "  Variance matrix:" _col(35) as err "NOT positive definite"
                display as err "  Note: Numerical Hessian may produce non-PD V; SEs may be approximate"
            }
        }
        else {
            display as txt "  Variance matrix:" _col(35) as err "Contains missing/extreme values"
        }

        * Parameter bounds check
        display _n as txt "Parameter bounds check:"

        if !`is_bivariate' {
            * Loadings are signed under the v_2=1 normalization. Read e(lambda_T)
            * if available (post-Stage-1 package), else fall back to the e(sigma_T)
            * alias. A "reasonable" loading is bounded in magnitude; the sign
            * indicates direction (negative => type 2 has lower hazard).
            *
            * Under factor(common) the same value is exposed as both lambda_T
            * and lambda_Y, so we display a single row labelled "lambda".
            capture local lambda_T = e(lambda_T)
            if missing(`lambda_T') local lambda_T = e(sigma_T)
            capture local lambda_Y = e(lambda_Y)
            if missing(`lambda_Y') local lambda_Y = e(sigma_Y)

            if "`factor_mode'" == "common" {
                if abs(`lambda_T') < 100 {
                    display as txt "  lambda (shared) = " %7.4f `lambda_T' ///
                        _col(35) as res "OK (signed; sign indicates direction)"
                }
                else {
                    display as txt "  lambda (shared) = " %7.4f `lambda_T' ///
                        _col(35) as err "WARNING: extreme magnitude"
                }
            }
            else {
                if abs(`lambda_T') < 100 {
                    display as txt "  lambda_T = " %7.4f `lambda_T' ///
                        _col(35) as res "OK (signed; sign indicates direction)"
                }
                else {
                    display as txt "  lambda_T = " %7.4f `lambda_T' ///
                        _col(35) as err "WARNING: extreme magnitude"
                }

                if abs(`lambda_Y') < 100 {
                    display as txt "  lambda_Y = " %7.4f `lambda_Y' ///
                        _col(35) as res "OK (signed; sign indicates direction)"
                }
                else {
                    display as txt "  lambda_Y = " %7.4f `lambda_Y' ///
                        _col(35) as err "WARNING: extreme magnitude"
                }
            }

            * Mixture probability check
            display _n as txt "Mixture probabilities:"
            local pi_ok = 1
            local spike_slab = 0
            forvalues k = 1/`K' {
                local pi_k = e(pi)[1, `k']
                if `pi_k' < 0.01 {
                    display as txt "  pi_`k' = " %6.4f `pi_k' ///
                        _col(35) as err "WARNING: near boundary (<1%)"
                    local pi_ok = 0
                    * Spike-and-slab corner: a mass point with vanishing
                    * probability paired with a non-trivial risk shift on
                    * either equation. The likelihood typically supports a
                    * flat ridge along (pi_k * lambda * v_k) ~ const, so
                    * BFGS slides into the corner.
                    local v_k = e(v)[1, `k']
                    capture local lam_T = e(lambda_T)
                    if missing(`lam_T') local lam_T = e(sigma_T)
                    capture local lam_Y = e(lambda_Y)
                    if missing(`lam_Y') local lam_Y = e(sigma_Y)
                    local shift_T = abs(`lam_T' * `v_k')
                    local shift_Y = abs(`lam_Y' * `v_k')
                    if `shift_T' > 1.5 | `shift_Y' > 1.5 {
                        local spike_slab = 1
                    }
                }
                else if `pi_k' > 0.99 {
                    display as txt "  pi_`k' = " %6.4f `pi_k' ///
                        _col(35) as err "WARNING: near boundary (>99%)"
                    local pi_ok = 0
                }
                else {
                    display as txt "  pi_`k' = " %6.4f `pi_k' _col(35) as res "OK"
                }
            }

            if `pi_ok' == 0 {
                display as err ///
                    "  Some mixture probabilities near boundary. Consider fewer mass points."
            }
            if `spike_slab' == 1 {
                display _n as err "  SPIKE-AND-SLAB CORNER detected:"
                display as err "  A mass point has pi_k < 1% paired with |lambda*v| > 1.5."
                display as err "  This is the classic NPMLE corner solution -- the data" ///
                    " do not"
                display as err "  identify a separated K-mass distribution at this sample" ///
                    " size."
                display as err "  Consider K-1 mass points, larger N, or a stronger covariate."
            }

            * Mass point separation check
            display _n as txt "Mass point separation:"
            local v_separated = 1
            forvalues k = 1/`K' {
                local v_`k' = e(v)[1, `k']
            }
            forvalues k = 2/`K' {
                local km1 = `k' - 1
                local diff = abs(`v_`k'' - `v_`km1'')
                if `diff' < 0.1 {
                    display as txt "  |v_`k' - v_`km1'| = " %6.4f `diff' ///
                        _col(35) as err "WARNING: poorly separated"
                    local v_separated = 0
                }
            }
            if `v_separated' == 1 {
                display as txt "  Mass points are well separated" _col(35) as res "OK"
            }
            else {
                display as err "  Poorly separated mass points suggest too many types"
            }
        }
        else {
            * Bivariate model checks
            if abs(e(v_T2)) > 0.1 {
                display as txt "  v_T2 = " %6.4f e(v_T2) _col(35) as res "OK"
            }
            else {
                display as txt "  v_T2 = " %6.4f e(v_T2) _col(35) as err "WARNING: near zero"
            }
            if abs(e(v_Y2)) > 0.1 {
                display as txt "  v_Y2 = " %6.4f e(v_Y2) _col(35) as res "OK"
            }
            else {
                display as txt "  v_Y2 = " %6.4f e(v_Y2) _col(35) as err "WARNING: near zero"
            }

            display _n as txt "Joint probabilities:"
            tempname pj
            matrix `pj' = e(pi_joint)
            forvalues i = 1/2 {
                forvalues j = 1/2 {
                    local pij = `pj'[`i', `j']
                    if `pij' < 0.01 | `pij' > 0.99 {
                        display as txt "  pi[`i',`j'] = " %6.4f `pij' ///
                            _col(35) as err "WARNING: near boundary"
                    }
                    else {
                        display as txt "  pi[`i',`j'] = " %6.4f `pij' ///
                            _col(35) as res "OK"
                    }
                }
            }
        }

        return scalar converged = `converged'
        return scalar grad_norm = `grad_norm'
    }

    *==========================================================================
    * POSTERIOR TYPE PROBABILITIES
    *==========================================================================
    if "`posterior'" != "" {
        display _n as txt "{hline 78}"
        display as txt "{bf:MIXTURE DISTRIBUTION SUMMARY}"
        display as txt "{hline 78}"

        if !`is_bivariate' {
            * Loadings are signed under the v_2=1 normalization. Read e(lambda_T)
            * if available, else fall back to the e(sigma_T) alias for back-compat.
            capture local lT_load = e(lambda_T)
            if missing(`lT_load') local lT_load = e(sigma_T)
            capture local lY_load = e(lambda_Y)
            if missing(`lY_load') local lY_load = e(sigma_Y)

            display _n as txt "Prior mixture distribution (population shares):"
            forvalues k = 1/`K' {
                local pi_k = e(pi)[1, `k']
                local v_k = e(v)[1, `k']
                local lT_v = `lT_load' * `v_k'
                local lY_v = `lY_load' * `v_k'
                if "`factor_mode'" == "common" {
                    display as txt "  Type `k': pi = " %5.3f `pi_k' ///
                        ", v = " %6.3f `v_k' ///
                        " (shared shift: " %7.3f `lT_v' ")"
                }
                else {
                    display as txt "  Type `k': pi = " %5.3f `pi_k' ///
                        ", v = " %6.3f `v_k' ///
                        " (treat shift: " %7.3f `lT_v' ", outcome shift: " %7.3f `lY_v' ")"
                }
            }

            * Interpretation. Under v_2=1, the type-k risk shift in the outcome
            * equation is lambda_Y * v_k. Negative shifts denote lower hazard
            * relative to type 1.
            display _n as txt "Interpretation:"
            local v_min = e(v)[1, 1]
            local v_max = e(v)[1, `K']
            forvalues k = 1/`K' {
                if e(v)[1, `k'] < `v_min' local v_min = e(v)[1, `k']
                if e(v)[1, `k'] > `v_max' local v_max = e(v)[1, `k']
            }

            if `v_min' < 0 & `v_max' > 0 {
                display as txt "  Types span v<0 and v>0 mass points"
            }
            else if `v_max' <= 0 {
                display as txt "  All free mass points are at v <= 0"
            }
            else {
                display as txt "  All free mass points are at v >= 0"
            }

            * Loading-ratio interpretation: meaningful only under
            * factor(separate). Under factor(common) the loadings are
            * forced equal so the ratio is identically 1; under
            * factor(separate) opposite signs are admissible and signal
            * negative selection on unobservables.
            if "`factor_mode'" == "common" {
                display as txt "  Loadings constrained equal under factor(common)."
            }
            else if `lT_load' != 0 & sign(`lT_load') == sign(`lY_load') {
                local ratio = `lY_load' / `lT_load'
                display as txt "  |lambda_Y / lambda_T| = " %5.2f abs(`ratio') ///
                    " (heterogeneity loads " ///
                    cond(abs(`ratio')>1, "more on outcome", "more on treatment") ")"
            }
            else if `lT_load' != 0 & sign(`lT_load') != sign(`lY_load') {
                display as txt "  lambda_T and lambda_Y have opposite signs:" ///
                    " heterogeneity drives treatment and outcome in"
                display as txt "  opposite directions. Consider hsmixture_bivariate" ///
                    " (separate v_T and v_Y) as a robustness check, or refit"
                display as txt "  with factor(common) to enforce the manuscript-style" ///
                    " one-factor restriction."
            }
        }
        else {
            display _n as txt "Implied correlation: rho = " %7.4f e(rho)
            if e(rho) > 0 {
                display as txt "  Positive: high-risk treatment types tend to be " ///
                    "high-risk outcome types"
            }
            else if e(rho) < 0 {
                display as txt "  Negative: high-risk treatment types tend to be " ///
                    "low-risk outcome types"
            }
        }

    }

    *==========================================================================
    * MODEL COMPARISON
    *==========================================================================
    if "`compare'" != "" | "`lrtest'" != "" {
        display _n as txt "{hline 78}"
        display as txt "{bf:MODEL COMPARISON}"
        display as txt "{hline 78}"

        * Current model info. BIC denominator is the person count (the IID
        * unit for this mixture model). Pre-v2.3.1 fits that did not store
        * e(N_persons) fall back to e(N).
        local aic = -2 * `ll' + 2 * `k_params'
        local bic = -2 * `ll' + `k_params' * ln(`N_p')

        display as txt "Current model (K=`K'):"
        display as txt "  Log likelihood:" _col(30) %12.4f `ll'
        display as txt "  Parameters:" _col(30) %12.0f `k_params'
        display as txt "  AIC:" _col(30) %12.2f `aic'
        display as txt "  BIC:" _col(30) %12.2f `bic' ///
            "  (N = " %7.0fc `N_p' " persons)"

        * IC-reliability flag. When the variance matrix is rank-deficient
        * the BIC penalty is still well defined (uses the design parameter
        * count; the v2.2.1 fix posts e(rank)=colsof(e(b)) so that
        * estimates stats agrees), but the log-likelihood itself was
        * obtained at a point where the gradient need not vanish. AIC/BIC
        * differences across models are then *not* reliable comparators.
        capture local v_pd_cur = e(v_pd)
        if _rc == 0 & `v_pd_cur' == 0 {
            display as err ///
                "  WARNING: variance matrix is not positive definite at this fit."
            display as err ///
                "  AIC/BIC may not be comparable across models. Re-fit from"
            display as err ///
                "  alternative starts or reduce K before comparing."
        }

        return scalar aic = `aic'
        return scalar bic = `bic'

        * LR test if specified
        if "`lrtest'" != "" {
            capture estimates describe `lrtest'
            if _rc != 0 {
                display as error "Model `lrtest' not found in stored estimates"
            }
            else {
                * Hold the current estimates so we can read e() from `lrtest'
                * and then restore the user's active model. The previous v2.3.0
                * implementation used `estimates restore .`, which is not
                * documented Stata syntax and silently left e() pointing at the
                * alternative model.
                tempname _hsj_hold
                quietly _estimates hold `_hsj_hold', restore
                quietly estimates restore `lrtest'
                local ll_alt = e(ll)
                local k_alt = e(k)
                capture local K_alt = e(K)
                local cmd_alt = e(cmd)
                local factor_alt = "`e(factor)'"
                quietly _estimates unhold `_hsj_hold'

                * Compatibility guard. The alternative must be the same command as
                * the active fit; a cross-command target (e.g. hsmixture_bivariate,
                * which stores no e(K)) does not nest the joint K-mixture and would
                * otherwise read an undefined e(K) into `K_alt'.
                if "`cmd_alt'" != "`e(cmd)'" {
                    display as error "lrtest() target is `cmd_alt' but the active fit is `e(cmd)';"
                    display as error "  the two are not nested and cannot be LR-compared."
                    exit 498
                }

                * Classify the test to label the reference distribution correctly.
                * Same K with a differing factor structure (factor(common) nested
                * in factor(separate) via lambda_T=lambda_Y) is a regular interior
                * restriction: the standard chi-square(df) reference is valid. A
                * K-vs-K' comparison adds a mass point whose probability sits on
                * the [0,1] boundary under the null, where chi-square is invalid.
                local same_k = 0
                if "`K_alt'" != "" & "`K_alt'" != "." {
                    if `K_alt' == `K' local same_k = 1
                }

                local lr_stat = 2 * (`ll' - `ll_alt')
                local df = `k_params' - `k_alt'

                if `df' > 0 {
                    if `is_bivariate' {
                        display _n as txt "LR test vs `lrtest':"
                    }
                    else {
                        display _n as txt "LR test vs `lrtest' (K=`K_alt'):"
                    }
                    display as txt "  LR statistic:" _col(30) %12.4f `lr_stat'
                    display as txt "  Degrees of freedom:" _col(30) %12.0f `df'

                    * Reference-distribution caveat depends on WHAT is tested.
                    * A same-K factor restriction (or any nested restriction on a
                    * fixed mixture structure, as in bivariate-vs-bivariate) is
                    * interior, so the standard chi-square(df) reference is valid.
                    * A K-vs-K' comparison adds a mass point whose probability
                    * lives on [0,1], so under the null 2*(L_alt - L_null) is
                    * *not* chi-square(df) but a mixture of chi-squares.
                    if `same_k' | `is_bivariate' {
                        display as txt "  Interior parameter restriction (same mixture"
                        display as txt "  structure, e.g. factor(common) vs factor(separate)"
                        display as txt "  at fixed K, or nested covariates); the standard"
                        display as txt "  chi-square(df) reference applies."
                    }
                    else {
                        display as err "  CAUTION: standard chi-square reference is invalid"
                        display as err "  for K-selection (mixture-on-boundary problem)."
                        display as err "  Treat the p-value below as a heuristic upper bound."
                        display as err "  Prefer BIC, or a parametric-bootstrap LR distribution."
                    }

                    if `lr_stat' < 0 {
                        * Negative LR means the larger-K fit landed at a
                        * worse local optimum than the smaller-K fit.
                        * The chi-square p-value is meaningless here.
                        display as txt "  p-value:" _col(30) %12s "(not computed)"
                        display as err "  WARNING: LR statistic is negative."
                        display as err "  The larger-K (K=`K') fit found a worse" ///
                            " log-likelihood than the smaller-K (K=`K_alt') fit."
                        display as err "  This indicates the larger model converged" ///
                            " to an inferior local optimum."
                        display as err "  Re-estimate K=`K' from multiple starting" ///
                            " values via from() and select the best log-likelihood" ///
                            " before comparing."
                        return scalar lr_stat = `lr_stat'
                        return scalar lr_df = `df'
                        return scalar lr_pval = .
                    }
                    else {
                        local p_val = chi2tail(`df', `lr_stat')
                        display as txt "  p-value:" _col(30) %12.4f `p_val'

                        if `p_val' < 0.05 {
                            display as txt "  Conclusion: K=`K' significantly better " ///
                                "than K=`K_alt' at 5% level"
                        }
                        else {
                            display as txt "  Conclusion: Cannot reject K=`K_alt' " ///
                                "in favor of K=`K'"
                        }

                        return scalar lr_stat = `lr_stat'
                        return scalar lr_df = `df'
                        return scalar lr_pval = `p_val'
                    }
                }
                else {
                    display as txt ///
                        "Note: Current model has fewer parameters; use reverse comparison"
                }
            }
        }

        * Guidance on K selection
        display _n as txt "Guidance on selecting K:"
        display as txt "  - Restrict comparisons to fits with e(converged)==1"
        display as txt "  - The chi-square LR test is non-regular for K vs K-1"
        display as txt "    (extra mass point sits on the boundary of pi in [0,1])"
        display as txt "  - AIC and BIC from non-converged fits are not interpretable"
        display as txt "  - Discard fits with spike-and-slab corner solutions"
        display as txt "  - Poorly separated mass points suggest K is too large"
    }

    display as txt "{hline 78}"

end
