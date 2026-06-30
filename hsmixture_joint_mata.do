*! Mata functions for hsmixture_joint
*! Uses Mata optimize() directly
*!
*! Authors: Jonghoon Park and R. Alan Seals (Auburn University)

version 14
mata:
mata set matastrict off

// ----------------------------------------------------------------
// Cache initialization
// ----------------------------------------------------------------
void _hsj_init_cache()
{
    external real scalar    _hsj_initialized
    external real matrix    _hsj_X_treat, _hsj_X_outcome
    external real colvector _hsj_treat_ev, _hsj_outcome_ev, _hsj_treat_ind, _hsj_riskset
    external real matrix    _hsj_info
    external real scalar    _hsj_N, _hsj_n_persons, _hsj_K
    external real scalar    _hsj_k_treat, _hsj_k_outcome

    string scalar   samp_var, id_name, treat_ev_name, outcome_ev_name
    string scalar   treat_name, riskset_name
    string rowvector treat_vars, outcome_vars
    real colvector  samp_vals, idx, pid

    // Read globals set by hsmixture_joint.ado
    _hsj_K           = strtoreal(st_global("HSJ_K"))
    id_name          = st_global("HSJ_id")
    treat_ev_name    = st_global("HSJ_treat_event")
    outcome_ev_name  = st_global("HSJ_outcome_event")
    treat_name       = st_global("HSJ_treat")
    riskset_name     = st_global("HSJ_riskset")
    treat_vars       = tokens(st_global("HSJ_treat_vars_mata"))
    outcome_vars     = tokens(st_global("HSJ_outcome_vars_mata"))

    // Identify estimation sample from Stata's touse variable
    samp_var       = st_local("touse")
    samp_vals = st_data(., samp_var)
    idx = selectindex(samp_vals :!= 0)
    _hsj_N = length(idx)

    // Cache event and treatment vectors
    _hsj_treat_ev   = st_data(idx, treat_ev_name)
    _hsj_outcome_ev = st_data(idx, outcome_ev_name)
    _hsj_treat_ind  = st_data(idx, treat_name)

    if (riskset_name != "") {
        _hsj_riskset = st_data(idx, riskset_name)
    }
    else {
        _hsj_riskset = J(_hsj_N, 1, 1)
    }

    // Build design matrices: [covariates, constant]
    _hsj_X_treat   = (st_data(idx, treat_vars), J(_hsj_N, 1, 1))
    _hsj_X_outcome = (st_data(idx, outcome_vars), J(_hsj_N, 1, 1))
    _hsj_k_treat   = cols(_hsj_X_treat)
    _hsj_k_outcome = cols(_hsj_X_outcome)

    // Panel structure for person-level aggregation
    pid = st_data(idx, id_name)
    _hsj_info = panelsetup(pid, 1)
    _hsj_n_persons = rows(_hsj_info)

    _hsj_initialized = 1

    if (st_global("HSJ_nolog") == "") {
        printf("{txt}Mata cache (joint): %g obs, %g persons, K=%g\n",
               _hsj_N, _hsj_n_persons, _hsj_K)
    }
}

// ----------------------------------------------------------------
// Core log-likelihood: takes parameter vector, returns scalar LL
// ----------------------------------------------------------------
real scalar _hsj_compute_ll(real rowvector b)
{
    external real matrix    _hsj_X_treat, _hsj_X_outcome
    external real colvector _hsj_treat_ev, _hsj_outcome_ev, _hsj_treat_ind, _hsj_riskset
    external real matrix    _hsj_info
    external real scalar    _hsj_N, _hsj_n_persons, _hsj_K
    external real scalar    _hsj_k_treat, _hsj_k_outcome

    real scalar     K, N, n_persons, k_treat, k_outcome, pos, n_loadings
    real scalar     factor_common
    real colvector  beta_T, beta_Y
    real scalar     delta, lambda_T, lambda_Y
    real rowvector  v, eta_raw, pi_k
    real scalar     sum_exp_eta, j, k, i
    real colvector  xb_treat, xb_outcome
    real matrix     person_ll_k
    real colvector  eta_T, eta_Y, h_T, h_Y
    real colvector  ll_obs_T, ll_obs_Y, ll_obs
    real scalar     r1, r2
    real matrix     a_k
    real colvector  max_a, ll_person
    real scalar     total_ll

    K = _hsj_K
    N = _hsj_N
    n_persons = _hsj_n_persons
    k_treat = _hsj_k_treat
    k_outcome = _hsj_k_outcome

    // Resolve factor structure. Empty global => "separate" (back-compat
    // with any caller that pre-dates the factor() option).
    factor_common = (st_global("HSJ_factor") == "common")
    n_loadings = factor_common ? 1 : 2

    // Extract parameters from b vector. v_2 = 1 normalization removes the
    // lambda * v rescaling redundancy (alpha-rescaling of v cancelled by
    // 1/alpha-rescaling of the loading(s)). With v_1=0 and v_2=1 fixed, the
    // remaining mass points v_3..v_K and the loading(s) are separately
    // identified. Loadings are real-line (signed).
    //
    // Layout under factor(separate):
    //   [beta_T (k_treat) | beta_Y (k_outcome) | delta |
    //    lambda_T | lambda_Y | v_3..v_K | eta_2..eta_K]
    // Layout under factor(common):
    //   [beta_T (k_treat) | beta_Y (k_outcome) | delta |
    //    lambda           | v_3..v_K | eta_2..eta_K]
    // For K=2 the v_3..v_K block is empty.
    beta_T   = b[1..k_treat]'
    beta_Y   = b[(k_treat+1)..(k_treat+k_outcome)]'
    pos      = k_treat + k_outcome
    delta    = b[pos+1]
    if (factor_common) {
        lambda_T = b[pos+2]
        lambda_Y = lambda_T
    }
    else {
        lambda_T = b[pos+2]
        lambda_Y = b[pos+3]
    }

    // Mass points: v_1=0 and v_2=1 fixed; v_3..v_K live at pos + n_loadings + j - 1.
    // (separate: n_loadings=2, v_3 at pos+4 = pos+1+j when j=3, matches v2.x;
    //  common:   n_loadings=1, v_3 at pos+3.)
    v = J(1, K, 0)
    if (K >= 2) v[2] = 1
    for (j = 3; j <= K; j++) {
        v[j] = b[pos + n_loadings + j - 1]
    }

    // Mixture weights via softmax: eta_1 = 0 normalized; eta_2..eta_K live at
    // pos + n_loadings + K + j - 2.
    // (separate: eta_2 at pos+K+2, matches v2.x; common: eta_2 at pos+K+1.)
    eta_raw = J(1, K, 0)
    for (j = 2; j <= K; j++) {
        eta_raw[j] = b[pos + n_loadings + K + j - 2]
    }
    sum_exp_eta = sum(exp(eta_raw))
    pi_k = exp(eta_raw) / sum_exp_eta

    // Base linear predictors
    xb_treat   = _hsj_X_treat * beta_T
    xb_outcome = _hsj_X_outcome * beta_Y

    // Type-specific person-level log-likelihoods
    person_ll_k = J(n_persons, K, 0)

    for (k = 1; k <= K; k++) {
        eta_T = xb_treat :+ (lambda_T * v[k])
        eta_T = rowmax((eta_T, J(N, 1, -20)))
        eta_T = rowmin((eta_T, J(N, 1, 10)))
        h_T = 1 :- exp(-exp(eta_T))
        h_T = rowmax((h_T, J(N, 1, 1e-20)))
        h_T = rowmin((h_T, J(N, 1, 1 - 1e-20)))
        ll_obs_T = _hsj_treat_ev :* ln(h_T) + (1 :- _hsj_treat_ev) :* (-exp(eta_T))
        ll_obs_T = ll_obs_T :* _hsj_riskset

        eta_Y = xb_outcome :+ delta * _hsj_treat_ind :+ (lambda_Y * v[k])
        eta_Y = rowmax((eta_Y, J(N, 1, -20)))
        eta_Y = rowmin((eta_Y, J(N, 1, 10)))
        h_Y = 1 :- exp(-exp(eta_Y))
        h_Y = rowmax((h_Y, J(N, 1, 1e-20)))
        h_Y = rowmin((h_Y, J(N, 1, 1 - 1e-20)))
        ll_obs_Y = _hsj_outcome_ev :* ln(h_Y) + (1 :- _hsj_outcome_ev) :* (-exp(eta_Y))

        ll_obs = ll_obs_T + ll_obs_Y
        for (i = 1; i <= n_persons; i++) {
            r1 = _hsj_info[i, 1]
            r2 = _hsj_info[i, 2]
            person_ll_k[i, k] = quadsum(ll_obs[|r1 \ r2|])
        }
    }

    // Log-sum-exp mixture
    a_k = J(n_persons, K, 0)
    for (k = 1; k <= K; k++) {
        a_k[., k] = max((ln(pi_k[k]), -700)) :+ person_ll_k[., k]
    }
    max_a = rowmax(a_k)
    ll_person = max_a + ln(rowsum(exp(a_k :- max_a)))

    // Replace missing with large penalty
    for (i = 1; i <= n_persons; i++) {
        if (missing(ll_person[i])) ll_person[i] = -1e10
    }

    total_ll = quadsum(ll_person)
    if (missing(total_ll)) total_ll = -1e20
    return(total_ll)
}

// ----------------------------------------------------------------
// optimize() evaluator
// ----------------------------------------------------------------
void _hsj_optim_eval(real scalar todo, real rowvector p,
                     real scalar lnf, real rowvector g, real matrix H)
{
    external real scalar _hsj_initialized
    if (_hsj_initialized != 1) _hsj_init_cache()

    lnf = _hsj_compute_ll(p)
}

// ----------------------------------------------------------------
// Run optimization (called from ado)
// ----------------------------------------------------------------
void _hsj_run_optimize(string scalar b0_name, real scalar max_iter)
{
    external real scalar _hsj_initialized

    real rowvector theta0, theta_hat, g
    real scalar    ll_hat, conv, rc
    real matrix    V
    string scalar  tracelevel

    // Initialize cache if needed
    if (_hsj_initialized != 1) _hsj_init_cache()

    // Get starting values from Stata matrix
    theta0 = st_matrix(b0_name)

    // Respect nolog option
    tracelevel = (st_global("HSJ_nolog") != "") ? "none" : "value"

    S = optimize_init()
    optimize_init_evaluator(S, &_hsj_optim_eval())
    optimize_init_evaluatortype(S, "d0")
    optimize_init_params(S, theta0)
    optimize_init_which(S, "max")
    optimize_init_technique(S, "bfgs")
    optimize_init_conv_maxiter(S, max_iter)
    // Convergence criteria. The v2.0.0 vtol=1e-6 + ignorenrtol="on" setting
    // let BFGS quit on relative-LL change alone (cert reported converged=1
    // with |grad|=6.77). With these tighter settings the cert achieves
    // parameter recovery, V positive definite, and |grad|/(1+|f|) ~ 3e-7
    // at termination. Tightening further (nrtol=1e-9) caused BFGS to enter
    // a contraction-reset loop near the flat optimum, so 1e-5 is the floor.
    optimize_init_conv_ptol(S, 1e-10)
    optimize_init_conv_vtol(S, 1e-10)
    optimize_init_conv_nrtol(S, 1e-5)
    optimize_init_singularHmethod(S, "hybrid")
    optimize_init_tracelevel(S, tracelevel)

    // Run optimization
    rc = _optimize(S)

    // Extract results
    theta_hat = optimize_result_params(S)
    ll_hat = optimize_result_value(S)
    conv = optimize_result_converged(S)

    if (!missing(ll_hat) & ll_hat > -1e15) {
        // Post results: always post when ll is finite (best-of-starts selection)
        st_matrix("__hsj_b", theta_hat)
        st_numscalar("__hsj_ll", ll_hat)
        st_numscalar("__hsj_has_result", 1)
        st_numscalar("__hsj_converged", conv)
        st_numscalar("__hsj_ic", optimize_result_iterations(S))

        // Hessian-based V. When the optimizer hits a flat spot (singular
        // Hessian) optimize_result_V_oim() returns a matrix with missing
        // values. We do NOT post that matrix — the .ado then falls through
        // to its scaffold-V branch so `ereturn post` doesn't fail with
        // "matrix has missing values". A clean V is posted only when the
        // numerical Hessian was invertible.
        V = optimize_result_V_oim(S)
        if (!hasmissing(V)) {
            st_matrix("__hsj_V", V)
        }

        // Gradient at the optimum (read by hsmixture_joint_postestimation
        // to report a finite gradient norm). Skip posting if the gradient
        // itself contains missings.
        g = optimize_result_gradient(S)
        if (!hasmissing(g)) {
            st_matrix("__hsj_g", g)
        }

        if (st_global("HSJ_nolog") == "") {
            printf("{txt}Completed in %g iterations (converged=%g). Log-likelihood = %12.4f\n",
                   optimize_result_iterations(S), conv, ll_hat)
        }
    }
    else {
        st_numscalar("__hsj_has_result", 0)
        st_numscalar("__hsj_converged", 0)
        if (st_global("HSJ_nolog") == "") {
            printf("{err}Optimization failed after %g iterations (rc=%g, ll=%g)\n",
                   max_iter, rc, ll_hat)
        }
    }
}

// ----------------------------------------------------------------
// ML evaluator (backward compatibility)
// ----------------------------------------------------------------
void _hsj_eval(string scalar b_name, string scalar lnf_name)
{
    external real scalar _hsj_initialized
    if (_hsj_initialized != 1) _hsj_init_cache()

    real rowvector b
    b = st_matrix(b_name)
    st_numscalar(lnf_name, _hsj_compute_ll(b))
}

// ----------------------------------------------------------------
// Cleanup
// ----------------------------------------------------------------
void _hsj_cleanup()
{
    external real scalar _hsj_initialized
    _hsj_initialized = 0
}

end

display as text "hsmixture_joint Mata functions compiled successfully."
