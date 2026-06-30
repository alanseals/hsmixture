*! Mata functions for hsmixture_bivariate
*! Uses Mata optimize() directly
*!
*! Authors: Jonghoon Park and R. Alan Seals (Auburn University)

version 14
mata:
mata set matastrict off

// ----------------------------------------------------------------
// Cache initialization
// ----------------------------------------------------------------
void _hsb_init_cache()
{
    external real scalar    _hsb_initialized
    external real matrix    _hsb_X_treat, _hsb_X_outcome
    external real colvector _hsb_treat_ev, _hsb_outcome_ev, _hsb_treat_ind, _hsb_riskset
    external real matrix    _hsb_info
    external real scalar    _hsb_N, _hsb_n_persons
    external real scalar    _hsb_k_treat, _hsb_k_outcome

    string scalar   samp_var, id_name, treat_ev_name, outcome_ev_name
    string scalar   treat_name, riskset_name
    string rowvector treat_vars, outcome_vars
    real colvector  samp_vals, idx, pid

    // Read globals
    id_name          = st_global("HSB_id")
    treat_ev_name    = st_global("HSB_treat_event")
    outcome_ev_name  = st_global("HSB_outcome_event")
    treat_name       = st_global("HSB_treat")
    riskset_name     = st_global("HSB_riskset")
    treat_vars       = tokens(st_global("HSB_treat_vars_mata"))
    outcome_vars     = tokens(st_global("HSB_outcome_vars_mata"))

    // Identify estimation sample from Stata's touse variable
    samp_var       = st_local("touse")
    samp_vals = st_data(., samp_var)
    idx = selectindex(samp_vals :!= 0)
    _hsb_N = length(idx)

    // Cache data vectors
    _hsb_treat_ev   = st_data(idx, treat_ev_name)
    _hsb_outcome_ev = st_data(idx, outcome_ev_name)
    _hsb_treat_ind  = st_data(idx, treat_name)

    if (riskset_name != "") {
        _hsb_riskset = st_data(idx, riskset_name)
    }
    else {
        _hsb_riskset = J(_hsb_N, 1, 1)
    }

    // Design matrices: [covariates, constant]
    _hsb_X_treat   = (st_data(idx, treat_vars), J(_hsb_N, 1, 1))
    _hsb_X_outcome = (st_data(idx, outcome_vars), J(_hsb_N, 1, 1))
    _hsb_k_treat   = cols(_hsb_X_treat)
    _hsb_k_outcome = cols(_hsb_X_outcome)

    // Panel structure
    pid = st_data(idx, id_name)
    _hsb_info = panelsetup(pid, 1)
    _hsb_n_persons = rows(_hsb_info)

    _hsb_initialized = 1

    if (st_global("HSB_nolog") == "") {
        printf("{txt}Mata cache (bivariate): %g obs, %g persons\n",
               _hsb_N, _hsb_n_persons)
    }
}

// ----------------------------------------------------------------
// Core log-likelihood: takes parameter vector, returns scalar LL
// ----------------------------------------------------------------
real scalar _hsb_compute_ll(real rowvector b)
{
    external real matrix    _hsb_X_treat, _hsb_X_outcome
    external real colvector _hsb_treat_ev, _hsb_outcome_ev, _hsb_treat_ind, _hsb_riskset
    external real matrix    _hsb_info
    external real scalar    _hsb_N, _hsb_n_persons
    external real scalar    _hsb_k_treat, _hsb_k_outcome

    real scalar     N, n_persons, k_treat, k_outcome, pos
    real colvector  beta_T, beta_Y
    real scalar     delta, v_T2, v_Y2
    real scalar     logit_pi_12, logit_pi_21, logit_pi_22, sum_exp
    real rowvector  pi_k, vT_type, vY_type
    real scalar     i, t, r1, r2
    real colvector  xb_treat, xb_outcome
    real matrix     person_ll_t
    real colvector  eta_T, eta_Y, h_T, h_Y
    real colvector  ll_obs_T, ll_obs_Y, ll_obs
    real matrix     a_t
    real colvector  max_a, ll_person
    real scalar     total_ll

    N = _hsb_N
    n_persons = _hsb_n_persons
    k_treat = _hsb_k_treat
    k_outcome = _hsb_k_outcome

    // Extract parameters
    beta_T = b[1..k_treat]'
    beta_Y = b[(k_treat+1)..(k_treat+k_outcome)]'
    pos    = k_treat + k_outcome
    delta  = b[pos+1]
    v_T2   = b[pos+2]
    v_Y2   = b[pos+3]
    logit_pi_12 = b[pos+4]
    logit_pi_21 = b[pos+5]
    logit_pi_22 = b[pos+6]

    // Softmax probabilities
    sum_exp = 1 + exp(logit_pi_12) + exp(logit_pi_21) + exp(logit_pi_22)
    pi_k = (1, exp(logit_pi_12), exp(logit_pi_21), exp(logit_pi_22)) / sum_exp

    // 4 joint types
    vT_type = (0, 0,    v_T2, v_T2)
    vY_type = (0, v_Y2, 0,    v_Y2)

    // Base linear predictors
    xb_treat   = _hsb_X_treat * beta_T
    xb_outcome = _hsb_X_outcome * beta_Y

    // Type-specific person-level log-likelihoods
    person_ll_t = J(n_persons, 4, 0)

    for (t = 1; t <= 4; t++) {
        eta_T = xb_treat :+ vT_type[t]
        eta_T = rowmax((eta_T, J(N, 1, -20)))
        eta_T = rowmin((eta_T, J(N, 1, 10)))
        h_T = 1 :- exp(-exp(eta_T))
        h_T = rowmax((h_T, J(N, 1, 1e-20)))
        h_T = rowmin((h_T, J(N, 1, 1 - 1e-20)))
        ll_obs_T = _hsb_treat_ev :* ln(h_T) + (1 :- _hsb_treat_ev) :* (-exp(eta_T))
        ll_obs_T = ll_obs_T :* _hsb_riskset

        eta_Y = xb_outcome :+ delta * _hsb_treat_ind :+ vY_type[t]
        eta_Y = rowmax((eta_Y, J(N, 1, -20)))
        eta_Y = rowmin((eta_Y, J(N, 1, 10)))
        h_Y = 1 :- exp(-exp(eta_Y))
        h_Y = rowmax((h_Y, J(N, 1, 1e-20)))
        h_Y = rowmin((h_Y, J(N, 1, 1 - 1e-20)))
        ll_obs_Y = _hsb_outcome_ev :* ln(h_Y) + (1 :- _hsb_outcome_ev) :* (-exp(eta_Y))

        ll_obs = ll_obs_T + ll_obs_Y
        for (i = 1; i <= n_persons; i++) {
            r1 = _hsb_info[i, 1]
            r2 = _hsb_info[i, 2]
            person_ll_t[i, t] = quadsum(ll_obs[|r1 \ r2|])
        }
    }

    // Log-sum-exp mixture
    a_t = J(n_persons, 4, 0)
    for (t = 1; t <= 4; t++) {
        a_t[., t] = max((ln(pi_k[t]), -700)) :+ person_ll_t[., t]
    }
    max_a = rowmax(a_t)
    ll_person = max_a + ln(rowsum(exp(a_t :- max_a)))

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
void _hsb_optim_eval(real scalar todo, real rowvector p,
                     real scalar lnf, real rowvector g, real matrix H)
{
    external real scalar _hsb_initialized
    if (_hsb_initialized != 1) _hsb_init_cache()

    lnf = _hsb_compute_ll(p)
}

// ----------------------------------------------------------------
// Run optimization (called from ado)
// ----------------------------------------------------------------
void _hsb_run_optimize(string scalar b0_name, real scalar max_iter)
{
    external real scalar _hsb_initialized

    real rowvector theta0, theta_hat, g
    real scalar    ll_hat, conv, rc
    real matrix    V
    string scalar  tracelevel

    // Initialize cache if needed
    if (_hsb_initialized != 1) _hsb_init_cache()

    // Get starting values from Stata matrix
    theta0 = st_matrix(b0_name)

    // Respect nolog option
    tracelevel = (st_global("HSB_nolog") != "") ? "none" : "value"

    S = optimize_init()
    optimize_init_evaluator(S, &_hsb_optim_eval())
    optimize_init_evaluatortype(S, "d0")
    optimize_init_params(S, theta0)
    optimize_init_which(S, "max")
    optimize_init_technique(S, "bfgs")
    optimize_init_conv_maxiter(S, max_iter)
    // Convergence criteria. Tightened from the v2.0.0 vtol=1e-6 +
    // ignorenrtol="on" setting that let BFGS quit on relative-LL change
    // alone. nrtol=1e-5 is the floor; tighter values cause BFGS to enter a
    // contraction-reset loop on the locally flat optimum.
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
        st_matrix("__hsb_b", theta_hat)
        st_numscalar("__hsb_ll", ll_hat)
        st_numscalar("__hsb_has_result", 1)
        st_numscalar("__hsb_converged", conv)
        st_numscalar("__hsb_ic", optimize_result_iterations(S))

        // Skip posting V/gradient if they contain missings (singular Hessian
        // at flat spot). The .ado falls through to scaffold-V; ereturn post
        // would otherwise fail with "matrix has missing values".
        V = optimize_result_V_oim(S)
        if (!hasmissing(V)) {
            st_matrix("__hsb_V", V)
        }

        g = optimize_result_gradient(S)
        if (!hasmissing(g)) {
            st_matrix("__hsb_g", g)
        }

        if (st_global("HSB_nolog") == "") {
            printf("{txt}Completed in %g iterations (converged=%g). Log-likelihood = %12.4f\n",
                   optimize_result_iterations(S), conv, ll_hat)
        }
    }
    else {
        st_numscalar("__hsb_has_result", 0)
        st_numscalar("__hsb_converged", 0)
        if (st_global("HSB_nolog") == "") {
            printf("{err}Optimization failed after %g iterations (rc=%g, ll=%g)\n",
                   max_iter, rc, ll_hat)
        }
    }
}

// ----------------------------------------------------------------
// ML evaluator (backward compatibility)
// ----------------------------------------------------------------
void _hsb_eval(string scalar b_name, string scalar lnf_name)
{
    external real scalar _hsb_initialized
    if (_hsb_initialized != 1) _hsb_init_cache()

    real rowvector b
    b = st_matrix(b_name)
    st_numscalar(lnf_name, _hsb_compute_ll(b))
}

// ----------------------------------------------------------------
// Cleanup
// ----------------------------------------------------------------
void _hsb_cleanup()
{
    external real scalar _hsb_initialized
    _hsb_initialized = 0
}

end

display as text "hsmixture_bivariate Mata functions compiled successfully."
