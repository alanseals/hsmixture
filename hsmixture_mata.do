*! Mata functions for hsmixture (single-equation)
*! Uses Mata optimize() directly
*!
*! Authors: Jonghoon Park and R. Alan Seals (Auburn University)

version 14
mata:
mata set matastrict off

// ----------------------------------------------------------------
// Cache initialization
// ----------------------------------------------------------------
void _hs_init_cache()
{
    external real scalar    _hs_initialized
    external real matrix    _hs_X
    external real colvector _hs_depvar
    external real matrix    _hs_info
    external real scalar    _hs_N, _hs_n_persons, _hs_K
    external real scalar    _hs_k_x

    string scalar   samp_var, id_name, depvar_name
    string rowvector xvars
    real colvector  samp_vals, idx, pid

    _hs_K         = strtoreal(st_global("HS_K"))
    id_name       = st_global("HS_id")
    depvar_name   = st_global("HS_depvar")
    xvars         = tokens(st_global("HS_xvars_mata"))

    // Identify estimation sample
    samp_var  = st_local("touse")
    samp_vals = st_data(., samp_var)
    idx = selectindex(samp_vals :!= 0)
    _hs_N = length(idx)

    // Cache data
    _hs_depvar = st_data(idx, depvar_name)
    _hs_X = (st_data(idx, xvars), J(_hs_N, 1, 1))
    _hs_k_x = cols(_hs_X)

    // Panel structure
    pid = st_data(idx, id_name)
    _hs_info = panelsetup(pid, 1)
    _hs_n_persons = rows(_hs_info)

    _hs_initialized = 1

    if (st_global("HS_nolog") == "") {
        printf("{txt}Mata cache (single): %g obs, %g persons, K=%g\n",
               _hs_N, _hs_n_persons, _hs_K)
    }
}

// ----------------------------------------------------------------
// Core log-likelihood
// Parameter layout: [beta (k_x) | lambda | v_3..v_K | eta_2..eta_K]
// (v_1 = 0 and v_2 = 1 fixed by scale normalization; lambda is signed real)
// ----------------------------------------------------------------
real scalar _hs_compute_ll(real rowvector b)
{
    external real matrix    _hs_X
    external real colvector _hs_depvar
    external real matrix    _hs_info
    external real scalar    _hs_N, _hs_n_persons, _hs_K
    external real scalar    _hs_k_x

    real scalar     K, N, n_persons, k_x, pos
    real colvector  beta
    real scalar     lambda, j, k, i
    real rowvector  v, eta_raw, pi_k
    real scalar     sum_exp_eta
    real colvector  xb
    real matrix     person_ll_k
    real colvector  eta_k, h_k, ll_obs
    real scalar     r1, r2
    real matrix     a_k
    real colvector  max_a, ll_person
    real scalar     total_ll

    K = _hs_K
    N = _hs_N
    n_persons = _hs_n_persons
    k_x = _hs_k_x

    // Extract parameters. v_2=1 normalization removes the lambda*v rescaling
    // redundancy. v_1=0, v_2=1 fixed; v_3..v_K and lambda are separately
    // identified. Lambda is real-line (signed) so type 2 can be either
    // higher- or lower-risk than type 1.
    //
    // Layout: [beta (k_x) | lambda | v_3..v_K | eta_2..eta_K]
    beta  = b[1..k_x]'
    pos   = k_x
    lambda = b[pos+1]

    // Mass points: v_1=0 and v_2=1 fixed; v_3..v_K read from b at pos+k-1.
    v = J(1, K, 0)
    if (K >= 2) v[2] = 1
    for (j = 3; j <= K; j++) {
        v[j] = b[pos + j - 1]
    }

    // Mixture weights via softmax: eta_1 = 0 normalized; eta_2..eta_K read
    // from b at pos+K+k-2. Max-shifted before exponentiating so a large
    // logit cannot overflow exp(); the shift cancels in the ratio.
    eta_raw = J(1, K, 0)
    for (j = 2; j <= K; j++) {
        eta_raw[j] = b[pos + K + j - 2]
    }
    eta_raw = eta_raw :- max(eta_raw)
    sum_exp_eta = sum(exp(eta_raw))
    pi_k = exp(eta_raw) / sum_exp_eta

    // Base linear predictor
    xb = _hs_X * beta

    // Type-specific person-level log-likelihoods
    person_ll_k = J(n_persons, K, 0)

    for (k = 1; k <= K; k++) {
        eta_k = xb :+ (lambda * v[k])
        eta_k = rowmax((eta_k, J(N, 1, -20)))
        eta_k = rowmin((eta_k, J(N, 1, 10)))
        h_k = 1 :- exp(-exp(eta_k))
        h_k = rowmax((h_k, J(N, 1, 1e-20)))
        h_k = rowmin((h_k, J(N, 1, 1 - 1e-20)))

        // ll_obs = d * ln(h) + (1-d) * (-exp(eta))  [cloglog stable form]
        ll_obs = _hs_depvar :* ln(h_k) + (1 :- _hs_depvar) :* (-exp(eta_k))

        for (i = 1; i <= n_persons; i++) {
            r1 = _hs_info[i, 1]
            r2 = _hs_info[i, 2]
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

    for (i = 1; i <= n_persons; i++) {
        if (missing(ll_person[i])) ll_person[i] = -1e10
    }

    total_ll = quadsum(ll_person)
    if (missing(total_ll)) total_ll = -1e20
    return(total_ll)
}

// ----------------------------------------------------------------
// Clip-hit diagnostic
// Counts linear-predictor evaluations outside the numerical bounds
// [-20, 10] at parameter vector b (one count per row x type). The
// likelihood truncates eta to that range for overflow safety; a nonzero
// count at the optimum means part of the fitted surface is the clipped
// (approximate) likelihood. Reported to the .ado as e(clip_hits).
// ----------------------------------------------------------------
real scalar _hs_count_clips(real rowvector b)
{
    external real matrix _hs_X
    external real scalar _hs_K, _hs_k_x

    real colvector beta, xb, eta_k
    real scalar    lambda, j, k, pos, nclip
    real rowvector v

    beta   = b[1.._hs_k_x]'
    pos    = _hs_k_x
    lambda = b[pos+1]

    v = J(1, _hs_K, 0)
    if (_hs_K >= 2) v[2] = 1
    for (j = 3; j <= _hs_K; j++) {
        v[j] = b[pos + j - 1]
    }

    xb = _hs_X * beta
    nclip = 0
    for (k = 1; k <= _hs_K; k++) {
        eta_k = xb :+ (lambda * v[k])
        nclip = nclip + sum((eta_k :< -20) :| (eta_k :> 10))
    }
    return(nclip)
}

// ----------------------------------------------------------------
// optimize() evaluator
// ----------------------------------------------------------------
void _hs_optim_eval(real scalar todo, real rowvector p,
                    real scalar lnf, real rowvector g, real matrix H)
{
    external real scalar _hs_initialized
    if (_hs_initialized != 1) _hs_init_cache()

    lnf = _hs_compute_ll(p)
}

// ----------------------------------------------------------------
// Run optimization (called from ado)
// ----------------------------------------------------------------
void _hs_run_optimize(string scalar b0_name, real scalar max_iter)
{
    external real scalar _hs_initialized

    real rowvector theta0, theta_hat, g
    real scalar    ll_hat, conv, rc
    real matrix    V
    string scalar  tracelevel

    if (_hs_initialized != 1) _hs_init_cache()

    theta0 = st_matrix(b0_name)

    tracelevel = (st_global("HS_nolog") != "") ? "none" : "value"

    S = optimize_init()
    optimize_init_evaluator(S, &_hs_optim_eval())
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

    rc = _optimize(S)

    theta_hat = optimize_result_params(S)
    ll_hat = optimize_result_value(S)
    conv = optimize_result_converged(S)

    if (!missing(ll_hat) & ll_hat > -1e15) {
        st_matrix("__hs_b", theta_hat)
        st_numscalar("__hs_ll", ll_hat)
        st_numscalar("__hs_has_result", 1)
        st_numscalar("__hs_converged", conv)
        st_numscalar("__hs_ic", optimize_result_iterations(S))
        st_numscalar("__hs_clip", _hs_count_clips(theta_hat))

        // Skip posting V/gradient if they contain missings (singular Hessian
        // at flat spot). The .ado falls through to scaffold-V; ereturn post
        // would otherwise fail.
        V = optimize_result_V_oim(S)
        if (!hasmissing(V)) {
            st_matrix("__hs_V", V)
        }

        g = optimize_result_gradient(S)
        if (!hasmissing(g)) {
            st_matrix("__hs_g", g)
        }

        if (st_global("HS_nolog") == "") {
            printf("{txt}Completed in %g iterations (converged=%g). Log-likelihood = %12.4f\n",
                   optimize_result_iterations(S), conv, ll_hat)
        }
    }
    else {
        st_numscalar("__hs_has_result", 0)
        st_numscalar("__hs_converged", 0)
        if (st_global("HS_nolog") == "") {
            printf("{err}Optimization failed (rc=%g, ll=%g)\n", rc, ll_hat)
        }
    }
}

// ----------------------------------------------------------------
// Cleanup
// ----------------------------------------------------------------
void _hs_cleanup()
{
    external real scalar _hs_initialized
    _hs_initialized = 0
}

end

display as text "hsmixture Mata functions compiled successfully."
