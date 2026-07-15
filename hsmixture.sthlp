{smcl}
{* *! version 2.4.0  14jul2026}{...}
{vieweralsosee "hsmixture_joint" "help hsmixture_joint"}{...}
{vieweralsosee "hsmixture_bivariate" "help hsmixture_bivariate"}{...}
{vieweralsosee "glm" "help glm"}{...}
{vieweralsosee "ml" "help ml"}{...}
{viewerjumpto "Syntax" "hsmixture##syntax"}{...}
{viewerjumpto "Description" "hsmixture##description"}{...}
{viewerjumpto "Methodology" "hsmixture##methodology"}{...}
{viewerjumpto "Options" "hsmixture##options"}{...}
{viewerjumpto "Examples" "hsmixture##examples"}{...}
{viewerjumpto "Stored results" "hsmixture##results"}{...}
{viewerjumpto "Methods and formulas" "hsmixture##methods"}{...}
{viewerjumpto "References" "hsmixture##references"}{...}
{viewerjumpto "Version history" "hsmixture##version"}{...}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col:{cmd:hsmixture} {hline 2}}Discrete-time hazard model with Heckman-Singer mixture heterogeneity{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:hsmixture} {depvar} [{indepvars}] {ifin}{cmd:,} {opth id(varname)} [{it:options}]


{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{p2coldent:* {opth id(varname)}}panel identifier (required){p_end}
{synopt:{opt k(#)}}number of mass points; default is {cmd:k(2)}{p_end}
{synopt:{opth time(varname)}}within-person time variable used by the data-contract checks{p_end}

{syntab:Maximization}
{synopt:{opt from(matname)}}initial values{p_end}
{synopt:{opt iter:ate(#)}}maximum iterations; default is {cmd:iterate(200)}{p_end}
{synopt:{opt nolog}}suppress iteration log{p_end}

{syntab:Reporting}
{synopt:{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}* {opt id()} is required.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:hsmixture} estimates a single-equation discrete-time hazard model with
{bf:Heckman-Singer (1984)} unobserved heterogeneity. This is a building block
for the joint timing-of-events model implemented in {help hsmixture_joint}.

{pstd}
The model assumes individuals belong to one of K latent types, each with a
different baseline hazard level. This discrete mixture approach, introduced
by {bf:Heckman and Singer (1984)} in their {it:Econometrica} paper, provides a
flexible, nonparametric approximation to arbitrary heterogeneity distributions.

{pstd}
The model uses the complementary log-log (cloglog) link function, which is
the discrete-time analog of the continuous-time proportional hazard model
(Prentice and Gloeckler 1978; Jenkins 1995).

{pstd}
{bf:Note:} Factor variables ({cmd:i.varname}) are not supported. Users should
create dummy variables manually, e.g., {cmd:tab period, gen(pd_)}.

{pstd}
{bf:Data contract.} The outcome is absorbing, and the likelihood treats every
estimation row as an at-risk period. Each person's rows must therefore stop at
the outcome event row: at most one event per id, and no rows after the event.
The command validates both and exits with an error on violations. Within-person
order for the check comes from {opt time()} when supplied, otherwise from the
current row order (so keep the data sorted by id and time, as required).


{marker methodology}{...}
{title:Methodology: Why Heckman-Singer?}

{pstd}
{bf:The Problem:} Duration models are highly sensitive to assumptions about
unobserved heterogeneity. Misspecifying the heterogeneity distribution can
severely bias coefficient estimates, including treatment effects.

{pstd}
{bf:The Heckman-Singer Solution:} Rather than assuming a parametric distribution
(e.g., gamma, normal), approximate the unknown distribution with K mass points:

{phang2}- v_1, v_2, ..., v_K are the mass point locations{p_end}
{phang2}- pi_1, pi_2, ..., pi_K are the population shares{p_end}

{pstd}
{bf:Key Results from Heckman and Singer (1984):}

{phang2}1. Discrete mixtures can consistently approximate continuous distributions{p_end}
{phang2}2. This approach avoids the specification bias from parametric assumptions{p_end}

{pstd}
{bf:Connection to Timing-of-Events:} This single-equation model is used within
{cmd:hsmixture_joint} to implement the {bf:Abbring and van den Berg (2003)}
timing-of-events framework, where correlated heterogeneity across treatment
and outcome processes is modeled via shared mass points.


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{opth id(varname)} specifies the variable identifying panels (individuals).
This is required. Data must be in person-period format.

{phang}
{opt k(#)} specifies the number of mass points (latent types). Default is
{cmd:k(2)}. Any {cmd:k(#)} >= 2 is accepted, but parameter recovery is
certified only at K=2; K>=3 is syntactically supported yet not validated
against a known data-generating process and is prone to boundary/spike
solutions. Whether a given K is identified
depends on the data, not the package. Report results only when
{cmd:e(converged) == 1}. The command warns when BFGS did not converge, the
relative gradient is too large, or the variance matrix is not positive definite.

{phang}
{opth time(varname)} supplies the within-person time variable used by the
data-contract checks (numeric, nonmissing, distinct within id). When omitted,
the checks use the current row order within id, which is correct when the data
are sorted by id and time. {opt time()} does not change the likelihood; the
estimator itself is order-invariant within person.

{dlgtab:Maximization}

{phang}
{opt from(matname)} specifies a matrix of initial values. If not specified,
starting values are obtained from a standard cloglog GLM.

{phang}
{opt iterate(#)} specifies the maximum number of iterations. Default is 200.
The Mata-accelerated optimizer uses BFGS with automatic tolerance settings.

{pstd}
{bf:Note on starting values:} {cmd:hsmixture} runs a single optimization start
from the GLM coefficients. Mixture surfaces are multimodal at any K, so a
single start can land at a local mode. The joint and bivariate estimators
({helpb hsmixture_joint} and {helpb hsmixture_bivariate}) use multistart with
6-7 configurations. For single-equation models, refit with at least one
alternative starting vector via {opt from()} and check that the log
likelihoods agree before reporting results.

{dlgtab:Reporting}

{phang}
{opt level(#)} specifies the confidence level, as a percentage, for confidence
intervals. The default is {cmd:level(95)} or as set by {helpb set level}.


{marker examples}{...}
{title:Examples}

{pstd}Setup (run {cmd:hsmixture_example.do} first to generate data){p_end}
{phang2}{cmd:. do hsmixture_example.do}{p_end}
{phang2}{cmd:. use hsmixture_example_data, clear}{p_end}

{pstd}Basic K=2 model{p_end}
{phang2}{cmd:. hsmixture outcome_event x1 pd_2-pd_20, id(id) k(2)}{p_end}

{pstd}K=3 model with more iterations{p_end}
{phang2}{cmd:. hsmixture outcome_event x1 pd_2-pd_20, id(id) k(3) iterate(200)}{p_end}

{pstd}Compare AIC/BIC across specifications (gated on strict convergence){p_end}
{phang2}{cmd:. local model_list}{p_end}
{phang2}{cmd:. forvalues k = 2/3 {c -(}}{p_end}
{phang2}{cmd:.     capture noisily hsmixture outcome_event x1 pd_2-pd_20, id(id) k(`k')}{p_end}
{phang2}{cmd:.     if _rc == 0 & e(converged) == 1 {c -(}}{p_end}
{phang2}{cmd:.         estimates store m_k`k'}{p_end}
{phang2}{cmd:.         local model_list `model_list' m_k`k'}{p_end}
{phang2}{cmd:.     {c )-}}{p_end}
{phang2}{cmd:.     else {c -(}}{p_end}
{phang2}{cmd:.         display as err "K=`k' did not strictly converge; not stored."}{p_end}
{phang2}{cmd:.     {c )-}}{p_end}
{phang2}{cmd:. {c )-}}{p_end}
{phang2}{cmd:. if "`model_list'" != "" {c -(}}{p_end}
{phang2}{cmd:.     estimates stats `model_list'}{p_end}
{phang2}{cmd:. {c )-}}{p_end}
{phang2}{cmd:. else {c -(}}{p_end}
{phang2}{cmd:.     display as err "No models strictly converged; nothing to compare."}{p_end}
{phang2}{cmd:. {c )-}}{p_end}

{pstd}{it:Note on BIC.} Stata's {cmd:estimates stats} computes BIC from the row
count {cmd:e(N)} (person-periods). {cmd:hsmixture} instead reports BIC on the
person count {cmd:e(N_persons)}, the independent unit for this mixture, so the
two BIC columns differ and can even rank K differently (the larger row count
inflates the per-parameter penalty). The base {cmd:hsmixture} command has no
person-count comparison path, so read the {cmd:estimates stats} BIC as
row-count based and compare it against the person-count BIC printed in each
fit's own output.{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:hsmixture} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(ll)}}log likelihood{p_end}
{synopt:{cmd:e(K)}}number of mass points{p_end}
{synopt:{cmd:e(converged)}}1 if {it:strictly} converged (BFGS converged AND |gradient|/(1+|LL|) < 1e-5 AND V positive definite), 0 otherwise{p_end}
{synopt:{cmd:e(rel_grad)}}relative gradient norm |gradient|/(1+|LL|) at the optimum{p_end}
{synopt:{cmd:e(converged_bfgs)}}BFGS optimizer's own convergence flag{p_end}
{synopt:{cmd:e(grad_norm)}}L2 norm of the gradient at the optimum{p_end}
{synopt:{cmd:e(v_pd)}}1 if variance matrix is positive definite, 0 otherwise{p_end}
{synopt:{cmd:e(v_scaffold)}}1 if {cmd:e(V)} is a placeholder posted after Hessian-inversion failure (SEs meaningless), 0 otherwise{p_end}
{synopt:{cmd:e(v_mineig)}}smallest eigenvalue of {cmd:e(V)} (the input to the {cmd:e(v_pd)} verdict){p_end}
{synopt:{cmd:e(clip_hits)}}linear-predictor evaluations at the numerical bounds [-20, 10] at the optimum (0 = exact likelihood everywhere){p_end}
{synopt:{cmd:e(n_aborted)}}1 if the first optimizer run aborted and was retried from a jittered start, 0 otherwise (see {helpb hsmixture_joint}){p_end}
{synopt:{cmd:e(N_persons)}}number of persons (IID unit; denominator for the person-count BIC){p_end}
{synopt:{cmd:e(ic)}}iteration count at the best optimum{p_end}
{synopt:{cmd:e(rank)}}rank posted for {cmd:e(V)} (design parameter count){p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom (design parameter count){p_end}
{synopt:{cmd:e(level)}}confidence level used{p_end}
{synopt:{cmd:e(lambda)}}factor loading (signed){p_end}
{synopt:{cmd:e(sigma)}}alias for {cmd:e(lambda)} (back-compat){p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:hsmixture}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}dependent variable{p_end}
{synopt:{cmd:e(idvar)}}panel variable{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix{p_end}
{synopt:{cmd:e(gradient)}}gradient at the optimum{p_end}
{synopt:{cmd:e(pi)}}mixture probabilities (1 x K){p_end}
{synopt:{cmd:e(v)}}mass points (1 x K){p_end}


{marker methods}{...}
{title:Methods and formulas}

{pstd}
The likelihood for individual i is a mixture over K types:

{pmore}
L_i = SUM_{k=1}^K pi_k * L_i(k)

{pstd}
where the type-specific contribution is:

{pmore}
L_i(k) = PROD_t [h_it(k)]^{d_it} * [1 - h_it(k)]^{1 - d_it}

{pstd}
The hazard function uses the complementary log-log link:

{pmore}
h_it(k) = 1 - exp(-exp(X_it' beta + lambda * v_k))

{pstd}
Parameters are constrained as:

{phang2}- v_1 = 0 (location normalization){p_end}
{phang2}- v_2 = 1 (scale normalization; removes lambda*v rescaling redundancy){p_end}
{phang2}- v_3,...,v_K free real-line parameters (when K >= 3){p_end}
{phang2}- lambda signed real-line factor loading (sign indicates direction of type-2 risk){p_end}
{phang2}- pi_k in (0,1), SUM pi_k = 1 (softmax parameterization){p_end}

{pstd}
The signed-loading parameterization replaces the older log-positive form
(sigma = exp(ln_sigma) > 0 with v_2 free), which had a one-dimensional
rescaling redundancy that interfered with identification. Under v_2=1,
lambda is the type-2 risk shift and the legacy {cmd:e(sigma)} alias points
to the same value.

{pstd}
{bf:Label switching in practice (K=2).} The likelihood is exactly invariant
under swapping the type labels: (pi_1 <-> pi_2, lambda -> -lambda) with the
constant absorbing the shift (cons -> cons + lambda). Both labelings are the
same mixture and the same maximum, and the optimizer can terminate at either
depending on starting values and the platform's floating-point path. The
{it:sign} of lambda, the identity of "type 2", and the specific value of
pi_2 are therefore not comparable across runs or machines; |lambda| and the
unordered shares {c -(}pi_1, pi_2{c )-} are. Report and compare those.


{marker references}{...}
{title:References}

{pstd}
{bf:Unobserved Heterogeneity:}

{phang}
Heckman, J.J. and B. Singer. 1984. A method for minimizing the impact of
distributional assumptions in econometric models for duration data.
{it:Econometrica} 52(2): 271-320.

{phang}
Heckman, J.J. and B. Singer. 1984. Econometric duration analysis.
{it:Journal of Econometrics} 24(1-2): 63-132.

{pstd}
{bf:Timing-of-Events Framework:}

{phang}
Abbring, J.H. and G.J. van den Berg. 2003. The nonparametric identification of
treatment effects in duration models. {it:Econometrica} 71(5): 1491-1517.

{pstd}
{bf:Discrete-Time Hazards:}

{phang}
Jenkins, S.P. 1995. Easy estimation methods for discrete-time duration models.
{it:Oxford Bulletin of Economics and Statistics} 57(1): 129-138.

{phang}
Prentice, R.L. and L.A. Gloeckler. 1978. Regression analysis of grouped survival
data with application to breast cancer data. {it:Biometrics} 34(1): 57-67.

{pstd}
{bf:General Duration Analysis:}

{phang}
Van den Berg, G.J. 2001. Duration models: Specification, identification and
multiple durations. In {it:Handbook of Econometrics}, Vol. 5, ed. J.J. Heckman
and E. Leamer, 3381-3460. Amsterdam: Elsevier.


{marker version}{...}
{title:Version history}

{phang}2.4.0  14jul2026  Data-contract hardening and numerical diagnostics.
    New absorbing-outcome row check errors when estimation rows follow a
    person's outcome event (previously such rows silently entered the
    likelihood as spurious at-risk periods). New time() option supplies
    explicit within-person order for the check; when it is omitted the check
    reads the current row order and the command notes when the data are not
    sorted by the id variable. Mixture probabilities now computed via
    max-shifted softmax (overflow-safe; algebraically identical, so values
    agree to floating-point rounding). An optimizer run that aborts
    (numeric-derivative failure near a flat optimum, a knife-edge event under
    Stata/MP) is now retried once from a jittered start; e(n_aborted) records
    it. The positive-definiteness verdict is now scale-relative rather than an
    absolute 1e-8 floor (see hsmixture_joint). New diagnostics e(clip_hits), e(v_scaffold),
    e(v_mineig). The likelihood formulas and estimates on contract-conforming
    data are unchanged.{p_end}
{phang}2.3.3  02jul2026  Coordinated package version bump. The strict-convergence
    gating of the hazard-ratio confidence interval applies to hsmixture_joint and
    hsmixture_bivariate; this single-equation command estimates no treatment
    effect and posts no CI, so it is unaffected. Mata numerical core
    unchanged.{p_end}
{phang}2.3.2  30jun2026  Coordinated package version bump. No functional change
    to this command; see hsmixture_bivariate for an e(rho) calculation fix.{p_end}
{phang}2.3.1  07may2026  Markout of covariates moved before data-contract
    asserts. e(N_persons) stored; BIC denominator switched to person count
    (the IID unit for this mixture model). Version aligned with package release.{p_end}
{phang}2.2.1  02may2026  Convergence fix: removed ignorenrtol; tightened to
    vtol=1e-10/nrtol=1e-5; gradient gate now relative |grad|/(1+|LL|) < 1e-5.
    Added e(rel_grad).{p_end}
{phang}2.2.0  22mar2026  Switched to Mata optimize(). Dropped factor variable
    support. Nobreak cleanup.{p_end}
{phang}2.1.0  30apr2026  v_2=1 + signed-lambda parameterization (replaced
    log-positive sigma). e(sigma) preserved as alias.{p_end}
{phang}2.0.0  21mar2026  Generalized for SSC distribution. Added level() option.{p_end}
{phang}1.1.0  09feb2026  Switched to mleval for parameter extraction.{p_end}
{phang}1.0.0  17jan2026  Initial release.{p_end}


{title:Citation}

{pstd}
If you use this package, please cite:

{phang}
Park, J. and R.A. Seals. 2026. hsmixture: Stata module for discrete-time
hazard models with Heckman-Singer mixture heterogeneity.


{title:Authors}

{pstd}
Jonghoon Park and R. Alan Seals{break}
Department of Economics{break}
Auburn University{break}
{browse "mailto:jzp0200@auburn.edu":jzp0200@auburn.edu} /
{browse "mailto:ras0029@auburn.edu":ras0029@auburn.edu}


{title:Also see}

{p 4 14 2}
{help hsmixture_joint:hsmixture_joint} (joint timing-of-events model),
{help hsmixture_bivariate:hsmixture_bivariate} (bivariate heterogeneity model),
{help glm:glm},
{help ml:ml}
