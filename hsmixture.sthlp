{smcl}
{* *! version 2.3.1  07may2026}{...}
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
{cmd:k(2)}. Any {cmd:k(#)} >= 2 is supported. Whether a given K is identified
depends on the data, not the package. Report results only when
{cmd:e(converged) == 1}. The command warns when BFGS did not converge, the
relative gradient is too large, or the variance matrix is not positive definite.

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
