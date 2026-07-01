{smcl}
{* *! version 2.3.2  30jun2026}{...}
{vieweralsosee "hsmixture" "help hsmixture"}{...}
{vieweralsosee "hsmixture_bivariate" "help hsmixture_bivariate"}{...}
{vieweralsosee "hsmixture_joint_postestimation" "help hsmixture_joint_postestimation"}{...}
{vieweralsosee "glm" "help glm"}{...}
{vieweralsosee "ml" "help ml"}{...}
{viewerjumpto "Syntax" "hsmixture_joint##syntax"}{...}
{viewerjumpto "Description" "hsmixture_joint##description"}{...}
{viewerjumpto "Timing-of-Events Framework" "hsmixture_joint##toe"}{...}
{viewerjumpto "Identification" "hsmixture_joint##identification"}{...}
{viewerjumpto "Options" "hsmixture_joint##options"}{...}
{viewerjumpto "Examples" "hsmixture_joint##examples"}{...}
{viewerjumpto "Diagnostics" "hsmixture_joint##diagnostics"}{...}
{viewerjumpto "Stored results" "hsmixture_joint##results"}{...}
{viewerjumpto "Methods and formulas" "hsmixture_joint##methods"}{...}
{viewerjumpto "References" "hsmixture_joint##references"}{...}
{viewerjumpto "Version history" "hsmixture_joint##version"}{...}
{title:Title}

{p2colset 5 23 25 2}{...}
{p2col:{cmd:hsmixture_joint} {hline 2}}Discrete-time timing-of-events model with Heckman-Singer heterogeneity{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:hsmixture_joint} {cmd:(}{it:treat_depvar} {cmd:=} {it:treat_indepvars}{cmd:)}
    {cmd:(}{it:outcome_depvar} {cmd:=} {it:outcome_indepvars}{cmd:,} {opt treat(varname)}{cmd:)}
    {ifin}{cmd:,} {opth id(varname)} [{it:options}]


{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{p2coldent:* {opth id(varname)}}panel identifier (required){p_end}
{synopt:{opt k(#)}}number of mass points; default is {cmd:k(2)}{p_end}
{synopt:{opth riskset(varname)}}treatment risk set indicator{p_end}
{synopt:{opt fact:or(name)}}heterogeneity factor structure: {cmd:common} or {cmd:separate} (default){p_end}

{syntab:Maximization}
{synopt:{opt from(matname)}}initial values{p_end}
{synopt:{opt iter:ate(#)}}maximum iterations; default is {cmd:iterate(100)}{p_end}
{synopt:{opt nst:arts(#)}}number of multistart configurations; default 7 (separate) or 6 (common){p_end}
{synopt:{opt nolog}}suppress iteration log{p_end}

{syntab:Reporting}
{synopt:{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}* {opt id()} is required.{p_end}
{p 4 6 2}* {opt treat()} is required in the second equation.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:hsmixture_joint} implements a {bf:discrete-time} version of the {bf:timing-of-events (ToE)}
framework developed by {bf:Abbring and van den Berg (2003)} in their seminal {it:Econometrica}
paper "The Nonparametric Identification of Treatment Effects in Duration Models."

{pstd}
The original Abbring-van den Berg framework uses continuous-time mixed proportional hazard
(MPH) models. This command adapts that framework to {bf:discrete-time hazards} using the
complementary log-log link, which is the discrete-time analog of the continuous-time
proportional hazard model (see Prentice and Gloeckler 1978; Jenkins 1995).

{pstd}
Unobserved heterogeneity is modeled using the {bf:Heckman-Singer (1984)} discrete mixture
approach, which approximates an arbitrary heterogeneity distribution with K mass points.
This nonparametric specification avoids distributional assumptions that can bias treatment
effect estimates in duration models.

{pstd}
The key output is {bf:delta}, the causal effect of treatment on the outcome hazard,
expressed as a log hazard ratio. The hazard ratio exp(delta) gives the multiplicative
effect of treatment on the instantaneous probability of the outcome.


{marker toe}{...}
{title:The Timing-of-Events Framework}

{pstd}
{bf:Core Idea:} The ToE approach identifies causal effects by exploiting {it:variation in
when} treatment occurs, not just {it:whether} it occurs. Two individuals with identical
observable and unobservable characteristics who experience treatment at different times
provide identifying variation.

{pstd}
{bf:The Selection Problem:} In observational studies, treated individuals typically differ
from untreated individuals in ways that also affect outcomes. Simple comparisons confound
treatment effects with selection.

{pstd}
{bf:ToE Solution:} By jointly modeling the treatment and outcome processes with correlated
unobserved heterogeneity, the ToE approach separates:

{phang2}1. {it:Selection effects}: Some individuals are more likely to both experience
treatment AND have the outcome (captured by correlated heterogeneity v_k){p_end}

{phang2}2. {it:Causal effects}: The direct impact of treatment on the outcome hazard
(captured by delta){p_end}

{pstd}
{bf:Relationship to Abbring-van den Berg:} Their Theorem 1 establishes nonparametric
identification of treatment effects in continuous-time MPH models under:

{phang2}(a) No-anticipation: Treatment does not affect outcomes before it occurs{p_end}
{phang2}(b) Proportional hazards with time-varying treatment indicator{p_end}
{phang2}(c) Sufficient variation in treatment timing{p_end}

{pstd}
This command implements a discrete-time analog of their framework. The complementary
log-log link ensures that discrete-time hazard ratios have the same interpretation as
continuous-time hazard ratios (Jenkins 1995).


{marker identification}{...}
{title:Identification}

{pstd}
{bf:Identification in the ToE Framework} relies on three key elements:

{pstd}
{bf:1. Variation in Treatment Timing}

{pmore}
Identification requires individuals who eventually receive treatment to do so at
different times. If all treated individuals were treated in period 1, we could not
distinguish treatment effects from selection. The variation in {it:when} treatment
occurs, conditional on observables and unobservables, provides identifying information.

{pstd}
{bf:2. The No-Anticipation Assumption}

{pmore}
Treatment must not affect outcomes {it:before} it occurs. If individuals change behavior
in anticipation of treatment, this contaminates the pre-treatment comparison group.
Formally: E[Y(0)_t | T_treatment = s] does not depend on s for t < s.

{pstd}
{bf:3. The Proportional Hazard Assumption}

{pmore}
Treatment affects the hazard multiplicatively: h_treated(t) = exp(delta) * h_untreated(t).
This means treatment shifts the entire hazard function up or down by a constant factor.

{pstd}
{bf:Role of the Heckman-Singer Mixture:}

{pmore}
The discrete mixture (K mass points) approximates the joint distribution of
unobserved factors affecting both treatment and outcome. Heckman and Singer
(1984) showed that under regularity conditions a discrete mixture can
consistently approximate a continuous heterogeneity distribution. Whether a
given K is identified on a particular dataset is a finite-sample question and
must be checked via {cmd:e(converged) == 1}.


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

{phang}
{opt treat(varname)} in the second equation specifies the treatment indicator.
This should be a time-varying dummy equal to 1 in all periods after treatment onset.
The treatment effect delta measures the shift in the outcome hazard when this
indicator switches from 0 to 1.

{phang}
{opth riskset(varname)} specifies an indicator variable for observations in the
treatment risk set. When specified, the treatment equation likelihood contribution
is computed only for observations where {it:riskset} == 1. The outcome equation
uses all observations. This is useful when treatment eligibility is conditional
(e.g., only certain individuals are at risk for the treatment in certain periods).

{phang}
{opt factor(name)} selects the heterogeneity factor structure. Two values:

{pmore}
{cmd:factor(separate)} (default). Two free signed loadings, lambda_T and
lambda_Y. Per-type shifts (lambda_T*v_k, lambda_Y*v_k) lie on a one-dimensional
locus through the origin, but the locus direction is data-determined. This
admits opposite signs (negative correlation between treatment-prone and
outcome-prone latent types).

{pmore}
{cmd:factor(common)} forces a single shared loading, lambda. Per-type shifts
(lambda*v_k, lambda*v_k) lie on the 45-degree line through the origin. This
is the classical Heckman-Singer one-factor MPH model in which the same
exp(eta_i) raises both hazards. Correlation between latent treatment and
outcome propensity is mechanically positive. This is the model Abbring and
van den Berg (2003) and most timing-of-events applications describe in their
identification arguments.

{pmore}
The two specifications nest: factor(common) is factor(separate) with the
constraint lambda_T = lambda_Y. Under the null that the constraint holds,
the two HRs should agree to numerical tolerance and a likelihood-ratio test
is well-defined (one degree of freedom). Under the alternative,
factor(separate) is the more flexible model and factor(common) is misspecified.

{dlgtab:Maximization}

{phang}
{opt from(matname)} specifies a matrix of initial values. If not specified,
starting values are obtained from separate cloglog GLMs for each equation.

{phang}
{opt iterate(#)} specifies the maximum number of iterations. Default is 100.
Joint models often require 200+ iterations. The Mata-accelerated optimizer
uses BFGS with automatic tolerance settings.

{dlgtab:Reporting}

{phang}
{opt level(#)} specifies the confidence level, as a percentage, for confidence
intervals. The default is {cmd:level(95)} or as set by {helpb set level}.


{marker examples}{...}
{title:Examples}

{pstd}{bf:Setup} (run {cmd:hsmixture_example.do} first to generate data){p_end}
{phang2}{cmd:. use hsmixture_example_data, clear}{p_end}

{pstd}{bf:Basic joint ToE with K=2}{p_end}
{phang2}{cmd:. hsmixture_joint (treat_event = pd_* x1) ///}{p_end}
{phang2}{cmd:.     (outcome_event = pd_* x1, treat(treated)) ///}{p_end}
{phang2}{cmd:.     , id(id) k(2)}{p_end}

{pstd}{bf:K=3 model}{p_end}
{phang2}{cmd:. hsmixture_joint (treat_event = pd_* x1) ///}{p_end}
{phang2}{cmd:.     (outcome_event = pd_* x1, treat(treated)) ///}{p_end}
{phang2}{cmd:.     , id(id) k(3) iterate(200)}{p_end}

{pstd}{bf:With treatment risk set restriction}{p_end}
{phang2}{cmd:. hsmixture_joint (treat_event = pd_* x1) ///}{p_end}
{phang2}{cmd:.     (outcome_event = pd_* x1, treat(treated)) ///}{p_end}
{phang2}{cmd:.     , id(id) k(2) riskset(at_risk)}{p_end}

{pstd}{bf:Model selection across K (gated on strict convergence)}{p_end}
{phang2}{cmd:. local model_list}{p_end}
{phang2}{cmd:. forvalues k = 2/4 {c -(}}{p_end}
{phang2}{cmd:.     capture noisily hsmixture_joint (...), id(id) k(`k') iterate(200)}{p_end}
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

{pstd}{bf:Diagnostics after estimation}{p_end}
{phang2}{cmd:. hsmixture_joint_postestimation, convergence}{p_end}
{phang2}{cmd:. hsmixture_joint_postestimation, all}{p_end}

{pstd}{bf:LR test: K=3 vs K=2 (both fits must strictly converge)}{p_end}
{phang2}{cmd:. hsmixture_joint (...), id(id) k(2)}{p_end}
{phang2}{cmd:. assert e(converged) == 1}{p_end}
{phang2}{cmd:. estimates store m_k2}{p_end}
{phang2}{cmd:. hsmixture_joint (...), id(id) k(3)}{p_end}
{phang2}{cmd:. assert e(converged) == 1}{p_end}
{phang2}{cmd:. hsmixture_joint_postestimation, lrtest(m_k2)}{p_end}


{marker diagnostics}{...}
{title:Diagnostics and Postestimation}

{pstd}
After estimation, use {cmd:hsmixture_joint_postestimation} for diagnostic checks:

{phang2}{cmd:. hsmixture_joint_postestimation, convergence}{p_end}

{pstd}
This displays:

{phang3}- Convergence status and iteration count{p_end}
{phang3}- Relative gradient norm |grad|/(1+|LL|) (should be < 1e-5){p_end}
{phang3}- Variance matrix check (should be positive definite){p_end}
{phang3}- Parameter bounds (|lambda| < 100, 0 < pi < 1){p_end}
{phang3}- Mass point separation (poorly separated points suggest K too large){p_end}

{pstd}
{bf:Key Diagnostic Checks:}

{phang}
{it:1. Convergence}: If the model did not converge, results are unreliable. Try
increasing {opt iterate()}, providing better starting values with {opt from()},
or reducing K.

{phang}
{it:2. Boundary parameters}: If lambda approaches 0 (no heterogeneity) or pi
approaches 0/1 (a type with vanishing share), the model may be poorly identified.
Consider fewer mass points.

{phang}
{it:3. Mass point separation}: If v_k values are nearly equal, the types are not
well-distinguished. Reduce K.

{phang}
{it:4. Sensitivity to K}: Compare treatment effects across K = 2, 3, 4. Large changes
suggest sensitivity to heterogeneity specification. Report the range.


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:hsmixture_joint} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(ll)}}log likelihood{p_end}
{synopt:{cmd:e(K)}}number of mass points{p_end}
{synopt:{cmd:e(converged)}}1 if {it:strictly} converged (BFGS converged AND |gradient|/(1+|LL|) < 1e-5 AND V positive definite), 0 otherwise{p_end}
{synopt:{cmd:e(converged_bfgs)}}BFGS optimizer's own convergence flag{p_end}
{synopt:{cmd:e(grad_norm)}}L2 norm of the gradient at the optimum{p_end}
{synopt:{cmd:e(rel_grad)}}relative gradient norm |gradient|/(1+|LL|) at the optimum{p_end}
{synopt:{cmd:e(v_pd)}}1 if variance matrix is positive definite, 0 otherwise{p_end}
{synopt:{cmd:e(ic)}}iteration count for the best starting configuration{p_end}
{synopt:{cmd:e(n_starts)}}number of starting configurations tried{p_end}
{synopt:{cmd:e(best_start)}}index of the starting configuration with best LL{p_end}
{synopt:{cmd:e(delta)}}treatment effect (log hazard ratio){p_end}
{synopt:{cmd:e(hr)}}hazard ratio exp(delta){p_end}
{synopt:{cmd:e(hr_ci_lo)}}lower CI bound for hazard ratio{p_end}
{synopt:{cmd:e(hr_ci_hi)}}upper CI bound for hazard ratio{p_end}
{synopt:{cmd:e(lambda)}}shared factor loading (only under {cmd:factor(common)}){p_end}
{synopt:{cmd:e(lambda_T)}}treatment-equation factor loading (signed). Under {cmd:factor(common)} this is an alias for {cmd:e(lambda)}.{p_end}
{synopt:{cmd:e(lambda_Y)}}outcome-equation factor loading (signed). Under {cmd:factor(common)} this is an alias for {cmd:e(lambda)}.{p_end}
{synopt:{cmd:e(sigma_T)}}alias for {cmd:e(lambda_T)} (v2.0.0 back-compat){p_end}
{synopt:{cmd:e(sigma_Y)}}alias for {cmd:e(lambda_Y)} (v2.0.0 back-compat){p_end}
{synopt:{cmd:e(sigma_P)}}alias for {cmd:e(lambda_T)} (v2.0.0 back-compat){p_end}
{synopt:{cmd:e(sigma_D)}}alias for {cmd:e(lambda_Y)} (v2.0.0 back-compat){p_end}
{synopt:{cmd:e(level)}}confidence level used{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:hsmixture_joint}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(treat_depvar)}}treatment dependent variable{p_end}
{synopt:{cmd:e(outcome_depvar)}}outcome dependent variable{p_end}
{synopt:{cmd:e(treat_var)}}treatment indicator variable{p_end}
{synopt:{cmd:e(idvar)}}panel variable{p_end}
{synopt:{cmd:e(riskset_var)}}risk set variable (if specified){p_end}
{synopt:{cmd:e(factor)}}factor structure ({cmd:common} or {cmd:separate}){p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix{p_end}
{synopt:{cmd:e(pi)}}mixture probabilities (1 x K){p_end}
{synopt:{cmd:e(v)}}mass points (1 x K){p_end}


{marker methods}{...}
{title:Methods and formulas}

{pstd}
{bf:Model Structure}

{pstd}
The population consists of K unobserved types. Individual i belongs to type k with
probability pi_k. Conditional on type, two discrete-time hazard processes are observed.
The factor structure governs how the type-k mass point v_k loads into the two hazards.

{pstd}
{bf:factor(separate)} (default). Each equation has its own signed loading.

{pstd}
{it:Treatment hazard}:

{pmore}
h^T_it(k) = 1 - exp(-exp(X^T_it' alpha + lambda_T * v_k))

{pstd}
{it:Outcome hazard}:

{pmore}
h^Y_it(k) = 1 - exp(-exp(X^Y_it' beta + delta * D_it + lambda_Y * v_k))

{pstd}
The pair (lambda_T, lambda_Y) defines the direction of the heterogeneity locus
in (treatment-shift, outcome-shift) space. Type-k shifts are (lambda_T*v_k,
lambda_Y*v_k). Negative correlation (opposite signs) between treatment-prone
and outcome-prone types is admissible.

{pstd}
{bf:factor(common)} forces a single shared loading. The treatment hazard becomes:

{pmore}
h^T_it(k) = 1 - exp(-exp(X^T_it' alpha + lambda * v_k))

{pstd}
and the outcome hazard:

{pmore}
h^Y_it(k) = 1 - exp(-exp(X^Y_it' beta + delta * D_it + lambda * v_k))

{pstd}
The same scalar lambda enters both equations, so per-type shifts (lambda*v_k,
lambda*v_k) lie on the 45-degree line through the origin. Correlation between
latent treatment and outcome propensity is mechanically positive and equal in
magnitude across equations. This is the parameterization usually written down
in the Heckman-Singer (1984) and Abbring-van den Berg (2003) frameworks.

{pstd}
where D_it is the treatment indicator (equal to 1 in all periods after treatment onset).

{pstd}
{bf:Likelihood}

{pstd}
The likelihood for individual i, integrating over unobserved types:

{pmore}
L_i = SUM_{k=1}^K pi_k * L^T_i(k) * L^Y_i(k)

{pstd}
where the type-specific contributions are:

{pmore}
L^T_i(k) = PROD_t [h^T_it(k)]^{y^T_it} * [1 - h^T_it(k)]^{1 - y^T_it}

{pmore}
L^Y_i(k) = PROD_t [h^Y_it(k)]^{y^Y_it} * [1 - h^Y_it(k)]^{1 - y^Y_it}

{pstd}
{bf:Parameterization}

{phang2}- v_1 = 0 (location normalization){p_end}
{phang2}- v_2 = 1 (scale normalization; removes the lambda*v rescaling redundancy){p_end}
{phang2}- v_3,...,v_K free real-line parameters (when K >= 3){p_end}
{phang2}- factor(separate): lambda_T, lambda_Y signed real-line factor loadings (two parameters){p_end}
{phang2}- factor(common): lambda signed real-line factor loading (one parameter){p_end}
{phang2}- pi_k in (0,1), SUM pi_k = 1 (softmax parameterization){p_end}

{pstd}
The signed-loading parameterization replaces the older log-positive parameterization
(where sigma_T = exp(ln_sigma_T) > 0 and v_2 was free). The two parameterizations
imply the same family of likelihood surfaces, but the new one is identified up to
label switching (whereas the old one had a one-dimensional rescaling redundancy:
multiplying all v_k by alpha and dividing both sigmas by alpha gave the same
likelihood). Researchers reading older code that interprets sigma_T as the
type-2 risk shift in the treatment equation should know that under v_2=1,
lambda_T and the type-2 risk shift coincide.

{pstd}
{bf:Connection to Continuous-Time MPH}

{pstd}
The complementary log-log link h(t) = 1 - exp(-exp(eta)) corresponds to a
continuous-time proportional hazard model with piecewise-constant baseline hazard.
This makes discrete-time cloglog estimates directly comparable to continuous-time
hazard ratios (Jenkins 1995).

{pstd}
{bf:Known Limitations}

{phang2}- K > 4 may cause convergence problems with moderate samples{p_end}
{phang2}- Large panels (>50k persons) may require substantial memory for the Mata cache{p_end}
{phang2}- Factor variables ({cmd:i.varname}) are not supported; create dummies manually with {cmd:tab} or indicator variables{p_end}


{marker references}{...}
{title:References}

{pstd}
{bf:Timing-of-Events Identification:}

{phang}
Abbring, J.H. and G.J. van den Berg. 2003. The nonparametric identification of
treatment effects in duration models. {it:Econometrica} 71(5): 1491-1517.

{phang}
Abbring, J.H. and G.J. van den Berg. 2005. Social experiments and instrumental
variables with duration outcomes. {it:Tinbergen Institute Discussion Paper} 05-047/3.

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
{bf:Discrete-Time Hazard Models:}

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

{phang}2.3.2  30jun2026  Coordinated package version bump. No functional change
    to this command; see hsmixture_bivariate for an e(rho) calculation fix.{p_end}
{phang}2.3.1  07may2026  Markout of covariates moved before data-contract
    asserts. e(N_persons) stored; BIC denominator switched to person count
    (the IID unit for this mixture model).{p_end}
{phang}2.3.0  04may2026  Added factor() option (common vs separate loadings).
    Added e(factor) and e(lambda). Multistart now adapts to factor mode.{p_end}
{phang}2.2.1  02may2026  Convergence fix: removed ignorenrtol; tightened to
    vtol=1e-10/nrtol=1e-5; gradient gate now relative |grad|/(1+|LL|) < 1e-5.
    Added e(rel_grad). Stage A 5/5 PASS; opposite-signs cert 7/7 PASS.{p_end}
{phang}2.1.0  30apr2026  v_2=1 + signed-lambda parameterization (replaced
    log-positive sigma). Added prisk()->riskset() rename (prisk kept as alias).
    Sigma aliases preserved for v2.0.0 callers. Equation labels now use
    actual depvar names. markout fix for missing riskset values.{p_end}
{phang}2.2.0  22mar2026  Switched to Mata optimize(). Dropped factor variable support. Nobreak cleanup.{p_end}
{phang}2.0.0  21mar2026  Generalized for SSC. Added riskset(), level(). Renamed
    equations to treat/outcome. Mata-accelerated. Fixed e(converged) bug.{p_end}
{phang}1.2.0  09feb2026  Added risk set option. Improved starting values.{p_end}
{phang}1.1.0  17jan2026  Switched to mleval. Added sortpreserve.{p_end}
{phang}1.0.0  15jan2026  Initial release.{p_end}


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

{psee}
Postestimation: {help hsmixture_joint_postestimation:hsmixture_joint_postestimation}

{psee}
Related: {help hsmixture:hsmixture} (single-equation version),
{help hsmixture_bivariate:hsmixture_bivariate} (bivariate heterogeneity),
{help glm:glm},
{help stcox:stcox},
{help streg:streg},
{help ml:ml}
