{smcl}
{* *! version 2.4.0  14jul2026}{...}
{vieweralsosee "hsmixture" "help hsmixture"}{...}
{vieweralsosee "hsmixture_joint" "help hsmixture_joint"}{...}
{vieweralsosee "hsmixture_joint_postestimation" "help hsmixture_joint_postestimation"}{...}
{viewerjumpto "Syntax" "hsmixture_bivariate##syntax"}{...}
{viewerjumpto "Description" "hsmixture_bivariate##description"}{...}
{viewerjumpto "Options" "hsmixture_bivariate##options"}{...}
{viewerjumpto "Examples" "hsmixture_bivariate##examples"}{...}
{viewerjumpto "Stored results" "hsmixture_bivariate##results"}{...}
{viewerjumpto "Methods and formulas" "hsmixture_bivariate##methods"}{...}
{viewerjumpto "References" "hsmixture_bivariate##references"}{...}
{viewerjumpto "Version history" "hsmixture_bivariate##version"}{...}
{title:Title}

{p2colset 5 27 29 2}{...}
{p2col:{cmd:hsmixture_bivariate} {hline 2}}Bivariate heterogeneity timing-of-events model{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:hsmixture_bivariate} {cmd:(}{it:treat_depvar} {cmd:=} {it:treat_indepvars}{cmd:)}
    {cmd:(}{it:outcome_depvar} {cmd:=} {it:outcome_indepvars}{cmd:,} {opt treat(varname)}{cmd:)}
    {ifin}{cmd:,} {opth id(varname)} [{it:options}]


{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{p2coldent:* {opth id(varname)}}panel identifier (required){p_end}
{synopt:{opth riskset(varname)}}treatment risk set indicator; required whenever treated persons remain observed after their treatment event (the usual ToE data shape){p_end}
{synopt:{opth prisk(varname)}}deprecated alias for {opt riskset()} (v2.0.0 back-compat){p_end}
{synopt:{opth time(varname)}}within-person time variable used by the data-contract checks{p_end}
{synopt:{opt nst:arts(#)}}number of starting configurations; default is {cmd:nstarts(6)}{p_end}

{syntab:Maximization}
{synopt:{opt from(matname)}}initial values (used for first starting config){p_end}
{synopt:{opt iter:ate(#)}}maximum iterations; default is {cmd:iterate(200)}{p_end}
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
{cmd:hsmixture_bivariate} estimates a joint timing-of-events model with
{bf:two-dimensional} unobserved heterogeneity. The treatment-side and
outcome-side mass points are separate (v_T and v_Y) and the joint distribution
over them is parameterized by a free 2x2 probability matrix.

{pstd}
{bf:Relationship to {cmd:hsmixture_joint}:}

{phang2}{cmd:hsmixture_joint, factor(common)}: one shared loading lambda.
Per-type shifts (lambda*v_k, lambda*v_k) lie on the 45-degree line through the
origin. This is the manuscript-style one-factor MPH model.{p_end}

{phang2}{cmd:hsmixture_joint, factor(separate)}: two free signed loadings
lambda_T, lambda_Y. Per-type shifts (lambda_T*v_k, lambda_Y*v_k) lie on a
one-dimensional locus through the origin, but the locus direction is
data-determined. Opposite signs (negative selection) are admissible.{p_end}

{phang2}{cmd:hsmixture_bivariate} (this command): per-type shifts (v_T*j, v_Y*k)
land on a 2x2 grid of corners (0,0), (0,v_Y2), (v_T2,0), (v_T2,v_Y2) with a
free joint probability matrix. As the off-diagonal probabilities go to zero the
model approaches {cmd:factor(separate)} (two free diagonal shifts); it matches
{cmd:factor(common)} only in the further knife-edge case v_T2 = v_Y2. Because
the softmax keeps every cell strictly positive, this is a limiting case rather
than a finite-parameter model.{p_end}

{pstd}
{bf:When to use this:}

{phang2}1. As a robustness check after {cmd:hsmixture_joint} returns
opposite-signs lambdas. The bivariate is the natural sequel because it does
not constrain the four type shifts to a 1-D locus.{p_end}

{phang2}2. When you want to estimate the implied correlation between
treatment-side and outcome-side heterogeneity directly, rather than infer it
from a sign on lambda_T versus lambda_Y.{p_end}

{phang2}3. When the data plausibly support a population structure with
substantial mass on all four corners of the (v_T, v_Y) grid, not just the
diagonal that the joint estimator's locus implies.{p_end}

{pstd}
The model uses a 2x2 grid of joint types: K_T = 2 treatment types crossed with
K_Y = 2 outcome types yields 4 joint types. The joint probabilities pi_{jk} are
freely estimated (via softmax), allowing arbitrary correlation between treatment
and outcome heterogeneity.


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{opth id(varname)} specifies the variable identifying panels (individuals).
This is required. Data must be in person-period format.

{phang}
{opt treat(varname)} in the second equation specifies the treatment indicator.
This should be a time-varying dummy equal to 1 in all periods after treatment onset.
Same-period and lagged switch-on conventions are both accepted; the command
validates that the indicator is absorbing (never reverts to 0) and never equals
1 before the person's treatment event (see {helpb hsmixture_joint} for the full
convention discussion, which applies here unchanged).

{phang}
{opth riskset(varname)} specifies an indicator variable for observations in the
treatment risk set. When specified, the treatment equation likelihood contribution
is computed only for observations where {it:riskset} == 1. The outcome equation
uses all observations. Effectively required: the command errors when
{opt riskset()} is omitted and any treated person has rows after the treatment
event (the normal shape of timing-of-events data), and validates that the risk
set is 1 on the event row and 0 on every row after it. See
{helpb hsmixture_joint} for the rationale.

{phang}
{opth time(varname)} supplies the within-person time variable used by the
data-contract checks (numeric, nonmissing, distinct within id). When omitted,
the checks use the current row order within id, which is correct when the data
are sorted by id and time, and the command prints a note when the data are not
sorted by {opt id()}. {opt time()} does not change the likelihood; the
estimator itself is order-invariant within person.

{phang}
{bf:Data contract.} The outcome is absorbing: at most one outcome event per id
and no estimation rows after the event row. The command validates both and
errors on violations (see {helpb hsmixture_joint} for details).

{phang}
{opt nstarts(#)} specifies the number of starting value configurations to try.
Default is 6, which is also the maximum (the command uses 6 pre-defined grids
of mass point initializations). The optimizer runs from each configuration and
selects the result with the highest log-likelihood.

{dlgtab:Maximization}

{phang}
{opt from(matname)} specifies a matrix of initial values. When provided,
it is used for the first starting configuration. Remaining configurations
use the default grid. If not specified, all starting values come from
separate cloglog GLMs with varying mass point initializations.

{phang}
{opt iterate(#)} specifies the maximum number of iterations. Default is 200.
The Mata-accelerated optimizer uses BFGS with automatic tolerance settings.

{phang}
{opt difficult}, {opt trace}, {opt gradient}, {opt hessian},
{opt technique(algorithm)}, {opt tolerance(#)}, {opt ltolerance(#)}, and
{opt nrtolerance(#)} are accepted for back-compatibility with v2.0.0
({cmd:ml model d0}-era) callers but have {bf:no effect}: estimation uses the
Mata {cmd:optimize()} BFGS routine with fixed tolerances. Supplying any of them
prints a note.

{dlgtab:Reporting}

{phang}
{opt level(#)} specifies the confidence level, as a percentage, for confidence
intervals. The default is {cmd:level(95)} or as set by {helpb set level}.


{marker examples}{...}
{title:Examples}

{pstd}{bf:Setup} (run {cmd:hsmixture_example.do} first to generate data){p_end}
{phang2}{cmd:. use hsmixture_example_data, clear}{p_end}

{pstd}{bf:Basic bivariate model}{p_end}
{phang2}{cmd:. hsmixture_bivariate (treat_event = pd_* x1) ///}{p_end}
{phang2}{cmd:.     (outcome_event = pd_* x1, treat(treated)) ///}{p_end}
{phang2}{cmd:.     , id(id)}{p_end}

{pstd}{bf:With fewer starting configurations (faster)}{p_end}
{phang2}{cmd:. hsmixture_bivariate (treat_event = pd_* x1) ///}{p_end}
{phang2}{cmd:.     (outcome_event = pd_* x1, treat(treated)) ///}{p_end}
{phang2}{cmd:.     , id(id) nstarts(3)}{p_end}

{pstd}{bf:Compare one-factor vs bivariate (gated on strict convergence)}{p_end}
{phang2}{cmd:. local model_list}{p_end}
{phang2}{cmd:. capture noisily hsmixture_joint (treat_event = pd_* x1) ///}{p_end}
{phang2}{cmd:.     (outcome_event = pd_* x1, treat(treated)) ///}{p_end}
{phang2}{cmd:.     , id(id) k(2)}{p_end}
{phang2}{cmd:. if _rc == 0 & e(converged) == 1 {c -(}}{p_end}
{phang2}{cmd:.     estimates store m_1factor}{p_end}
{phang2}{cmd:.     local model_list `model_list' m_1factor}{p_end}
{phang2}{cmd:. {c )-}}{p_end}
{phang2}{cmd:. capture noisily hsmixture_bivariate (treat_event = pd_* x1) ///}{p_end}
{phang2}{cmd:.     (outcome_event = pd_* x1, treat(treated)) ///}{p_end}
{phang2}{cmd:.     , id(id)}{p_end}
{phang2}{cmd:. if _rc == 0 & e(converged) == 1 {c -(}}{p_end}
{phang2}{cmd:.     estimates store m_bivariate}{p_end}
{phang2}{cmd:.     local model_list `model_list' m_bivariate}{p_end}
{phang2}{cmd:. {c )-}}{p_end}
{phang2}{cmd:. if "`model_list'" != "" {c -(}}{p_end}
{phang2}{cmd:.     estimates stats `model_list'}{p_end}
{phang2}{cmd:. {c )-}}{p_end}
{phang2}{cmd:. else {c -(}}{p_end}
{phang2}{cmd:.     display as err "No models strictly converged; nothing to compare."}{p_end}
{phang2}{cmd:. {c )-}}{p_end}

{pstd}{bf:Diagnostics}{p_end}
{phang2}{cmd:. hsmixture_joint_postestimation, all}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:hsmixture_bivariate} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(ll)}}log likelihood{p_end}
{synopt:{cmd:e(converged)}}1 if {it:strictly} converged (BFGS converged AND |gradient|/(1+|LL|) < 1e-5 AND V positive definite), 0 otherwise{p_end}
{synopt:{cmd:e(converged_bfgs)}}BFGS optimizer's own convergence flag{p_end}
{synopt:{cmd:e(grad_norm)}}L2 norm of the gradient at the optimum{p_end}
{synopt:{cmd:e(rel_grad)}}relative gradient norm |gradient|/(1+|LL|) at the optimum{p_end}
{synopt:{cmd:e(v_pd)}}1 if variance matrix is positive definite, 0 otherwise{p_end}
{synopt:{cmd:e(v_scaffold)}}1 if {cmd:e(V)} is a placeholder posted after Hessian-inversion failure (SEs meaningless), 0 otherwise{p_end}
{synopt:{cmd:e(v_mineig)}}smallest eigenvalue of {cmd:e(V)} (the input to the {cmd:e(v_pd)} verdict){p_end}
{synopt:{cmd:e(clip_hits)}}linear-predictor evaluations at the numerical bounds [-20, 10] at the optimum (0 = exact likelihood everywhere){p_end}
{synopt:{cmd:e(N_persons)}}number of persons (IID unit; denominator for the person-count BIC){p_end}
{synopt:{cmd:e(k)}}number of estimated parameters (= colsof(e(b))){p_end}
{synopt:{cmd:e(rank)}}rank posted for {cmd:e(V)} (design parameter count){p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom (design parameter count){p_end}
{synopt:{cmd:e(ic)}}iteration count for the best starting configuration{p_end}
{synopt:{cmd:e(n_finite)}}starting configurations that returned a finite log-likelihood{p_end}
{synopt:{cmd:e(n_bfgs_conv)}}starting configurations whose BFGS run converged on its own criterion{p_end}
{synopt:{cmd:e(n_aborted)}}starting configurations whose first optimizer run aborted and was retried from a jittered start (see {helpb hsmixture_joint}){p_end}
{synopt:{cmd:e(delta)}}treatment effect (log hazard ratio){p_end}
{synopt:{cmd:e(hr)}}hazard ratio exp(delta){p_end}
{synopt:{cmd:e(hr_ci_lo)}}lower CI bound for hazard ratio (missing if not strictly converged){p_end}
{synopt:{cmd:e(hr_ci_hi)}}upper CI bound for hazard ratio (missing if not strictly converged){p_end}
{synopt:{cmd:e(se_delta)}}standard error of delta{p_end}
{synopt:{cmd:e(v_T2)}}second treatment mass point{p_end}
{synopt:{cmd:e(v_Y2)}}second outcome mass point{p_end}
{synopt:{cmd:e(rho)}}implied correlation between v_T and v_Y{p_end}
{synopt:{cmd:e(n_params)}}number of estimated parameters{p_end}
{synopt:{cmd:e(n_starts)}}number of starting configurations tried{p_end}
{synopt:{cmd:e(best_start)}}index of best starting configuration{p_end}
{synopt:{cmd:e(level)}}confidence level used{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:hsmixture_bivariate}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(treat_depvar)}}treatment dependent variable{p_end}
{synopt:{cmd:e(outcome_depvar)}}outcome dependent variable{p_end}
{synopt:{cmd:e(treat_var)}}treatment indicator variable{p_end}
{synopt:{cmd:e(idvar)}}panel variable{p_end}
{synopt:{cmd:e(riskset_var)}}risk set variable (if specified){p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix{p_end}
{synopt:{cmd:e(gradient)}}gradient at the optimum{p_end}
{synopt:{cmd:e(pi_joint)}}joint probability matrix (2 x 2){p_end}


{marker methods}{...}
{title:Methods and formulas}

{pstd}
{bf:Model Structure}

{pstd}
The bivariate model uses separate mass points for the treatment and outcome
equations. The 2x2 joint grid creates 4 types:

{pmore}
Type (1,1): v_T = 0, v_Y = 0 with probability pi_{11}{break}
Type (1,2): v_T = 0, v_Y = v_{Y2} with probability pi_{12}{break}
Type (2,1): v_T = v_{T2}, v_Y = 0 with probability pi_{21}{break}
Type (2,2): v_T = v_{T2}, v_Y = v_{Y2} with probability pi_{22}

{pstd}
The hazard functions are:

{pmore}
h^T_it(j,k) = 1 - exp(-exp(X^T_it' alpha + v_{Tj}))

{pmore}
h^Y_it(j,k) = 1 - exp(-exp(X^Y_it' beta + delta * D_it + v_{Yk}))

{pstd}
The four joint probabilities are parameterized via softmax (with pi_{11} as
the reference category) to ensure they sum to 1 and remain positive.

{pstd}
{bf:Label switching.} Each dimension's type labels can swap independently
(v_T2 -> -v_T2 with the treatment constant absorbing the shift and the
pi rows relabeled; likewise v_Y2 and the pi columns), leaving the likelihood
unchanged. The sign of the v_T2 or v_Y2 {it:parameter} and the position of
any single cell in {cmd:e(pi_joint)} are therefore not comparable across
runs or machines. Quantities that describe the latent joint distribution
itself are unaffected by relabeling and are comparable: delta and the hazard
ratio, |v_T2| and |v_Y2| (the spacing of each dimension's mass points), and
{cmd:e(rho)} {it:including its sign} -- the relabeling is a location shift
of the latent variable, not a reflection, so the implied correlation is
identical under either labeling. The certification script scores only
label-invariant quantities.

{pstd}
{bf:Implied Correlation}

{pstd}
The implied correlation between treatment and outcome heterogeneity is:

{pmore}
rho = Cov(v_T, v_Y) / sqrt(Var(v_T) * Var(v_Y))

{pstd}
where the moments are computed from the joint probability matrix. A positive
rho indicates that individuals with higher treatment risk also tend to have
higher outcome risk (positive selection).


{marker references}{...}
{title:References}

{phang}
Abbring, J.H. and G.J. van den Berg. 2003. The nonparametric identification of
treatment effects in duration models. {it:Econometrica} 71(5): 1491-1517.

{phang}
Heckman, J.J. and B. Singer. 1984. A method for minimizing the impact of
distributional assumptions in econometric models for duration data.
{it:Econometrica} 52(2): 271-320.

{phang}
Jenkins, S.P. 1995. Easy estimation methods for discrete-time duration models.
{it:Oxford Bulletin of Economics and Statistics} 57(1): 129-138.

{phang}
Van den Berg, G.J. 2001. Duration models: Specification, identification and
multiple durations. In {it:Handbook of Econometrics}, Vol. 5, ed. J.J. Heckman
and E. Leamer, 3381-3460. Amsterdam: Elsevier.


{marker version}{...}
{title:Version history}

{phang}2.4.0  14jul2026  Data-contract hardening and numerical diagnostics,
    matching hsmixture_joint: absorbing-outcome row check, riskset()
    effectively required with post-event validation, treat() indicator
    consistency checks, new time() option (with a note when the data are not
    sorted by the id variable and time() is omitted), max-shifted softmax for
    the cell probabilities (likelihood and posted e(pi_joint); algebraically
    identical, so values agree to floating-point rounding), and new diagnostics
    e(clip_hits), e(v_scaffold), e(v_mineig), e(n_finite), e(n_bfgs_conv),
    e(n_aborted). A starting configuration whose optimizer run aborts
    (numeric-derivative failure near a flat optimum, a knife-edge event under
    Stata/MP) is now retried once from a jittered start instead of silently
    killing the remaining configurations (see hsmixture_joint for details).
    The positive-definiteness verdict is now scale-relative rather than an
    absolute 1e-8 floor (see hsmixture_joint).
    The likelihood formulas and estimates on contract-conforming data are
    unchanged.{p_end}
{phang}2.3.3  02jul2026  Report the hazard-ratio confidence interval and the
    positive-definite verdict only for strictly converged fits. e(hr_ci_lo) and
    e(hr_ci_hi) are missing and the printed CI shows "not available" otherwise,
    and the postestimation positive-definite check honors e(v_pd) or a >1e-8
    floor. Mata numerical core unchanged, so point estimates are unaffected.{p_end}
{phang}2.3.2  30jun2026  Fixed e(rho) implied-correlation calculation. A
    variance term squared a macro holding a possibly-negative mean without
    parentheses, so `E_vT'^2 parsed as -(E_vT^2) when the mean was negative
    (exponentiation outranks unary minus), inflating the variance and deflating
    rho. Parenthesized the squared terms. Affects the reported e(rho) only;
    likelihood, coefficients, hazard ratio, and standard errors unchanged.
    Certification now passes 6/6.{p_end}
{phang}2.3.1  07may2026  Markout of covariates moved before data-contract
    asserts. e(N_persons) stored; BIC denominator switched to person count.
    Bivariate parameter-recovery cert added (hsmixture_certification_bivariate.do).{p_end}
{phang}2.3.0  04may2026  Updated description to reflect three-mode framing
    (factor(common), factor(separate), bivariate). V-fallback honest failure.{p_end}
{phang}2.2.1  02may2026  Convergence fix (vtol=1e-10, nrtol=1e-5,
    relative-gradient gate).{p_end}
{phang}2.2.0  22mar2026  Switched to Mata optimize(). Dropped factor variable support. Nobreak cleanup.{p_end}
{phang}2.0.0  21mar2026  Generalized for SSC. Renamed to treat/outcome.
    Added riskset(), level(). Fixed n_params computation.{p_end}
{phang}1.0.0  09feb2026  Initial release.{p_end}


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
Related: {help hsmixture_joint:hsmixture_joint} (one-factor under
{cmd:factor(common)} or 1-D locus under {cmd:factor(separate)}),
{help hsmixture:hsmixture} (single-equation),
{help hsmixture_joint_postestimation:hsmixture_joint_postestimation} (diagnostics)
