{smcl}
{* *! version 2.3.2  30jun2026}{...}
{vieweralsosee "hsmixture_joint" "help hsmixture_joint"}{...}
{vieweralsosee "hsmixture_bivariate" "help hsmixture_bivariate"}{...}
{viewerjumpto "Syntax" "hsmixture_joint_postestimation##syntax"}{...}
{viewerjumpto "Description" "hsmixture_joint_postestimation##description"}{...}
{viewerjumpto "Options" "hsmixture_joint_postestimation##options"}{...}
{viewerjumpto "Examples" "hsmixture_joint_postestimation##examples"}{...}
{title:Title}

{p2colset 5 40 42 2}{...}
{p2col:{cmd:hsmixture_joint_postestimation} {hline 2}}Postestimation diagnostics for hsmixture_joint and hsmixture_bivariate{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:hsmixture_joint_postestimation} [{cmd:,} {it:options}]


{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt conv:ergence}}display convergence diagnostics{p_end}
{synopt:{opt post:erior}}display posterior type probability information{p_end}
{synopt:{opt comp:are}}display model comparison statistics{p_end}
{synopt:{opt lrtest(name)}}LR test against stored estimates{p_end}
{synopt:{opt all}}display all diagnostics{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:hsmixture_joint_postestimation} provides diagnostic tests and model comparison
tools after estimation with {cmd:hsmixture_joint} or {cmd:hsmixture_bivariate}.
These diagnostics help assess:

{phang2}1. Whether the model converged properly{p_end}
{phang2}2. Whether parameters are at or near boundary values{p_end}
{phang2}3. Whether the number of mass points K is appropriate{p_end}
{phang2}4. How the current model compares to alternatives{p_end}


{marker options}{...}
{title:Options}

{phang}
{opt convergence} displays convergence diagnostics including:

{pmore}- Convergence status (converged or not){p_end}
{pmore}- Number of iterations{p_end}
{pmore}- Relative gradient norm |grad|/(1+|LL|) (should be < 1e-5 at solution){p_end}
{pmore}- Variance matrix check (should be positive definite){p_end}
{pmore}- Parameter bounds check (|lambda| < 100, 0 < pi < 1){p_end}
{pmore}- Mass point separation check{p_end}

{phang}
{opt posterior} displays the estimated mixture distribution: the {it:prior}
(unconditional) population type shares pi_k with each type's mass point and
factor-implied hazard shifts. Note this reports prior population shares, not
per-person {it:posterior} type probabilities; no per-observation Bayes-rule
posterior is computed. The option name is retained for back-compatibility.

{phang}
{opt compare} displays AIC, BIC, and guidance on model selection.

{phang}
{opt lrtest(name)} performs a likelihood ratio test comparing the current
model to the stored estimates in {it:name}. Use this to test K=3 vs K=2, etc.

{phang}
{opt all} displays all available diagnostics.


{marker examples}{...}
{title:Examples}

{pstd}{bf:Basic convergence check}{p_end}
{phang2}{cmd:. hsmixture_joint (...), id(id) k(2)}{p_end}
{phang2}{cmd:. hsmixture_joint_postestimation, convergence}{p_end}

{pstd}{bf:All diagnostics}{p_end}
{phang2}{cmd:. hsmixture_joint_postestimation, all}{p_end}

{pstd}{bf:LR test for K selection (both fits must strictly converge)}{p_end}
{phang2}{cmd:. hsmixture_joint (...), id(id) k(2)}{p_end}
{phang2}{cmd:. assert e(converged) == 1}{p_end}
{phang2}{cmd:. estimates store m_k2}{p_end}
{phang2}{cmd:. hsmixture_joint (...), id(id) k(3)}{p_end}
{phang2}{cmd:. assert e(converged) == 1}{p_end}
{phang2}{cmd:. hsmixture_joint_postestimation, lrtest(m_k2)}{p_end}

{pstd}{bf:Systematic K selection (gated on strict convergence)}{p_end}
{phang2}{cmd:. local model_list}{p_end}
{phang2}{cmd:. forvalues k = 2/4 {c -(}}{p_end}
{phang2}{cmd:.     capture noisily hsmixture_joint (...), id(id) k(`k') iterate(200)}{p_end}
{phang2}{cmd:.     local rc = _rc}{p_end}
{phang2}{cmd:.     if `rc' == 0 {c -(}}{p_end}
{phang2}{cmd:.         hsmixture_joint_postestimation, convergence}{p_end}
{phang2}{cmd:.         if e(converged) == 1 {c -(}}{p_end}
{phang2}{cmd:.             estimates store m_k`k'}{p_end}
{phang2}{cmd:.             local model_list `model_list' m_k`k'}{p_end}
{phang2}{cmd:.         {c )-}}{p_end}
{phang2}{cmd:.         else {c -(}}{p_end}
{phang2}{cmd:.             display as err "K=`k' did not strictly converge; not stored."}{p_end}
{phang2}{cmd:.         {c )-}}{p_end}
{phang2}{cmd:.     {c )-}}{p_end}
{phang2}{cmd:.     else {c -(}}{p_end}
{phang2}{cmd:.         display as err "K=`k' estimation failed (rc=`rc'); not stored."}{p_end}
{phang2}{cmd:.     {c )-}}{p_end}
{phang2}{cmd:. {c )-}}{p_end}
{phang2}{cmd:. if "`model_list'" != "" {c -(}}{p_end}
{phang2}{cmd:.     estimates stats `model_list'}{p_end}
{phang2}{cmd:. {c )-}}{p_end}
{phang2}{cmd:. else {c -(}}{p_end}
{phang2}{cmd:.     display as err "No models strictly converged; nothing to compare."}{p_end}
{phang2}{cmd:. {c )-}}{p_end}

{pstd}{it:Note on BIC.} {cmd:estimates stats} computes BIC from the row count
{cmd:e(N)} (person-periods), whereas this package reports BIC on the person
count {cmd:e(N_persons)}. The two can rank K differently because the larger row
count inflates the per-parameter penalty. For the person-count BIC, run
{cmd:hsmixture_joint_postestimation, compare} after each fit and read its BIC
line.{p_end}

{pstd}{bf:After hsmixture_bivariate}{p_end}
{phang2}{cmd:. hsmixture_bivariate (...), id(id)}{p_end}
{phang2}{cmd:. hsmixture_joint_postestimation, all}{p_end}


{title:Interpreting Diagnostics}

{pstd}
{bf:Convergence Problems:}

{pmore}
If the model did not converge, try: (1) increase {opt iterate()},
(2) provide better starting values with {opt from()}, or (3) reduce K.

{pstd}
{bf:Boundary Parameters:}

{pmore}
If lambda approaches 0, heterogeneity may not be needed. If pi approaches 0 or 1,
some types have negligible population share, suggesting K is too large.

{pstd}
{bf:Mass Point Separation:}

{pmore}
If mass points v_k are nearly equal (|v_k - v_j| < 0.1), the types are not
well-distinguished. This typically indicates K is too large for the data.

{pstd}
{bf:Model Selection:}

{pmore}
AIC and BIC are reported for each fit. Restrict comparisons to models that
strictly converged ({cmd:e(converged) == 1}); information criteria from a
non-converged fit are not interpretable.


{title:Stored results}

{pstd}
{cmd:hsmixture_joint_postestimation} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{synopt:{cmd:r(converged)}}1 if model converged{p_end}
{synopt:{cmd:r(grad_norm)}}gradient norm at solution{p_end}
{synopt:{cmd:r(aic)}}Akaike information criterion{p_end}
{synopt:{cmd:r(bic)}}Bayesian information criterion{p_end}
{synopt:{cmd:r(lr_stat)}}LR test statistic (if lrtest specified){p_end}
{synopt:{cmd:r(lr_df)}}LR test degrees of freedom{p_end}
{synopt:{cmd:r(lr_pval)}}LR test p-value{p_end}


{title:Version history}

{phang}2.3.2  30jun2026  Coordinated package version bump. No functional change
    to this command; see hsmixture_bivariate for an e(rho) calculation fix.{p_end}
{phang}2.3.1  07may2026  lrtest() path now uses _estimates hold/unhold to
    preserve the user's active e() (the previous restore . was non-standard
    syntax). BIC denominator switched to person count via e(N_persons), with
    fallback to e(N) for pre-v2.3.1 fits.{p_end}
{phang}2.3.0  05may2026  Added factor() awareness for hsmixture_joint output;
    common-loading display rows.{p_end}
{phang}2.2.1  02may2026  Strict-convergence diagnostic uses relative gradient
    |grad|/(1+|LL|) < 1e-5; spike-and-slab corner detection.{p_end}
{phang}2.0.0  21mar2026  Added support for hsmixture_bivariate. Generalized for SSC.{p_end}
{phang}1.0.0  17jan2026  Initial release.{p_end}


{title:Authors}

{pstd}
Jonghoon Park and R. Alan Seals{break}
Department of Economics{break}
Auburn University{break}
{browse "mailto:jzp0200@auburn.edu":jzp0200@auburn.edu} /
{browse "mailto:ras0029@auburn.edu":ras0029@auburn.edu}


{title:Also see}

{psee}
Estimation: {help hsmixture_joint:hsmixture_joint},
{help hsmixture_bivariate:hsmixture_bivariate}
