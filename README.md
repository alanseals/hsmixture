# hsmixture: Discrete-time hazard models with Heckman-Singer mixture heterogeneity

Stata package implementing three estimators for discrete-time duration models with nonparametric (Heckman-Singer 1984) unobserved heterogeneity.

- **`hsmixture`** -- Single-equation discrete-time hazard with K-type mixture
- **`hsmixture_joint`** -- Joint timing-of-events model (Abbring and van den Berg 2003), with `factor(common)` or `factor(separate)` heterogeneity structures
- **`hsmixture_bivariate`** -- Bivariate heterogeneity joint model (2x2 grid with free joint probabilities)
- **`hsmixture_joint_postestimation`** -- Convergence diagnostics and model comparison

## Installation

Install a tagged release directly from this GitHub repository. Pinning to a
tag gives a reproducible install, and every tag is a state that passed the
full certification suite (see `CERTIFICATION.md`); `main` may be ahead of the
last certified tag.

```stata
net install hsmixture, from("https://raw.githubusercontent.com/alanseals/hsmixture/v2.4.0") replace
```

Then verify with `which hsmixture`. To track the development tip instead,
replace `v2.4.0` with `main`.

To install from a local clone instead, copy the `.ado`, `.sthlp`, and `_mata.do` files in this directory to your Stata adopath (e.g., `~/ado/personal/h/` on macOS or `c:\ado\personal\h\` on Windows).

## Quick start

```stata
* Generate period dummies (factor variables not supported)
tab period, gen(pd_)
drop pd_1

* Single-equation hazard with K=2 mass points
hsmixture outcome_event x1 pd_*, id(id) k(2)

* Joint timing-of-events with separate loadings (default)
hsmixture_joint (treat_event = x1 pd_*) ///
    (outcome_event = x1 pd_*, treat(treated)) ///
    , id(id) k(2) riskset(treat_at_risk)

* Same model with one shared loading (manuscript-style one-factor MPH)
hsmixture_joint (treat_event = x1 pd_*) ///
    (outcome_event = x1 pd_*, treat(treated)) ///
    , id(id) k(2) factor(common) riskset(treat_at_risk)

* Postestimation diagnostics (always run before interpreting HRs)
hsmixture_joint_postestimation, all
```

## Heterogeneity factor structures (`factor()`)

`hsmixture_joint` supports two heterogeneity factor structures.

- **`factor(separate)`** (default). Two free signed loadings, `lambda_T` and `lambda_Y`. Per-type shifts `(lambda_T*v_k, lambda_Y*v_k)` lie on a 1-D locus through the origin, but the locus direction is data-determined. Opposite signs (negative selection between treatment-prone and outcome-prone latent types) are admissible.
- **`factor(common)`**. One free shared loading `lambda`. Per-type shifts `(lambda*v_k, lambda*v_k)` lie on the 45-degree line. This is the classical Heckman-Singer one-factor MPH model in which the same `exp(eta_i)` raises both hazards, with positive correlation between latent treatment and outcome propensity. This is the parameterization used in most Abbring-van den Berg identification arguments.

The two specifications nest. `factor(common)` is `factor(separate)` with the constraint `lambda_T = lambda_Y`. When the constraint is true the two HRs agree to numerical tolerance; when it is false `factor(separate)` is the more flexible model.

For data with substantial off-diagonal mass on the (v_T, v_Y) grid, neither joint variant is enough. Use `hsmixture_bivariate`, which estimates a free joint probability matrix on a 2x2 grid of corners.

## Data contract (validated, not assumed)

The estimators check the person-period panel before touching the likelihood and error on violations. The contract:

- **Absorbing outcome.** At most one outcome event per id, and a person's rows must stop at the event row. Rows after the event would enter the likelihood as spurious at-risk periods, so the commands refuse them.
- **One-time treatment (joint/bivariate).** At most one treatment event per id, and post-event person-periods must be excluded from the treatment equation via `riskset()`. `riskset()` is therefore effectively required: the commands error when it is omitted and any treated person remains observed after the treatment event (the normal shape of ToE data). The risk set must be 1 on the event row and 0 on every row after it; riskset = 0 before the event (delayed eligibility) is fine.
- **Treatment indicator (joint/bivariate).** `treat()` must be absorbing (never revert to 0) and must never equal 1 before the person's treatment event. Same-period and lagged switch-on conventions are both accepted; delta applies exactly where the indicator equals 1. Persons who enter already treated are accepted; persons who switch to treated with no event row in sample trigger a warning.
- **Ordering.** The row checks order person-periods by the current row order (keep the data sorted by id and time) or by an explicit `time(varname)` if you pass one. If the data are not sorted by the id variable the commands print a note, since row order then carries no time information and the checks cannot be trusted; pass `time()` in that case. `time()` does not change the likelihood, and the estimator is order-invariant within person, so this affects only the checks.

## When to use `riskset()`

For one-time treatment timing (e.g., program entry, first marriage, job displacement), supply a 0/1 indicator that equals 1 only when the person is *at risk* for the treatment event. Under the lagged treatment convention this is `treat_at_risk = (treated == 0)` where `treated` is the post-treatment dummy.

The package validates that `treat_event = 1` only occurs when `riskset = 1`, that no `riskset = 1` rows follow the treatment event, and that each person has at most one treatment event and at most one outcome event in the estimation sample.

## Multistart

`hsmixture_joint` and `hsmixture_bivariate` use multiple starting configurations and select the best log-likelihood. Mixture surfaces are multimodal at any K, so single-start is unsafe. Defaults are 7 starts (separate), 6 starts (common or bivariate). Override with `nstarts()`.

`hsmixture` (single-equation) runs a single start from the GLM coefficients. Users running `hsmixture` alone should compare results against runs with user-supplied starting values via `from()` to guard against local modes.

## Convergence

All three estimators expose a strict-convergence flag.

- `e(converged) == 1` requires (a) BFGS converged on its own criterion, (b) relative gradient `|grad|/(1+|LL|) < 1e-5`, and (c) variance matrix positive definite.
- `e(converged_bfgs)` is BFGS's own flag alone.
- `e(grad_norm)`, `e(rel_grad)`, `e(v_pd)` expose the components.

Always gate model comparison and HR interpretation on `e(converged) == 1`. Information criteria from a non-converged fit are not interpretable.

## Information criteria denominator

AIC and BIC are computed using `e(N_persons)`, not row count `e(N)`. The IID unit in this mixture model is the person, since the latent type is integrated out at the person level. Using the person-period count would inflate the BIC penalty. `e(N)` follows Stata's row-count convention and is preserved for compatibility; `e(N_persons)` is the appropriate denominator for cross-K comparisons. The package's Display routines and `hsmixture_joint_postestimation, compare` both use `e(N_persons)`.

## Requirements

- Stata 14 or later
- Data in person-period (long) format, sorted by id and time (or pass `time()`)
- Covariates as explicit variables (no factor variable notation)
- The data contract above. In short: rows stop at the outcome event; `riskset()` supplied for the joint estimators; `treat()` absorbing and never ahead of its event. A fully censored id with no outcome event is allowed. The single-equation `hsmixture` takes no treatment variable and requires only the outcome contract.

## Documentation

After installation, type `help hsmixture`, `help hsmixture_joint`, `help hsmixture_bivariate`, or `help hsmixture_joint_postestimation` in Stata for full documentation including syntax, methodology, and examples.

## Examples and certification

- `hsmixture_example.do` -- smoke test that exercises all three estimators end-to-end on a small synthetic DGP. Not a parameter-recovery test (the DGP is at the edge of identifiability for K=2).
- `hsmixture_certification.do` -- parameter-recovery test for `hsmixture_joint` with same-sign loadings.
- `hsmixture_certification_common.do` -- parameter-recovery test for `factor(common)`. Includes a falsifiability test that refits `factor(separate)` on common-loading data and checks that no spurious gap is manufactured.
- `hsmixture_certification_opposite_signs.do` -- parameter-recovery test for the opposite-signs heterogeneity case (negative selection between treatment-prone and outcome-prone types).
- `hsmixture_certification_bivariate.do` -- parameter-recovery test for `hsmixture_bivariate` on a DGP with substantial off-diagonal mass on the (v_T, v_Y) grid.

Each cert script asserts on a recovery threshold and exits non-zero if the criteria fail. `hsmixture_certification.do` additionally asserts that the data-contract validations reject malformed panels (post-outcome rows, omitted riskset, reverting or mistimed treatment indicators). Each script opens with a provenance block (date, Stata version, OS, resolved .ado path and version).

`CERTIFICATION.md` records which certification runs back each released version, including the commit tested. `RELEASING.md` is the release checklist. `DEVELOPMENT.md` holds maintainer notes: why the current version changed what it did, known platform gotchas, and open items.

## Citation

If you use this package, please cite:

> Park, J. and R.A. Seals. 2026. hsmixture: Stata module for discrete-time hazard models with Heckman-Singer mixture heterogeneity.

## References

- Abbring, J.H. and G.J. van den Berg. 2003. The nonparametric identification of treatment effects in duration models. *Econometrica* 71(5): 1491-1517.
- Heckman, J.J. and B. Singer. 1984. A method for minimizing the impact of distributional assumptions in econometric models for duration data. *Econometrica* 52(2): 271-320.
- Jenkins, S.P. 1995. Easy estimation methods for discrete-time duration models. *Oxford Bulletin of Economics and Statistics* 57(1): 129-138.

## License

Released under the MIT License. See the `LICENSE` file for the full text.

## Authors

Jonghoon Park and R. Alan Seals
Department of Economics, Auburn University
