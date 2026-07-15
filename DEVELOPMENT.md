# Development notes

Maintainer-facing notes for the hsmixture package. `README.md` is for users,
`CERTIFICATION.md` records what has been certified and on what code, and
`RELEASING.md` is the release checklist. This file explains the reasoning
behind the current work so a maintainer picking the package up cold does not
have to rediscover it.

Version history per command lives in the `.sthlp` files. This file covers
*why*, not *what*.

---

## Current state (2026-07-15)

**Working tree is v2.4.0 and is NOT released.** `main` on GitHub is at v2.3.3
(commit `2c584b0`), tag `v2.3.3` at `c1ffed6`. Every v2.4.0 change is
uncommitted. All four certification scripts passed on 2026-07-15; the
changelog-rationale fix that day re-ran the two scripts it touched. See
`CERTIFICATION.md` for the scoreboard, the provenance blocks, and the rule
governing what a given edit resets.

Do not commit, tag, or push until `CERTIFICATION.md` shows all four scripts
passing on the exact code being released. `README.md`'s install line already
names the `v2.4.0` tag, which does not exist yet; that is expected in a
pre-release tree and is why the commit and the tag must be pushed together
(`RELEASING.md` step 3).

---

## What v2.4.0 changes, and why

v2.4.0 came out of an external audit (2026-07-14) that found the estimators
would silently accept person-period panels violating the assumptions their
likelihood is built on. Certifying the fixes then surfaced two further
defects that only appear at runtime. All five items below are independent.

### 1. Data-contract validation (the reason for the release)

The likelihood sums a contribution over **every** surviving estimation row.
It therefore assumes, without ever checking, that:

- the outcome is absorbing, so a person's rows stop at the outcome event; and
- treatment is one-time, so post-event rows are excluded from the treatment
  equation via `riskset()`; and
- `treat()` is an absorbing indicator that never precedes its own event.

A panel violating any of these produced a quietly biased fit rather than an
error. v2.4.0 validates all three and errors out. `riskset()` is now
effectively required for the joint estimators, because without it every row
enters the treatment equation as at-risk.

New `time(varname)` supplies explicit within-person ordering for these checks.
When it is omitted the checks read physical row order, which is what the
documented contract (sorted by id and time) delivers; the commands print a
note when the data are not sorted by the id variable, since row order then
carries no time information and the checks cannot be trusted.

**These checks changed no estimate.** They reject data that was always
invalid. The certification DGPs and the empirical pipeline that motivated the
package both satisfy the contract already.

### 2. Optimizer-abort recovery

Mata's `optimize()` can abort outright ("could not calculate numerical
derivatives -- discontinuous region with missing values encountered") when a
numeric-derivative probe degenerates near a flat optimum. Before v2.4.0 a
single abort killed the **entire multistart**: observed 2026-07-14, config 1
aborted and configs 2-7 then died instantly without a single iteration,
returning `r(430)` and no results.

The abort is not a property of the likelihood, which is finite everywhere
(overflowed linear predictors clip to a finite range; `ln(0)` floors at
-700; Mata's `min`/`max`/`rowmax`/`rowmin` ignore missing rather than
propagate it -- verified directly). It is a knife-edge event in Stata's
derivative machinery: see "MP nondeterminism" below.

v2.4.0 captures the return code per start. On an abort it resets the Mata
cache and retries **that start once** from a deterministically jittered
vector (`b0 * 1.001` -- enough to put BFGS on a different arithmetic path,
not enough to change which basin the start explores). A second abort marks
the start failed and the multistart continues. `e(n_aborted)` counts the
retries. The single-equation command has one start, so it gets the same
retry, where an abort would otherwise kill the whole fit.

This is defensive code that cannot be triggered on demand. It is documented
rather than demonstrated; see "Known gotchas" for why.

### 3. Scale-relative positive-definiteness test

`e(v_pd)` (and therefore `e(converged)`) previously required the smallest
eigenvalue of `e(V)` to exceed an **absolute** `1e-8`. That threshold was
scale-dependent, and it conflated two different questions: is this matrix
numerically positive definite, and is this matrix the `I * 1e-20` placeholder
posted when Hessian inversion fails.

v2.4.0 separates them, as the audit recommended. `e(v_pd)` now requires the
smallest eigenvalue to be strictly positive **and** to exceed `1e-12` of the
largest -- scale-free, with roughly 100x margin over eigenvalue rounding
error at this problem size. The placeholder case is excluded explicitly by
`e(v_scaffold)`, which is the only thing rejecting it now (the placeholder is
perfectly conditioned, so a purely relative test would otherwise pass it).

**Correction worth recording.** This change was originally justified by a
belief that an eigenvalue was hovering near the `1e-8` floor and MP noise was
pushing it across. That was wrong. Measured minimum eigenvalue on the fit in
question is `1.85e-04`, four orders of magnitude clear of the old floor. The
`v_pd = 0` verdict that prompted the investigation was the abort bug in
disguise: when the optimizer aborts, the `.ado` posts the placeholder V and
correctly reports `v_pd = 0`. One root cause, not two. The relative test is
still correct on its own merits and still what the audit asked for, but it
did not fix the symptom that motivated it -- item 2 did.

### 4. Numerical hardening and diagnostics

- Mixture probabilities use a max-shifted softmax in all three likelihoods
  and in every `e()`-posting block. Algebraically identical; values agree to
  floating-point rounding (verified by Monte Carlo: max relative difference
  ~7e-16, roughly 13 orders of magnitude below the tightest certification
  margin). Overflow-safe for large mixture logits.
- `e(clip_hits)` counts linear-predictor evaluations at the `[-20, 10]`
  numerical bounds at the optimum. Zero means the fitted surface is the exact
  cloglog likelihood everywhere. The estimators print a caution when nonzero.
- `e(v_mineig)`, `e(v_scaffold)` expose the inputs to the `v_pd` verdict, so
  a `v_pd = 0` is diagnosable instead of opaque.
- `e(n_finite)`, `e(n_bfgs_conv)`, `e(n_aborted)` report per-start outcomes.
  The multistart log prints the counts, and the "no usable result" error now
  says that plainly instead of "all starting configurations failed to
  converge", which fired only when no start returned a finite likelihood.

Start selection is unchanged and deliberate: the best **finite**
log-likelihood wins regardless of its convergence flag. A higher-likelihood
point is the better MLE candidate even when BFGS stopped short there;
preferring a converged-but-lower-likelihood start would report a local mode
as the optimum. The selected start's status is reported honestly via
`e(converged)` / `e(converged_bfgs)`.

### 5. What the certification scripts score, and why

The K=2 mixture's type labels are **not identified**. The likelihood is
exactly invariant under

    (pi_1 <-> pi_2, lambda_T -> -lambda_T, lambda_Y -> -lambda_Y,
     cons_T -> cons_T + lambda_T, cons_Y -> cons_Y + lambda_Y)

because `v_1 = 0` / `v_2 = 1` fixes the mass-point **scale** but not which
type carries the shift. Verified numerically against a replica of the Mata
likelihood: log-likelihood difference exactly 0.0. Both labelings are the
MLE, and the optimizer can terminate at either depending on starting values
and the platform's floating-point path.

Consequently the certification scripts score only identified quantities:
delta, loading **magnitudes**, the **relative** sign of the two loadings (the
swap flips both together, so same-sign vs opposite-sign structure survives),
and the **unordered** mixture shares. The opposite-signs script always did
this; the same-sign and `factor(common)` scripts did not, and passed for
years only because the optimizer happened to reach the non-swapped labeling.
Both were fixed on 2026-07-14 after a run reached the mirror labeling and
"failed" on a criterion that tested a coin flip.

The same-sign script gained a criterion in the process (relative sign), so
this made it stricter, not looser.

**Rule for maintainers: never score, report, or compare across runs a signed
loading, an individual mass point, or a specific `pi_k`.** Only delta/HR,
magnitudes, relative signs, and unordered shares are comparable. This applies
to downstream users' tables as much as to the certification scripts.

---

## Known gotchas

### MP nondeterminism

Under Stata/MP the order of parallel arithmetic reductions varies run to run,
so results can differ in the last bits between otherwise identical runs. This
is normally invisible. It is not invisible here, because mixture likelihoods
have flat optima where a last-bit difference decides which path BFGS takes.

Observed on identical code and identical data (2026-07-14, 4-core Mac): the
opposite-signs certification's first start completed in 73 iterations, then
aborted at iteration 62, then completed in ~85 iterations, on three
consecutive runs. All four certification scripts print `MP / cores` in their
provenance block for this reason.

Practical consequences:

- A run that aborts is not evidence of a bug in the data or the model.
- Which labeling the optimizer lands on is not stable across platforms.
- Do not expect bit-identical reproduction across machines or MP settings.

### The unscored misspecification demo is expensive

`hsmixture_certification_opposite_signs.do` Part 3b deliberately fits
`factor(common)` to opposite-signs data to document what misspecification
looks like. It is informational and **not a scored criterion**.

That fit stalls on a flat region where each BFGS iteration degenerates into
repeated failed line searches plus a full numeric gradient. Measured at ~48
seconds per iteration on a 4-core Mac -- five times the pace of the scored
fit. At its former budget (6 default starts x 200 iterations) it could run
**16 hours**, dwarfing the rest of the suite. It is now capped at
`nstarts(1) iterate(30)`, which is ample: the diagnostic signature is that
the fit does not converge and delta collapses toward the naive cloglog, and
that is visible immediately.

A 2026-07-02 Windows log recorded this fit "stopping on its own near
iteration 55". That was the abort ending it -- a machine-dependent accident,
not a property to rely on. Hence the explicit cap.

`hsmixture_example.do`'s K=3 and bivariate fits have the same character on
its small DGP (documented in that file's header as expected not to converge)
and are capped at `iterate(100)` for the same reason.

### Certification runtime is dominated by the estimator, not the tests

The scored fits run ~10 seconds per iteration on 45-80k rows because the d0
evaluator computes ~46 numeric derivatives per iteration. Budget roughly
20-40 minutes per certification script, longer for the opposite-signs and
bivariate scripts. This is inherent to the design and is not something to
"fix" under release pressure.

---

## Running the certification suite

From the package directory, so that `.` resolves the working-tree `.ado`
files ahead of any installed copy (Stata searches `.` before PERSONAL and
PLUS):

```stata
cd <package directory>
do hsmixture_certification.do                  // 15/15
do hsmixture_certification_common.do           //  7/7
do hsmixture_certification_opposite_signs.do   //  7/7
do hsmixture_certification_bivariate.do        //  6/6
do hsmixture_example.do                        // smoke test, not a criterion
```

Each script exits non-zero on failure and opens with a provenance block
(date, Stata version, OS, MP/cores, and `which` output showing the resolved
`.ado` path and version). **Read the provenance block before trusting a
pass**: if it does not name the package directory and the expected version,
the run tested the wrong code.

Two operational rules, both learned the hard way:

1. **Never edit a `.do` file while Stata is running it.** Stata reads
   do-files incrementally; an edit shifts bytes under the reader and the run
   dies with a parse error (`r(133)`) partway through. Freeze the tree for
   the duration of a sweep.
2. **Run the certifications before the example.** The example's capped
   corner fits are slow and prove the least; putting them first delays every
   verdict that matters.

---

## Open items

- **Finish v2.4.0 certification** and fill in `CERTIFICATION.md`. See that
  file for the live status.
- **Correct the changelog rationale for the PD test — DONE 2026-07-15.**
  `hsmixture_joint.sthlp`, `hsmixture_certification.do`, and
  `hsmixture_certification_opposite_signs.do` attributed the `v_pd = 0`
  verdict to an eigenvalue near the absolute threshold. All three now give
  the real reason from item 3: the floor was scale-dependent and conflated a
  genuinely ill-conditioned V with the `I*1e-20` placeholder. The
  opposite-signs comment additionally claimed the `v_pd = 0` fit had "real
  standard errors"; it did not, because the posted V was the placeholder, and
  that sentence is gone. `hsmixture.sthlp` and `hsmixture_bivariate.sthlp`
  needed no change — they state only what the test became, not why, and both
  defer to `hsmixture_joint`. The two edited certification scripts were
  re-run afterward; see `CERTIFICATION.md` for the scoping rule and the
  resulting logs.
- **`hsmixture_joint_postestimation.ado` and the `1e-8` floor — examined
  2026-07-15, no code change needed.** The earlier note here implied the
  command recomputes a PD verdict against the stale absolute floor. It does
  not. It reads `e(v_pd)` when present, which every estimator from v2.3.3 on
  posts, so no fit from a current version can disagree with the
  estimation-time verdict. The `1e-8` recompute is reachable only for
  estimates stored before `e(v_pd)` existed, and for those it applies the
  floor those versions were produced under, which is the right answer for
  them. What was actually wrong was the comment above it, which still
  described the pre-v2.4.0 gate; corrected 2026-07-15. Left as a deliberate
  legacy path. `hsmixture_joint_postestimation.sthlp`'s "no functional change
  to this command" remains true.
- **Base `hsmixture` has no parameter-recovery certification.** It is
  exercised end-to-end by the example and by the contract-rejection tests
  (C5), but no script asserts that it recovers a known DGP. The joint certs
  do not transfer: `hsmixture_mata.do` is an independent hand-coded
  likelihood.
- **The bivariate certification detects structure, not recovery.** Its
  criteria check for nonzero off-diagonal mass and a bounded rho; they do not
  assert cell-wise recovery of the 2x2 probability matrix, and the label swap
  flips which cells are "off-diagonal". Do not cite it as a recovery
  certificate.
- **K >= 3 is accepted but not certified** against a known DGP and is prone
  to spike-and-slab corners.
- **No automated CI.** Certification requires a licensed Stata, so the gate
  is the manual checklist in `RELEASING.md`. A workflow that only checked
  file presence would be theater and should not be added.
