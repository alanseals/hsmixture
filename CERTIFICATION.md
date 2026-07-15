# Certification record

This file records which certification runs back each released version of
hsmixture. The certification logs themselves are deliberately untracked
(see `.gitignore`); this file is the durable, reviewable summary.

A version is **certified** only when all four scripts pass on the exact
commit being tagged. Passing runs on earlier code do not carry forward.
`RELEASING.md` is the checklist; `DEVELOPMENT.md` explains the reasoning
behind the current work.

When recording a run, copy the provenance block each script prints (date,
Stata version, OS, MP/cores, and the `which` output showing the resolved
`.ado` path and version) together with the git commit SHA that was tested.

## What an edit resets (scoping rule, recorded 2026-07-15)

If any `.ado`, any `_mata.do`, or a script's **own** `.do` file changes, the
affected boxes reset and those scripts must run again on the new code.

A `.sthlp` is documentation. Stata never reads it during estimation, and the
provenance block reports the resolved `.ado` starbang, not the help file. A
`.sthlp`-only edit therefore cannot change a number and does not reset
certification.

This supersedes the earlier wording, which reset every box on any `.sthlp`
change. That rule was written as a blunt guard against `.ado`/`_mata.do`
edits and swept in help text by accident. It is recorded here rather than
quietly dropped, because the 2026-07-15 documentation fix relied on it.

| Script | What it certifies | Criteria |
|---|---|---|
| `hsmixture_certification.do` | `hsmixture_joint` same-sign recovery of identified quantities, plus data-contract rejection tests across all three estimators including the `time()` false-rejection guards | 15 |
| `hsmixture_certification_common.do` | `factor(common)` recovery + `factor(separate)` falsifiability | 7 |
| `hsmixture_certification_opposite_signs.do` | opposite-signs recovery + an unscored `factor(common)` misspecification record | 7 |
| `hsmixture_certification_bivariate.do` | `hsmixture_bivariate` structure detection (not cell-wise recovery) | 6 |

`hsmixture_example.do` is a smoke test, not a certification. It proves the
commands run end-to-end and that the postestimation tools can consume their
output. It is not scored and its DGP is at the edge of identifiability by
design.

---

## v2.4.0 — all four scripts PASSED 2026-07-15 (commit SHA pending)

Code under test: working tree, uncommitted, ahead of `main` (`2c584b0`).
Every run on macOS, Stata 17 MP, 4 cores, batch mode from the package
directory.

| Script | Result | Retries | Time | Run start |
|---|---|---|---|---|
| `hsmixture_certification.do` | **PASSED 15/15** | 0 | 23 min | 2026-07-15 11:38:22 |
| `hsmixture_certification_common.do` | **PASSED 7/7** | 0 | 33 min | 2026-07-15 04:18:53 |
| `hsmixture_certification_opposite_signs.do` | **PASSED 7/7** | 0 | 31 min | 2026-07-15 12:00:58 |
| `hsmixture_certification_bivariate.do` | **PASSED 6/6** | 0 | 79 min | 2026-07-15 07:47:11 |

- [x] `hsmixture_certification.do` — 15/15
- [x] `hsmixture_certification_common.do` — 7/7
- [x] `hsmixture_certification_opposite_signs.do` — 7/7
- [x] `hsmixture_certification_bivariate.do` — 6/6

**Certified code: the commit tagged `v2.4.0`.** Resolve the SHA with
`git rev-parse v2.4.0^{commit}`. This file ships inside the commit it
certifies, so it cannot name that commit's own hash; the tag is the
identifier instead, and since a published tag is never moved
(`RELEASING.md` step 3) it is as durable as the hash.

**Why `main` and `opposite_signs` ran twice.** Both first passed earlier on
2026-07-15 (`main` at 04:00:55, `opposite_signs` at 07:16 and again at
09:06:16). The 2026-07-15 documentation fix then edited comments in both
scripts, so under the scoping rule above their logs no longer came from the
code on disk and both re-ran. The runs recorded in the table are the ones
that stand. `common` and `bivariate` were untouched by that fix and were not
re-run.

**One `.ado` was edited after the sweep, and it resets nothing.** The same fix
corrected a stale comment in `hsmixture_joint_postestimation.ado` that still
described the pre-v2.4.0 PD gate. `git diff` against `2c584b0` confirms only
comment lines and the starbang differ in that file; no executable line
changed. No certification script loads it (verified by grep across all four),
so no box it could affect exists. Recorded here because an `.ado` mtime later
than the sweep would otherwise look like a violation of the rule above.

### Provenance blocks

All four resolve to the package directory at the version being released. No
`.ado`, `_mata.do`, or `.pkg` file has changed since 2026-07-14 21:40:24,
which precedes every run above.

```
hsmixture_certification.do
  Run date/time:  15 Jul 2026 11:38:22
  Stata version:  17 (born 21 May 2024)
  OS / machine:   Unix / Mac (Apple Silicon)
  MP / cores:     1 / 4
  which hsmixture_joint -> ./hsmixture_joint.ado
                           *! version 2.4.0  14jul2026

hsmixture_certification_common.do
  Run date/time:  15 Jul 2026 04:18:53
  Stata version:  17 (born 21 May 2024)
  OS / machine:   Unix / Mac (Apple Silicon)
  MP / cores:     1 / 4
  which hsmixture_joint -> ./hsmixture_joint.ado
                           *! version 2.4.0  14jul2026

hsmixture_certification_opposite_signs.do
  Run date/time:  15 Jul 2026 12:00:58
  Stata version:  17 (born 21 May 2024)
  OS / machine:   Unix / Mac (Apple Silicon)
  MP / cores:     1 / 4
  which hsmixture_joint -> ./hsmixture_joint.ado
                           *! version 2.4.0  14jul2026

hsmixture_certification_bivariate.do
  Run date/time:  15 Jul 2026 07:47:11
  Stata version:  17 (born 21 May 2024)
  OS / machine:   Unix / Mac (Apple Silicon)
  MP / cores:     1 / 4
  which hsmixture_bivariate -> ./hsmixture_bivariate.ado
                               *! version 2.4.0  14jul2026
```

### Evidence recorded during this cycle

**The `v_pd` question is resolved.** The opposite-signs fit that failed strict
convergence on 2026-07-14 was measured on rerun at minimum eigenvalue
`1.85e-04` with `v_scaffold = 0` and `v_pd = 1`, passing 7/7. The earlier
`v_pd = 0` was the optimizer-abort defect posting a placeholder variance
matrix, not a threshold artifact. The bivariate fit independently measured
`1.86e-04` with `v_scaffold = 0` and `v_pd = 1` on a different DGP. The
scale-relative PD test shipped in v2.4.0 is correct on its own merits, but it
did not fix the symptom that motivated it; the abort retry did. The changelog
rationale that misstated this was corrected on 2026-07-15 (see
`DEVELOPMENT.md`).

**The abort retry is documented, not demonstrated.** No run in this cycle
tripped the optimizer abort (`retries = 0` on all four scripts), which is
expected: it is a knife-edge event under MP arithmetic variation, not a
reproducible one. Direct evidence of the failure it fixes exists in the
archived 2026-07-14 `r(430)` log. The recovery path is exercised by
inspection and documented in `DEVELOPMENT.md` item 2.

**The bivariate's configuration 5 was an ordinary local-mode failure, not an
abort.** Start 5 of 6 exited `converged = 0` at 361 iterations after repeated
"BFGS stepping has contracted, resetting BFGS Hessian" messages, at
log-likelihood `-27698.91`. This returns `_rc = 0` through Stata's normal
optimizer path and is not the numeric-derivative abort the retry handles,
which is why `retries = 0` is correct for that run. Five of six starts
converged to `-27691.7387`; the reported fit is the best of them
(configuration 2). Recorded because the resemblance is misleading in the raw
log.

**Reproducibility.** `hsmixture_certification_opposite_signs.do` ran three
times on identical code (07:16, 09:06, 12:00). All three produced identical
scored output (delta `1.0498`, SE `0.0871`, min pi `0.4625`), and the first
two logs are byte-identical apart from timestamps. This is a useful bound on
the MP nondeterminism documented in `DEVELOPMENT.md`, which is real but did
not surface on this DGP.

---

## v2.3.3 (02jul2026) — certified

All four certification scripts passed on Windows Stata on 2026-07-02 (logs
retained locally: same-sign 5/5, common 7/7, opposite-signs 7/7, bivariate
6/6). Those runs predate this provenance process, so the exact commit SHA was
not recorded in the logs; the code tested matches the numerical content of
tag `v2.3.3` (commit `c1ffed6`; the one later commit on `main`, `2c584b0`, is
documentation-only, verified by diff).

Known scope limits of the v2.3.3 certification, carried forward:

- The same-sign and `factor(common)` scripts scored the **signed** loadings,
  which are not identified. Those 2026-07-02 passes were therefore partly
  luck of which labeling the optimizer reached. Fixed in v2.4.0; see
  `DEVELOPMENT.md` item 5.
- The bivariate cert detects bivariate structure (nonzero off-diagonal mass,
  bounded rho) but does not certify cell-wise recovery of the 2x2
  probability matrix.
- The single-equation `hsmixture` has no dedicated parameter-recovery cert.
- K >= 3 is accepted syntactically but not certified against a known DGP.
