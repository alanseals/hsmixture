# Release checklist

There is no hosted Stata CI runner, so certification is a manual gate.
This checklist is the release process; do not skip steps. A release that
skips certification is not a release.

## 1. Pre-flight (any machine)

- [ ] Working tree clean except the intended release changes; `git status`.
- [ ] Version stamps agree everywhere: the `*! version X.Y.Z  DDmonYYYY`
      starbang in all four `.ado` files, the `{* *! version ...}` header in
      all four `.sthlp` files, and `Distribution-Date` in `hsmixture.pkg`.
      `grep -n "version 2\." *.ado *.sthlp hsmixture.pkg` and read the output.
- [ ] Every functional change has a version-history entry in the affected
      `.sthlp` files. "No functional change to this command" only when true.
- [ ] `hsmixture.pkg` lists every file a `net install` must deliver
      (all `.ado`, all `_mata.do`, all `.sthlp`, example + cert scripts,
      LICENSE).
- [ ] The README install line names the version about to be released.
      It goes in the release COMMIT, because the tag is cut from that
      commit and must contain a README that points at itself. Until the
      tag exists the line refers to a tag that 404s; that is expected in a
      pre-release working tree and is the reason step 3 pushes the commit
      and the tag together, and the reason `main` is never pushed alone.
- [ ] No private paths, credentials, or unpublished estimates in any
      tracked file.

## 2. Certification (any licensed Stata)

Run all four scripts from the package directory on the EXACT code being
released. Two rules, both learned by breaking them:

- **Freeze the tree for the duration of the sweep.** Stata reads do-files
  incrementally, so editing a script while it is running shifts bytes under
  the reader and kills the run with a parse error partway through.
- **Read each provenance block before trusting a pass.** It prints the
  resolved `.ado` path and version. If it does not name this package
  directory and the version being released, the run tested the wrong code
  and the pass is worthless.

Budget several hours. Certification runtime is dominated by the estimator
(~10 s/iteration on 45-80k rows), not by the tests. See `DEVELOPMENT.md`.

- [ ] `do hsmixture_certification.do` — PASSED (15/15)
- [ ] `do hsmixture_certification_common.do` — PASSED (7/7)
- [ ] `do hsmixture_certification_opposite_signs.do` — PASSED (7/7)
- [ ] `do hsmixture_certification_bivariate.do` — PASSED (6/6)

Each script exits non-zero on failure and prints a provenance block at the
top of its log. If any script fails, fix the code, and restart this
checklist from step 1.

- [ ] Copy the provenance blocks and pass counts into `CERTIFICATION.md`
      under the new version's entry and check its boxes. Identify the
      certified code by the tag about to be cut, not by a SHA: the file ships
      inside the commit it certifies and cannot name that commit's own hash.
      The tag is never moved, so it resolves the SHA durably
      (`git rev-parse vX.Y.Z^{commit}`).

## 3. Commit, tag, release

- [ ] Commit the release (including the updated `CERTIFICATION.md` and the
      README install line from step 1).
- [ ] Tag: `git tag vX.Y.Z` on that commit. Never move or rewrite a
      published tag; a follow-up fix gets a new patch tag.
- [ ] Push commit and tag together: `git push origin main --tags` (no force
      flags, ever). Pushing `main` without the tag publishes a README whose
      install command 404s.
- [ ] Create a GitHub Release on the tag:
      `gh release create vX.Y.Z --title "hsmixture vX.Y.Z" --notes-file <notes>`,
      with the version-history entry and the certification summary as notes.
- [ ] Confirm the install path resolves before telling anyone it works:
      `curl -sI https://raw.githubusercontent.com/alanseals/hsmixture/vX.Y.Z/hsmixture.pkg`
      must return 200.

## If certification fails

Do not push. The working tree carries version stamps for a release that
does not exist yet; either fix the code and restart at step 1, or revert
the stamps (`.ado`/`.sthlp` starbangs, `.pkg` Distribution-Date, README
install line, `CERTIFICATION.md` entry) back to the last certified
version. A pre-release working tree is not a release.

## 4. After release

- [ ] Re-sync any production copies (e.g., a project's `ado/` directory)
      from the tagged commit, byte-identical, and archive the previous
      copies rather than deleting them.
- [ ] Sanity-check the public install path in a fresh Stata session:
      `net install hsmixture, from("https://raw.githubusercontent.com/alanseals/hsmixture/vX.Y.Z") replace`
      then `which hsmixture` and confirm the starbang version.
