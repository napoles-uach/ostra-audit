# Reporting protocol: from a finding to an upstream issue

This expands step 6 of [`PROTOCOL.md`](PROTOCOL.md). It describes how a verified finding becomes an issue (or a
comment, or a pull request) in the upstream project, so that any contributor can follow the same steps and every
report has the same quality. It was distilled from the first reports (libxc#646, pyscf#3669, multicharge#72,
pyscf#3674, libcint).

**Golden rule:** agents draft; a human reviews and files. No agent opens issues, comments or pull requests upstream.

---

## 1. Gate: is the finding ready to report?

All of the following must hold. If any is missing, the finding stays in the repository.

- [ ] **Status** is at least `attributed`: minimal reproducer, exact source line, and a minimal change that makes
      the discrepancy disappear (three evidence levels: anomaly reproduced → cause identified → defect confirmed).
- [ ] **Impact** is measured (level I–III) or self-evident (crash, out-of-bounds access, silently wrong output).
      Level IV findings are not reported; they are kept for the paper.
- [ ] **Independent verification:** a second agent or person, who did not write the original harness, reproduced
      it with their own test (see "Verification by a second agent" in `PROTOCOL.md`). Verification is done
      **blind**: the verifier receives the location and the claim, not our numbers or patch; the patch is reviewed
      in a second phase. Record the result in `verification/`.
- [ ] **Not design and not documented:** checked against the documentation and the original paper (including
      errata). Deliberate choices that the authors very likely know about (e.g. HSE s-scaling) are not reported.
- [ ] **Not already known:** searched open *and closed* issues, pull requests, changelog and recent commits of
      upstream (search by file name, function name and symptom).
- [ ] **Still present on the latest code:** re-checked on the current default branch (record the commit hash),
      not only on the release where it was found.

## 2. Choose where and how to report

- **Report where the defect lives.** If library A bundles or wraps library B and the bug is in B, report to B; to A
  only if A's copy/build is what is broken (e.g. PySCF's bundled libcint vs upstream libcint) or A exposes it.
- **One issue per repository and per root cause.** Related findings in the same repository with the same cause are
  consolidated (e.g. the gCP tables in simple-dftd3). Different causes → separate issues, cross-referenced.
- **Unmaintained projects:** a single consolidated issue; consider also notifying the main downstream users.
- **Existing issue on the same topic:** add a comment with the new evidence instead of opening a new issue
  (e.g. multicharge#72).
- **Issue vs pull request:** default to an issue with the patch attached. Open a PR only if the maintainers'
  guidelines prefer it or they ask for it.
- **Order:** severity and silence first (wrong results without warning > crash > consistency > niche), and at most
  one or two open reports per project at a time, so as not to overload maintainers.

## 3. Prepare the `report/` folder

Inside `findings/OSTRA-NNN-*/report/`:

| File | Content |
|---|---|
| `issue.md` | The issue body, ready to paste |
| `issue_template_format.md` | The same text rearranged into the project's issue template, if it has one (`.github/ISSUE_TEMPLATE/`) |
| `*_repro.py` (or `.c`) | Self-contained minimal reproducer(s), runnable with only the library installed |
| `*.patch` | Minimal patch against the current default branch (`git diff` format, verified to apply) |

### Structure of the issue body

1. **Title**: what fails, where, under which condition, consequence. Example:
   *"Out-of-bounds read of TURNOVER_POINT in fmt.c for Boys order m >= 40 corrupts SR integrals with high angular momentum"*.
2. **Summary** (2–4 sentences): the defect and who is affected.
3. **Where**: file, function, lines, with a permalink to the commit checked.
4. **Reproducer**: short script (collapsible `<details>` if long) and its expected vs actual output.
5. **Evidence**: one compact table of numbers against a reference (analytic, mpmath, finite differences,
   another implementation). Say what the reference is.
6. **Cause**: one paragraph explaining *why*, pointing to the line.
7. **Impact**: which inputs/methods are affected and by how much; say also what is *not* affected.
8. **Suggested fix**: the patch, plus how it was verified: the failing case now passes, unaffected cases are
   bit-identical, and a physics check (e.g. analytic gradient vs finite differences) passes.
9. **Environment**: library versions and commit, compiler, OS **where the tests were actually run**
   (if it was a Linux container, say so; if confirmed on another platform, add it).
10. **Signature** (below).

### Tone

- Factual and neutral. Describe the code, not the authors. Avoid "careless", "wrong design", "should have".
  Prefer "fragile numerical decision", "the table is not extended for…", "the threshold is reused…".
- Acknowledge uncertainty explicitly: "we could not determine whether…".
- Short. Maintainers read many issues; details go into collapsible blocks.

### Signature (standard text)

> Found and verified as part of OSTRA (Open-Source elecTRonic-structure Audit), a physics-based audit of
> electronic-structure libraries. Analysis assisted by AI models (*list only the ones actually involved in
> this finding, e.g. Claude, Anthropic; Gemini, Google*); all results were checked and the report reviewed by
> a human before submission.

No institution in the signature (the project is open to external collaborators). The filer's GitHub account
identifies the person.

## 4. Human review and filing

1. The reviewer reads the whole draft and re-runs the reproducer.
2. Re-check that upstream has not changed (latest commit) and that no new related issue appeared.
3. Paste into GitHub and use **Preview**: tables render, `<details>` blocks open, code is highlighted,
   links work.
4. File it from a human account.

## 5. After filing

1. In `finding.md`, set `report` to `Reported: [owner/repo#N](URL) (YYYY-MM-DD)` and `status` to `reported`.
2. Regenerate `findings/INDEX.md` and commit (`OSTRA-NNN reported as owner/repo#N`); push.
3. Follow-up:
   - Answer maintainer questions with evidence (new runs go into `verification/` or `evidence/`).
   - No answer after ~3 weeks: one polite ping. No more after that.
   - If the maintainers say it is by design or will not fix: accept, record it as `reclassified` with their
     reasoning. That is a result too.
4. When fixed: record the fixing commit/PR and version, re-run the reproducer on it, set `status: resolved`.

## Report values used in `finding.md` / INDEX

`Pending` · `Pending (issue drafted in report/)` · `Reported: repo#N` · `Commented: repo#N` · `Fixed: repo#N` ·
`Not needed` (already fixed upstream) · `Not to be reported` (design, documented, or no impact).
