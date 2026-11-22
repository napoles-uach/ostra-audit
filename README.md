# OSTRA — Open-Source elecTRonic-structure Audit

OSTRA is a physics-based, AI-assisted audit of the shared numerical libraries that electronic-structure codes
depend on (exchange-correlation functionals, Gaussian integrals, dispersion and charge models). Its target is
**silent errors**: the code does not crash, but returns plausible, wrong numbers. Such errors have no natural
oracle, so ordinary testing and security tooling do not see them. OSTRA uses physics as the oracle: exact
identities, invariances, energy–potential consistency, analytic references and independent implementations.

Project at the Universidad Autónoma de Chihuahua (UACH). AI models (Claude, Anthropic, and a second agent) build
and run the verification harnesses; every result is checked, and every report reviewed, by a human before submission.

> **Status: private working repository.** Findings not yet reported to maintainers are kept here until they have
> passed independent verification. Each upstream report links back to its folder in this repository.

## Layout

| Path | Contents |
|---|---|
| `PROTOCOL.md` | The six-step audit protocol, finding classes, impact levels, status workflow |
| `findings/INDEX.md` | Table of all findings with class, impact level and status |
| `findings/OSTRA-NNN-*/` | One folder per finding: `finding.md` (record), `environment.txt` (exact versions), `repro/` (minimal reproducers), `evidence/` (raw outputs, figures, data), `verification/` (independent checks), `report/` (draft or submitted issue text) |
| `candidates/` | Signals not yet verified, or not reproduced |
| `discarded/` | Signals removed by the protocol, and why (also evidence that the filter works) |
| `catalog/` | Fragile numerical decisions: recurring patterns across libraries, with linked cases |
| `probes/` | General harnesses and raw scan results |
| `agents/` | Log of which agent proposed what and how verification changed it |

## Rules

1. Nothing is a finding until it passes the protocol; nothing is reported without human approval.
2. Record the state **at detection time**: library version and commit, environment, exact inputs, raw outputs.
   Later corrections are appended, never overwritten (see the "Corrections" sections).
3. Before any report, check that it is not already known (issue trackers, changelogs, the original paper).
4. Prefer collecting and verifying over filing issues one by one; issues are filed once a finding is mature.

## Summary (2026-11-22)

26 findings in nine libraries (libxc, libcint, PySCF, multicharge, dftd4, XCFun, simple-dftd3, libecpint, and
PySCF's XCFun interface). Through OSTRA-017 (2026-11-21):
12 new defects (one may be outdated preprint parameters, OSTRA-014), 2 known cases rediscovered and quantified,
1 already fixed upstream before our report, 1 documented limitation, 1 case with no impact. One fix merged upstream (libxc#646), two reports filed
(pyscf#3669, multicharge#72 comment). OSTRA-018 to OSTRA-026 (2026-11-22, explorer 2 batch 3, independently verified) add
seven new defects and the confirmed PK09 anomaly. See `findings/INDEX.md`.
