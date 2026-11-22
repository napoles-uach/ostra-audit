# Agents and verification log

| Agent | Role | Environment |
|---|---|---|
| Explorer 1 | Main session (Claude Opus 5.5, Anthropic): exploration, triage, verification of explorer 2 | Cloud container + Linux VM with libxc devel |
| Explorer 2 | Gemini in Google Antigravity ("Neon"), working in `PhysAud/` (macOS), given the context prompt | macOS, PySCF 2.14.0 / libxc 7.0.0 / XCFun 2.1.1 |
| Reviewer | ChatGPT 6.2 (OpenAI): critical review of reports and of the publication plan; reads public code, does not execute | — |
| Human | J. M. Nápoles-Duarte (UACH): decisions, reporting, final review | — |

## Batch 1 — explorer 2, libxc threshold inventory (2026-11-21)
Nine claims. Verifier outcome: 3 confirmed (OSTRA-010, OSTRA-011, and OSTRA-009 with correction: not a devel
regression, already in 7.1.2), 1 confirmed but minor and partly a definition comparison (s = 15 / w = 14 in
OSTRA-006), 1 not reproduced (pk09 NaN), 3 low priority / likely design (candidates). Overclaims corrected: "done at
100 %" and "control group confirms the pattern" (selection bias: boundary tests were placed where errors were expected).

## Batch 2 — explorer 2, XCFun and PySCF (2026-11-21)
Seven claims. Verifier outcome: 6 confirmed (OSTRA-012 to OSTRA-017), 1 consistent with a known case (TDDFT
gradients, folded into OSTRA-006). Open question for the paper: whether explorer 2 received the batch-1 corrections
before batch 2 (affects interpretation of the higher confirmation rate).

## Review — ChatGPT 6.2 on batch 2 (2026-11-21)
Code-reading review (no execution). Agrees on REVTPSSX (registration and both formula errors), REVTPSSC (z^2) and
R4SCANX (parameters) as strong evidence; VV10 segfault very likely (adds: `xcfun_itrf.c` passes the name from
`xcfun_enumerate_parameters` to `xcfun_set` without a NULL check); TFK needs clarification of the polarization
contract; SCAN edge needs one-sided limits; TDDFT causality pending. Proposals adopted: three evidence levels
(anomaly reproduced / cause identified / defect scientifically confirmed = minimal fix matches the publication and
removes the anomaly without side effects); distinguish the local bound on F_x from the global Lieb-Oxford
inequality; judge impact on energy differences (reaction/atomization), not only total energies; pause exploration
for a confirmation phase.
