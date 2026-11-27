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

## Batch 3 — explorer 2 (Gemini), catalog second cases (2026-11-22)
First batch with explicit verification feedback and read access to this repository. Eight candidates in five
libraries plus the answer on lda_c_pk09. Independent verification (explorer 1, own tests and code reading at the
same commits: simple-dftd3 84cb7cd, dftd4 4a5bbf0, multicharge 4bb6c96, libecpint d7966f1, PySCF master 484fd7c):

| # | Candidate | Verdict | Notes |
|---|---|---|---|
| pk09 | NaN band of lda_c_pk09 near r_s = 22.739 | Confirmed (polarized only) | zeta = 0.4: NaN for |dr_s| < ~5e-5, eps_c = -1287 Ha nearby (7.1.2). Model pole vs implementation pending the paper |
| 01 | simple-dftd3 gCP: zeff beyond 35/36-entry tables | Confirmed (code) | zeff = Z - 50 for Fr/Ra/Ac (37-39), zeff = Z for Z > 94; emiss_hf_minis has 35 entries (Kr/Xe/Rn -> 36); no element validation |
| 02 | PySCF ECP: z > 16 vs table [0,16); CART_CUM | A confirmed, B crash confirmed but only partly attributed | floor(16.0/0.04) = 400 (exactly z = 16 only). ECPscalar_ipnuc segfaults (exit 139) for basis l = 6 with any lc, and for l = 5, lc = 7 (li+1+lc = 13 > 12 fits CART_CUM); the l = 6 / any-lc crash needs another cause |
| 03 | dftd4 accepts Z = 112-118, multicharge max_elem = 103 | Confirmed (released dftd4 4.3.1, Python) | Cn-O and Og-O: Cn/Og get q = -0.87, O +0.87, no error |
| 04 | libecpint upper_bound is a lower bound on the rising side | Confirmed (code); impact downgraded | No caller inside libecpint: affects only external users of the public header |
| 05 | simple-dftd3 gCP: 12-term series vs closed form switched on abs(za-zb) | Partly verified | Series error 3.5e-4 to 1.8e-3 at x = 5; closed form 4e-3 to 0.92 at x = 0.05-0.1 for k = 7-8; impact on gCP energies not re-measured |
| 06 | libecpint calculate(z <= 0) writes only K_0; scalar overload returns 1 | Confirmed (code); impact downgraded | Internal callers guard x > 0 or are safe by zero-initialization; affects external users |
| 07 | simple-dftd3 gCP gradient requires gradient .and. sigma | Confirmed (code) | C API with NULL sigma passes an unallocated array (absent), gradient stays zero; disp.f90 uses .or. |
| 08 | multicharge get_charges: intent(out) dqdr unwritten without dqdL | Confirmed (code) | grad/cpq require both; dqdr left undefined |

Corrections to the explorer's summary: "8 confirmed defects" holds for the code defects, but two (04, 06) have no
impact inside their own library and 02-B is not fully attributed. Catalog gains: fixed-size tables (01, 02-A, 03),
unverified numerical comments/contracts in comments (04, 05, 06), implicit API contracts (06, 07, 08).

## Batch 4 — explorer 2, isolation of OSTRA-020 part B (2026-11-22)
Single assigned task. Verifier outcome: cause confirmed by reading the code (local-channel `omega_nuc` in
`type1_rad_ang`), practical impact confirmed with real built-in basis/ECP pairs (AgH cc-pV5Z-PP + def2-SVP: ECP
gradient integrals segfault). Two additional ASan-only claims (cache-size underestimate, ECP_LMAX) not reproduced as
crashes in an ordinary build; pending an ASan build. Verifier error recorded: a first attempt called second-derivative
integrals without `comp=9` and crashed for that reason; corrected before drawing conclusions.

## 2026-11-23 — Blind verification of OSTRA-008 by explorer 2 (Gemini 4 Neon, Antigravity)
First finding by explorer 1 verified blind by another agent (location and claim only; patch given after phase 1).
Result: anomaly, cause and defect confirmed on macOS arm64 with four binaries, ASan and mpmath / full - LR oracles.
Corrections to explorer 1: the m = 71 bound is loose (reachable max 47 in SR, 62-63 overall) and a guard is still
needed for direct calls with m >= 72. New side observations moved to candidates: double-precision loss above t0
(partly re-checked by explorer 1) and an lrys_wheeler singularity at R = 1.0 on arm64 without quadmath.
Two 503 capacity errors from the model API during the run; no effect on results.

## 2026-11-26 — Guided search round 1 (catalog-driven)
Explorer 1: XCFun vs libxc over all 104 comparable names (new OSTRA-027 TW, OSTRA-028 SPBEC; 4 known re-detected) and
equivalent-cell invariance (D3, PySCF Ewald, tblite clean; dftd4 gradient inconsistency extends OSTRA-005).
Explorer 2: implicit-API-contract sweep of tblite/dftd4/simple-dftd3/mctc-lib, 56 routines in ~69 min, 5 root causes
with tested patches. Verification by explorer 1: dftd4 C API sigma-only crash reproduced with the released 4.3.1
(gradient-only works); mctc-lib dcndr-only reproduced with an independent Fortran program (sentinel kept; rescaled
when cut > 0); tblite, simple-dftd3 model and dftd4 model conditions confirmed by reading main. Observation: the same
`present(a) .and. present(b)` convention recurs in mctc-lib, multicharge, dftd4, simple-dftd3 and tblite, i.e. a
pattern propagated across one software ecosystem.

## 2026-11-26 — Blind verification of OSTRA-027 and 028 by explorer 2
Both confirmed at all three evidence levels. Explorer 2 corrected explorer 1's TW patch: the total-density form
gnn/(8n) breaks the exact spin scaling (open-shell Li -45.7 mEh vs libxc), which explorer 1 missed by testing only the
H atom; reproduced by explorer 1 and the patch replaced. Explorer 2 also found that XCFun already has a correct
von Weizsäcker functional (VWK), and added the rapidly-varying limit as a physical oracle for SPBEC. Second time a
blind verification corrects the proposer (after OSTRA-008's m = 71 bound).

