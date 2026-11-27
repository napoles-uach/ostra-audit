# OSTRA-027: XCFun TW (von Weizsäcker kinetic energy) squares the squared gradient: T_W of the H atom is 0.0075 Eh instead of 0.5

| Field | Value |
|---|---|
| library | XCFun |
| location | src/functionals/tw.cpp:21 (`1./8.*pow(d.gaa + d.gbb, 2.0)/d.n`) |
| versions | 2.1.1 (a89b783, master, unchanged since 2015) and the XCFun bundled with PySCF 2.15.1 |
| detected | 2026-11-26 |
| detected_by | Explorer 1 (Claude), guided search (pattern: formula transcription slips); blind verification by explorer 2 (Gemini), 2026-11-26 |
| oracle | Independent implementation (libxc GGA_K_VW) and exact limit (T_W = T_s for one orbital) |
| libraries | XCFun |
| pattern | transcription-slip |
| oracle_class | independent-implementation, exact-identity |
| class | New bug (formula) |
| impact_level | III (wrong by orders of magnitude when used; kinetic functional, used in orbital-free DFT and embedding) |
| status | Reported; Verified (blind, explorer 2), attributed, fix tested; ready to report |
| report | Reported: [xcfun#170](https://github.com/dftlibs/xcfun/issues/170) (2026-11-27, consolidated issue) |

## Observation
`d.gaa` and `d.gbb` are already squared gradients |∇ρ_a|², |∇ρ_b|², so the code computes (|∇ρ_a|²+|∇ρ_b|²)²/(8n)
instead of |∇n|²/(8n), and drops the cross term `gab`. H atom (aug-cc-pV5Z, UHF density, grid level 6):
<T> of the orbital 0.49998 Eh; libxc GGA_K_VW 0.49998 Eh; XCFun TW 0.0075 Eh. In the guided XCFun-vs-libxc scan the
integrated difference is ~1e10 mEh on molecular densities (dimensionally wrong formula).

## Documentation and history
Added in 17334f0 (2015-08-31). The placeholder test was removed in f81cd87 (2015-09-07, "The test was not
implemented / no test data given"), so the functional has never been tested.

## Cause and fix
One line: `return 1. / 8. * (d.gaa + 2.0 * d.gab + d.gbb) / d.n;` (or `d.gnn`). With the fix, H-atom T_W = libxc to
all printed digits; all other functionals bit-identical to the unpatched build (random-point regression).
`report/tw.patch`, `verification/` (scripts in `searches/2026-11-26-guided/xcfun_libxc/`).

## Impact
Only users of the TW kinetic functional (orbital-free DFT, frozen-density embedding, kinetic-energy studies) through
XCFun (PySCF exposes it as `TW`; also any XCFun host program). Silent when used.

## Blind verification and correction of the fix (explorer 2, 2026-11-26)
Confirmed at all three evidence levels (H atom 0.007505 Eh, equal to the analytic value of the wrong formula,
2/(27 pi^2); H2O 1.3e7 Eh vs 57.3 Eh). **Correction to explorer 1's patch:** `gnn/(8n)` is the von Weizsäcker
functional of the total density and does not satisfy the exact spin scaling of the kinetic energy,
T[n_a, n_b] = (T[2 n_a] + T[2 n_b])/2; explorer 1 had tested only the H atom (n_b = 0), where both agree. Open-shell
Li (UKS PBE/cc-pVTZ): libxc GGA_K_VW 7.178193 Eh, explorer 1's patch 7.132449 Eh (-45.7 mEh), reproduced by explorer 1.
The correct minimal fix is `0.125 * (d.gaa / d.a + d.gbb / d.b)` (now `report/tw.patch`): equal to libxc to 1e-15
in e, vrho, vsigma; 80/80 other functionals bit-identical; testall passes. XCFun already contains a correct,
spin-scaled von Weizsäcker functional, VWK (vonw.cpp), equal to libxc; TW could also be made an alias of it.
Material: `verification/gemini-2026-11-26/`.

