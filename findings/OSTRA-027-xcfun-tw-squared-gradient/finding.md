# OSTRA-027: XCFun TW (von Weizsäcker kinetic energy) squares the squared gradient: T_W of the H atom is 0.0075 Eh instead of 0.5

| Field | Value |
|---|---|
| library | XCFun |
| location | src/functionals/tw.cpp:21 (`1./8.*pow(d.gaa + d.gbb, 2.0)/d.n`) |
| versions | 2.1.1 (a89b783, master, unchanged since 2015) and the XCFun bundled with PySCF 2.15.1 |
| detected | 2026-11-26 |
| detected_by | Explorer 1 (Claude), guided search (pattern: formula transcription slips) |
| oracle | Independent implementation (libxc GGA_K_VW) and exact limit (T_W = T_s for one orbital) |
| libraries | XCFun |
| pattern | transcription-slip |
| oracle_class | independent-implementation, exact-identity |
| class | New bug (formula) |
| impact_level | III (wrong by orders of magnitude when used; kinetic functional, used in orbital-free DFT and embedding) |
| status | Verified (run), attributed, fix tested; needs blind verification by a second agent |
| report | Pending |

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
