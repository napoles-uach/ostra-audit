# OSTRA-013: REVTPSSC uses z^2 instead of z^3 in eps_c = eps_PKZB (1 + d eps_PKZB z^3)

| Field | Value |
|---|---|
| library | XCFun |
| location | src/functionals/revtpssc_eps.hpp:107-109 |
| versions | 2.1.1 (PySCF 2.15.1 bundled; also macOS PySCF 2.14.0); same code on master a89b783 |
| detected | 2026-11-21 |
| detected_by | Explorer 2 (XCFun vs libxc differential audit); verified by explorer 1 |
| oracle | Code comparison with TPSSC; differential vs libxc |
| libraries | XCFun |
| pattern | transcription-slip |
| oracle_class | independent-implementation |
| class | New bug (formula), silent |
| impact_level | I-II (explorer: SCF shifts 1.5-2.1 kcal/mol; not re-measured) |
| status | Confirmed (level 3): minimal fix matches the publication and libxc to 4e-11 |
| report | Pending |

`revtpssc_eps.hpp` uses `pow2(...)` where `tpssc_eps.hpp` (same formula, same d = 2.8) uses `pow3(...)`.
Max relative difference vs libxc 2 % on 200 random points (control TPSSC vs libxc: 3e-11). Silent: passes the
uniform gas (z = 0) and one-electron densities (z = 1).

## Confirmation phase (2026-11-21) — level 3: defect scientifically confirmed
- **Publication check.** revTPSS preprint (arXiv:0906.5029): "Aside from our use of Eq. (3), we keep the form of TPSS
  correlation unchanged", i.e. eps_c = eps_c^revPKZB [1 + d eps_c^revPKZB z^3], d = 2.8; new beta(rs) and C(zeta, 0)
  (XCFun's coefficients match the paper). Erratum PRL 106, 179902 does not touch the formulas.
- **Minimal fix:** z^3 in the final expression only (the PKZB part keeps z^2, correctly).
- **Result** vs libxc 7.1.2: eps 2.7e-2 -> 1.3e-5 (derivatives 69 -> 1.4e-3). The remainder is a convention: libxc uses
  the exact PBE beta (0.06672455060314922, commented in gga_c_regtpss.mpl) instead of the paper's 0.066725. With it:
  eps 4.0e-11, derivatives 1.6e-8, identical to the TPSSC control (3.9e-11).
- **Status:** confirmed defect (formula). The beta rounding is a convention, not a defect.

## Impact on observables (2026-11-24)
`studies/impact-xcfun-2026-11-24/`: first-order change of E_c on self-consistent revTPSS densities, 26 molecules,
def2-TZVP. Atomization energies: MAD 0.54, max 1.35 kcal/mol (C2H6); reaction energies: MAD 0.59, max 1.32 kcal/mol
(N2 + 3 H2 -> 2 NH3). Silent. Impact level I confirmed (order of chemical accuracy).
