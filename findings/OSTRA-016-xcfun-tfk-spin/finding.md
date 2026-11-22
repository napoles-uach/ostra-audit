# OSTRA-016: TFK (Thomas-Fermi kinetic) ignores spin polarization

| Field | Value |
|---|---|
| library | XCFun |
| location | src/functionals/tfk.cpp |
| versions | 2.1.1 (PySCF 2.15.1 bundled; also macOS PySCF 2.14.0); same code on master a89b783 |
| detected | 2026-11-21 |
| detected_by | Explorer 2 (XCFun vs libxc differential audit); verified by explorer 1 |
| oracle | Exact spin scaling T[ra, rb] = (T[2ra] + T[2rb]) / 2 |
| class | New bug (formula) |
| impact_level | III (spin-polarized kinetic-energy use, e.g. embedding) |
| status | Verified |
| report | none yet |

Code is `CF * pow(d.n, 5/3)`. Error -8.21 % at zeta = 0.4 (explorer: -26.8 % at 0.8, -37 % at 1).
