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
| status | Confirmed (level 3): exact spin scaling after minimal fix |
| report | none yet |

Code is `CF * pow(d.n, 5/3)`. Error -8.21 % at zeta = 0.4 (explorer: -26.8 % at 0.8, -37 % at 1).

## Confirmation phase (2026-11-21)
- **Polarization contract.** Sibling kinetic functionals in XCFun are spin-resolved (BTK: 0.5*(f(2 rho_a) + f(2 rho_b));
  PW91K evaluates per spin), so TFK is inconsistent with its own library. Its built-in self-test uses rho_a = 1,
  rho_b = 0.8 and the reference 7.6475577 equals C_F (1.8)^(5/3), the spin-unaware formula (spin-scaled value 7.699):
  the regression test froze the error.
- **Minimal fix:** C_F 2^(2/3) (rho_a^(5/3) + rho_b^(5/3)). Result: exact spin scaling at zeta = 0.4, 0.8, 1.0 (<=1e-16).
  The self-test reference values must be updated together with the fix.
- **Status:** confirmed defect.
