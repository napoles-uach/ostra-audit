# Impact study 1: XCFun meta-GGA defects on atomization and reaction energies (2026-11-24)

Findings: OSTRA-014 (R4SCANX damping parameters 0.232/0.232 instead of the published 0.802/0.178) and OSTRA-013
(REVTPSSC uses z^2 instead of z^3). Context: OSTRA-012 (REVTPSSX registration and formula).

## Design
- Same system, same basis, same grid; only the defect changes: released XCFun 2.1.1 (bundled with PySCF 2.15.1)
  vs the same XCFun source with the minimal fix (`findings/xcfun-confirmation-2026-11-21/` patches), loaded with
  LD_PRELOAD. A third build with libxc conventions ("conv") and libxc 7.1.2 are controls.
- 26 small molecules and 10 atoms (G2-like set; fixed, approximately experimental geometries: the comparison is
  between implementations at identical geometry, so geometry accuracy does not enter). def2-TZVP, grid level 5,
  conv_tol 1e-10, RKS/UKS. 10 reaction energies built from the same species.
- r4SCAN (r4scanx + r4scanc): fully self-consistent with each library.
- revTPSS correlation: released REVTPSSX cannot run (see below), so OSTRA-013 is isolated as the first-order
  change of E_c[REVTPSSC] on the self-consistent revTPSS density (libxc). Validated on r4SCAN: first-order vs
  self-consistent differences agree within 1-10 % (results/nsc.jsonl vs sc.jsonl).

## Results (released - fixed, kcal/mol)
| | MAD AE | max AE | MAD reaction | max reaction |
|---|---|---|---|---|
| r4SCAN (OSTRA-014) | 0.63 | 1.99 (C2H6) | 0.58 | 2.03 (2 H2O2 -> 2 H2O + O2) |
| revTPSS correlation (OSTRA-013) | 0.54 | 1.35 (C2H6) | 0.59 | 1.32 (N2 + 3 H2 -> 2 NH3) |

Errors grow with the number of bonds (C2H6, CH3OH, H2CO) and are of the order of chemical accuracy (1 kcal/mol),
i.e. the size of differences used to rank functionals in benchmarks. No warning is produced.

## Controls
- Released XCFun from our build == PySCF bundled library (identical values).
- Fixed XCFun reproduces libxc: revTPSS total energies to 1e-11 Eh for all molecules; r4SCAN molecular energies to
  < 0.005 kcal/mol. Residual 0.03 / 0.07 kcal/mol only for the open-shell O and F atoms (different SCF solutions
  for p-shell atoms), which therefore also appears in AEs of O/F molecules; reaction energies agree to 0.001.
- revTPSS with the released XCFun (REVTPSSX without XC_KINETIC, OSTRA-012): the SCF does not finish
  (LinAlgError: singular matrix or > 15 min) for NH3, H2O, HF, C2H2, CH4, C2H4, C2H6. Loud, not silent.

## Side observation (new candidate)
XCFun TPSSX and REVTPSSX return NaN for the H atom (fully spin-polarized one-electron density) in both released and
fixed builds; libxc returns finite values. Correlation is fine (~1e-12). `scripts/hatom.py`.

Figure: results/fig_xcfun_impact.png (.svg).
