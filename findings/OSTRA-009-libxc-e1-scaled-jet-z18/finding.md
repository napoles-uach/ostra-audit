# OSTRA-009: xc_E1_scaled_jet switches to the asymptotic series at z = 18, where it is still inaccurate

| Field | Value |
|---|---|
| library | libxc |
| location | src/expint_e1.c (xc_E1_scaled_jet, z > 18.0) |
| versions | 7.1.2 and devel 3c7794c50 (not in 7.0.0) |
| detected | 2026-11-21 |
| detected_by | Explorer 2 (threshold inventory); verified by explorer 1 |
| oracle | Boundary test of branch switch; mpmath reference; energy-potential consistency |
| libraries | libxc |
| pattern | hard-branch-switch |
| oracle_class | analytic-reference, branch-boundary |
| class | New numerical bug |
| impact_level | II (to be measured in molecules) |
| status | Verified; attribution of the LC-wPBE inconsistency pending |
| report | Pending |

## Verified
- Just above z = 18 the optimally truncated asymptotic series gives relative errors 8.4e-8 (value), 1.6e-6 (1st
  derivative), 1.5e-5 (2nd), 9.5e-5 (3rd) vs mpmath; machine precision needs z >~ 35.
- LC-wPBE energy-potential inconsistency 3.2e-5 (relative, best FD step) in both 7.1.2 and devel; HSE06 4.5e-7;
  wPBEh with omega = 0 and PBE are clean.

## Corrections to the explorer's report
- Not a devel regression: the 7.1.2 release already contains the function with the z = 18 switch (introduced in
  7.1.1 to fix an overflow, ChangeLog #629 / !811); the explorer compared against 7.0.0.

## Pending
Build a copy with the switch moved and check that the LC-wPBE inconsistency disappears; propose a cut that does not
degrade the recurrence branch (which loses digits as z grows for high derivative orders).
