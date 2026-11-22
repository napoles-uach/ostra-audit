# Fragile numerical decisions

Recurring patterns behind the findings, independent of the library. Each pattern needs at least two independent
cases, ideally in different libraries, before it is presented as a pattern.

| Pattern | Cases | Cheap defence | Cases so far |
|---|---|---|---|
| Hard branch switches without checking continuity of value and derivatives (or switched on the wrong variable) | OSTRA-006 (s = 1, s = 15, w = 14), OSTRA-009 (z = 18), OSTRA-010 (PC07 cut), OSTRA-015 (alpha = 1), OSTRA-023 (B_k switched on abs(za-zb), not x) | Automatic test on both sides of every switch, for the value and every published derivative; switch on the variable that controls accuracy | 5 (3 libraries) |
| Thresholds or parameters valid for one parameter set, reused or left stale | OSTRA-010 (pc07_thr for b = 3 used with b = 0.258), OSTRA-014 (r4SCAN damping 0.232) | Derive thresholds from the tolerance; tie parameters to the cited publication in tests | 2 (2 libraries) |
| Comments or documented guarantees asserting an unverified numerical fact | OSTRA-010 ("~1.6e-17", really 4.6e-5), OSTRA-022 ("upper bound" that is a lower bound), OSTRA-023 ("12 terms seem accurate enough"), OSTRA-024 (comment says K_l = 0, code does not write it) | Turn the claim into an assertion or a test | 4 (3 libraries) |
| Fixed-size tables decoupled from the largest index the code can request | OSTRA-008 (TURNOVER_POINT[40] vs MXRYSROOTS = 32), OSTRA-019 (gCP zeff vs 36-entry tables), OSTRA-020 (z = 16 row; CART_CUM), OSTRA-021 (dftd4 max 118 vs multicharge 103) | `static_assert` linking both sizes; validate inputs against table limits at the API boundary; tests under sanitizers (ASan) | 4 (5 libraries) |
| Convergence parameters fixed for the typical case | OSTRA-005 (Ewald +-2 images), OSTRA-001 (tail check) | Adaptive convergence, or detect and warn | 2 (2 libraries) |
| Configuration or registration inconsistent across siblings | OSTRA-011 (FHC setter not propagated), OSTRA-012 (REVTPSSX missing XC_KINETIC), OSTRA-017 (placeholder ids passed to C) | Test that every setter/flag behaves identically for pure and composite objects; validate ids at the boundary | 3 (3 libraries) |
| Implicit API contracts | OSTRA-004 (ECP kernels leave output unwritten), OSTRA-024 (Bessel z <= 0), OSTRA-025 (gCP gradient needs sigma), OSTRA-026 (intent(out) dqdr unwritten) | Document postconditions; initialize outputs; treat each optional output independently | 4 (4 libraries) |
| Formula transcription slips that pass standard limits | OSTRA-013 (z^2 for z^3; passes z = 0 and z = 1), OSTRA-012 (missing square), OSTRA-016 (spin scaling) | Test with intermediate values of every variable, not only known limits; differential test against an independent implementation | 3 (1 library) |
| Regression tests against the code's own previous output | Common to all: no existing test caught these cases | Independent oracles: invariances, exact identities, analytic references | all |
