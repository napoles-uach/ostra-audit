# Findings index

Impact levels: I = changes publishable observables; II = breaks numerical consistency (energy/force/potential); III = niche but severe or silent when hit; IV = no impact. See `../PROTOCOL.md`.

| ID | Library | Finding | Class | Impact | Status |
|---|---|---|---|---|---|
| [OSTRA-001](OSTRA-001-libxc-xc_integrate-tail-check/finding.md) | libxc | Tanh-sinh quadrature tail check in xc_integrate truncates integrals with mass away from the midpoint | New numerical bug | III | Resolved upstream |
| [OSTRA-002](OSTRA-002-libxc-hse-vsigma-noise/finding.md) | libxc | Noisy v_sigma of HSE-family exchange at small reduced gradient | Numerical precision | II | Already fixed in devel by the developers |
| [OSTRA-003](OSTRA-003-libcint-highl-eri-precision/finding.md) | libcint | Four-centre ERIs lose precision at high angular momentum (Rys roots) | Documented limitation | IV | Documented; not reported |
| [OSTRA-004](OSTRA-004-pyscf-ecp-iprinv-uninitialized/finding.md) | PySCF (libcgto) | ECP-centre derivative kernels return uninitialized memory when the rinv atom has no ECP | API contract, new | III | Reported |
| [OSTRA-005](OSTRA-005-multicharge-ewald-fixed-images/finding.md) | multicharge / dftd4 | EEQ Ewald sums with a fixed +-2 image range: charges depend on the equivalent unit cell | Convergence, known | I | Known issue; data contributed |
| [OSTRA-006](OSTRA-006-libxc-wpbeh-s-scaling-discontinuity/finding.md) | libxc | HSE03/06/12, LC-wPBE and wPBEh energy is discontinuous at s = 1 (s_scaling_2) | Known design (Henderson et al. 2009, eq. 8 vs eq. 9) | II | Classified; not to be reported (deliberate fidelity to HSE06) |
| [OSTRA-007](OSTRA-007-libxc-tpss-c-alpha0-derivative/finding.md) | libxc | Derivative of TPSS/TM-family correlation at tau = tau_W (z = 1) depends on rounding | No impact | IV | Documented; not reported |
| [OSTRA-008](OSTRA-008-libcint-turnover-point-oob/finding.md) | libcint | Out-of-bounds read of TURNOVER_POINT[40] for Boys order m >= 40 corrupts SR integrals at high l | Memory safety, new | III | Verified and attributed; issue drafted |
| [OSTRA-009](OSTRA-009-libxc-e1-scaled-jet-z18/finding.md) | libxc | xc_E1_scaled_jet switches to the asymptotic series at z = 18, where it is still inaccurate | New numerical bug | II | Verified; attribution of the LC-wPBE inconsistency pending |
| [OSTRA-010](OSTRA-010-libxc-pc07-threshold-inherited/finding.md) | libxc | pc07_thr = 1/40 inherited by PC07-opt (SCAN-L, r2SCAN-L): a 4.6e-5 indicator value is truncated to zero | New numerical bug | II | Verified analytically |
| [OSTRA-011](OSTRA-011-libxc-fhc-enforcement-not-propagated/finding.md) | libxc | xc_func_set_fhc_enforcement does not propagate to auxiliary functionals of mixed meta-GGAs | API contract, new | III | Verified |
| [OSTRA-012](OSTRA-012-xcfun-revtpssx/finding.md) | XCFun | REVTPSSX registered without XC_KINETIC and with two formula errors | New bug (registration + formula) | III | Verified by code reading; NaN reproduced |
| [OSTRA-013](OSTRA-013-xcfun-revtpssc-z2/finding.md) | XCFun | REVTPSSC uses z^2 instead of z^3 in eps_c = eps_PKZB (1 + d eps_PKZB z^3) | New bug (formula), silent | I-II | Verified |
| [OSTRA-014](OSTRA-014-xcfun-r4scanx-damping-params/finding.md) | XCFun | R4SCANX uses d_p4 = d_a4 = 0.232 instead of 0.802 / 0.178 | New bug (parameters) or outdated preprint values | I-II | Verified vs libxc; published values to be confirmed in Furness et al., JCP 156, 034109 (2022) |
| [OSTRA-015](OSTRA-015-xcfun-scan-alpha1-s0/finding.md) | XCFun | SCANX/SCANC return +-inf at alpha = 1 (uniform gas) and SCANX v_sigma is halved at s = 0 | Unguarded edge case | III-IV | Verified |
| [OSTRA-016](OSTRA-016-xcfun-tfk-spin/finding.md) | XCFun | TFK (Thomas-Fermi kinetic) ignores spin polarization | New bug (formula) | III | Verified |
| [OSTRA-017](OSTRA-017-pyscf-xcfun-vv10-segfault/finding.md) | PySCF | Segmentation fault for VV10/rVV10 functionals with the XCFun backend | New bug (crash) | III | Verified |
