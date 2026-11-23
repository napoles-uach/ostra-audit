# Candidates (not yet verified, or not reproduced)

Signals reported by an explorer that have not passed independent verification. Each needs the explorer's exact
conditions before it can move to `findings/` or `discarded/`.

| Proposed by | Library | Signal | Verifier status (2026-11-21) |
|---|---|---|---|
| Explorer 2 | libxc | `lda_c_pk09`: NaN band around the pole r_s = 22.739412 (rho ~ 2e-5) | Promoted to **OSTRA-018** (2026-11-22): reproduces only with spin polarization (0 < abs(zeta) < 1) |
| Explorer 2 | libxc | TPSS / TPSSh: v_sigma_bb jump near 1 - abs(zeta) = 1e-12 | Not verified. Likely the zeta threshold (by design). Note: the explorer quoted "+873 Eh" for v_sigma, which has the wrong units |
| Explorer 2 | libxc | `mgga_c_lak`: `m_max(rs*a, eps)` algebraic error for alpha < 0 | Not verified. alpha < 0 means tau < tau_W, unphysical and normally removed by FHC enforcement; check whether reachable |
| Explorer 2 | libxc | LB94: 0.32 % jump at x_s = 300 | Not verified. Probably the asymptotic switch of the method (design) |
| Explorer 1 | libxc | MCML / VCML small flags in devel consistency scan; BEEF noise | Open since 2026-11-19 |
| Explorer 1 | dftd3 (simple-dftd3) | Triclinic FD deviations in periodic D3 | Unclassified (runs timed out); likely hard cutoffs |
| Explorer 1 | libxc | LTA / TH / TLDA family differ at tau = 0 with rho > 0 (7.1.2 `mgga_x_lta` returned 1e273) | Degenerate input; low priority |
| Explorer 2 | libxc | Inventory of 84 numerical decisions and 143 boundary tests | Material lives in the explorer's folder (`PhysAud/`); to be imported verbatim |
| Explorer 2 | libcint | Double-precision Boys F_m loses digits just above TURNOVER_POINT (pre-existing: 1e-8 at m = 39; 1e-4 at m = 63); matters where long double = double and no quadmath (macOS arm64) | Partly re-checked by explorer 1 (2026-11-23, values in OSTRA-008 finding.md); impact on integrals not measured |
| Explorer 2 | libcint | lrys_wheeler singularity at R = 1.0 on macOS arm64 without libquadmath (also with m = 39) | Not verified; seen during OSTRA-008 verification |
