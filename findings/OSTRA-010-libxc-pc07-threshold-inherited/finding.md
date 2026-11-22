# OSTRA-010: pc07_thr = 1/40 inherited by PC07-opt (SCAN-L, r2SCAN-L): a 4.6e-5 indicator value is truncated to zero

| Field | Value |
|---|---|
| library | libxc |
| location | maple/mgga_exc/mgga_k_pc07.mpl:24-27 |
| versions | devel 3c7794c50 |
| detected | 2026-11-21 |
| detected_by | Explorer 2; verified by explorer 1 |
| oracle | Boundary test; analytic evaluation |
| class | New numerical bug |
| impact_level | II (zk jump ~1.8e-6 relative per point; molecular impact not measured) |
| status | Verified analytically |
| report | Pending |

The comment justifies the cut with "pc07_fab0(a/40) ~ 1.6e-17"; that number omits the power b. With the original
PC07 parameters (a = 0.5389, b = 3) the value is 4.2e-51, but SCAN-L and r2SCAN-L (exchange and correlation) use
PC07-opt (a = 1.784720, b = 0.258304), where f(a/40) = 4.59e-5 and is clamped to zero. The upper cut is fine
(1 - f(39a/40) ~ 1e-17).
