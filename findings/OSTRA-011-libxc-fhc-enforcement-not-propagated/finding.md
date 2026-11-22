# OSTRA-011: xc_func_set_fhc_enforcement does not propagate to auxiliary functionals of mixed meta-GGAs

| Field | Value |
|---|---|
| library | libxc |
| location | src/functionals.c:511-518 |
| versions | 7.1.2 and devel 3c7794c50 (function added in 7.0.0) |
| detected | 2026-11-21 |
| detected_by | Explorer 2; verified by explorer 1 |
| oracle | API contract (setter must act on composite functionals like its siblings) |
| class | API contract, new |
| impact_level | III (codes that disable FHC, e.g. PAW, with hybrid/mixed meta-GGAs) |
| status | Verified |
| report | Pending |

Every other setter (dens, zeta, sigma, tau thresholds) loops over `p->func_aux`; this one does not. Test with
tau < tau_W: disabling FHC changes eps for pure TPSS exchange but has no effect for TPSSh and TPSS0.
The explorer's count of 38 affected functionals was not re-verified.
