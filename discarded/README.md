# Discarded signals

Signals removed by the protocol. They are evidence that the filter works and part of the funnel reported in the paper.

| Library | Signal | Step that removed it | Reason |
|---|---|---|---|
| libxc | SCAN energy–potential deviations near alpha = 1 | 2 (artifact) | Finite-difference artifact: deviation moves with the step size |
| libxc | PBE-GX kink | 3 (design) | The kink is in the functional's published definition |
| libxc | 34 functionals differ spin-polarized vs unpolarized | 2 (artifact) | Disappear when density/sigma/tau thresholds are lowered: thresholds act per spin channel by design |
| libxc | `gga_c_epc19`: v_rho_up != v_rho_down at equal densities | 3 (design) | Electron–proton correlation uses the two channels for electrons and protons, not spins |
| libxc | Kinetic vW family and GG99 at s = 0 exactly | 2 (artifact) | Sigma threshold clamp; kinetic values ~1e-77 |
| libxc | Several s = 100, minimum-density spin discrepancies (1e-6 to 1e-5) | 5 (impact) | Precision noise at extreme inputs; several absent in devel |
| libcint | Overlap screening threshold | 3 (design) | Documented screening cut |
| libcint | `int2e_erfc` l = 5 with synthetic exponents (alpha = 0.5, R = 0.7, omega = 20) reaching 2.6e6 (explorer 2) | 5 (impact) | With a real basis (Ne2 cc-pV5Z, omega 0.11–20) erfc + erf = full holds to 8.6e-11. To revisit together with OSTRA-008 |
| dftd4 | Net-force flags for Si and NaCl | 2 (artifact) | Metric artifact; forces ~1e-18 by symmetry |
| PySCF PBC | Skewed unimodular Si cells, GDF/RSDF (explorer 2) | 2 (artifact) | Ewald, RHF(GDF), E_coul invariant; RKS drift comes from Becke grids and shrinks with grid level |
| PySCF | VV10 and implicit-solvent analytic gradients (explorer 2) | — (clean) | Analytic vs FD4 agree to 6e-10–1.5e-9 Eh/bohr |
| XCFun | M05X / M052XX vs libxc | 3 (convention) | Exact (1 - a_HF) scaling convention |
| libxc | wPBEh "0.10–0.20 kcal/mol vs CP2K" (explorer 2) | 3 (definition) | Compares two different rescaling definitions (Heyd 2004 vs Henderson 2009), not a defect |
