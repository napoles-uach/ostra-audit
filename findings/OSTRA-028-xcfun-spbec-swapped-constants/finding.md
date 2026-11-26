# OSTRA-028: XCFun SPBEC (sPBE correlation for SSB) has the PBE beta and gamma constants swapped: E_c 63 % too negative

| Field | Value |
|---|---|
| library | XCFun |
| location | src/functionals/spbec.cpp:22-27 (gamm = 0.066725, beta = 0.031091; G uses beta_gamma*phi3) |
| versions | 2.1.1 (a89b783, master; constants unchanged since the import from Dalton, 2010-2011) and PySCF 2.15.1 bundle |
| detected | 2026-11-26 |
| detected_by | Explorer 1 (Claude), guided search (pattern: formula transcription slips) |
| oracle | Independent implementation (libxc GGA_C_SPBE) and re-implementation from the published formula |
| libraries | XCFun |
| pattern | transcription-slip |
| oracle_class | independent-implementation, analytic-reference |
| class | New bug (parameters) |
| impact_level | I-III (I for anyone computing SSB/sPBE correlation with XCFun; usage to be assessed) |
| status | Verified (run), attributed, fix tested; needs blind verification and a check against the paper |
| report | Pending |

## Observation
H2O (def2-SVP, PBE density): sPBE correlation from the paper formula with PBE constants (beta = 0.066725,
gamma = (1 - ln 2)/pi² = 0.031091) gives -0.368141 Eh with VWN5 LDA (as XCFun) and -0.365746 Eh with PW92 LDA;
libxc GGA_C_SPBE -0.365745 Eh; XCFun SPBEC -0.597502 Eh (63 % more correlation; PBE correlation is -0.331 Eh).

## Cause and fix
The variable named `gamm` holds beta and `beta` holds gamma, so beta/gamma is inverted (0.466 instead of 2.146) and
the prefactor of H uses beta instead of gamma; G also divides by beta_gamma*phi3 instead of gamma*phi3. Swapping the
two values and using `gamm * phi3` in G gives -0.368140 Eh (= the paper formula with VWN5 to 1e-6 Eh). Other
functionals bit-identical. `report/spbec.patch`.

## Pending
Check the constants and the form of H and A against Swart, Solà, Bickelhaupt, JCP 131, 094103 (2009). Assess who uses
XCFun's SPBEC (SSB-D in Dalton/DIRAC/ADF-like hosts?). No test reference data in the source.
