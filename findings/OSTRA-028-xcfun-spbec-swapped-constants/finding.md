# OSTRA-028: XCFun SPBEC (sPBE correlation for SSB) has the PBE beta and gamma constants swapped: E_c 63 % too negative

| Field | Value |
|---|---|
| library | XCFun |
| location | src/functionals/spbec.cpp:22-27 (gamm = 0.066725, beta = 0.031091; G uses beta_gamma*phi3) |
| versions | 2.1.1 (a89b783, master; constants unchanged since the import from Dalton, 2010-2011) and PySCF 2.15.1 bundle |
| detected | 2026-11-26 |
| detected_by | Explorer 1 (Claude), guided search (pattern: formula transcription slips); blind verification by explorer 2 (Gemini), 2026-11-26 |
| oracle | Independent implementation (libxc GGA_C_SPBE) and re-implementation from the published formula |
| libraries | XCFun |
| pattern | transcription-slip |
| oracle_class | independent-implementation, analytic-reference |
| class | New bug (parameters and formula) |
| impact_level | I-III (I for anyone computing SSB/sPBE correlation with XCFun; usage to be assessed) |
| status | Reported; Verified (blind, explorer 2), attributed, fix tested, checked against the paper; ready to report |
| report | Reported: [xcfun#170](https://github.com/dftlibs/xcfun/issues/170) (2026-11-27, consolidated issue) |

## Observation
H2O (def2-SVP, PBE density): sPBE correlation from the paper formula with PBE constants (beta = 0.066725,
gamma = (1 - ln 2)/pi² = 0.031091) gives -0.368141 Eh with VWN5 LDA (as XCFun) and -0.365746 Eh with PW92 LDA;
libxc GGA_C_SPBE -0.365745 Eh; XCFun SPBEC -0.597502 Eh (63 % more correlation; PBE correlation is -0.331 Eh).

## Cause and fix
The variable named `gamm` holds beta and `beta` holds gamma, so beta/gamma is inverted (0.466 instead of 2.146) and
the prefactor of H uses beta instead of gamma; G also divides by beta_gamma*phi3 instead of gamma*phi3. Swapping the
two values and using `gamm * phi3` in G gives -0.368140 Eh (= the paper formula with VWN5 to 1e-6 Eh). Other
functionals bit-identical. `report/spbec.patch`.

## Check against the paper (2026-11-26)
Swart, Solà, Bickelhaupt, J. Chem. Phys. 131, 094103 (2009), Eq. (9): H_sPBE = gamma phi^3 ln[1 + (beta/gamma) t^2
/ (1 + G t^2)], G = (beta/gamma) / (exp(-eps_c^unif / (gamma phi^3)) - 1), and the text keeps beta and gamma at the
PBE values "~0.066 725, ~0.031 091". This is exactly the form re-implemented in the verification and the one the
patched XCFun computes; the released code has beta and gamma swapped and uses beta/gamma in place of gamma inside G.
(The PDF is not stored in the repository.)

## Blind verification (explorer 2, 2026-11-26)
Confirmed at all three levels: H2O (cc-pVTZ) -0.595478 vs libxc -0.362040 Eh (-64 %); vsigma off by a factor ~3 at a
polarized point. Additional physical evidence: the released code violates the exact rapidly-varying limit
eps_c -> 0 for t -> infinity (n = 0.1, |grad n| = 1e4: -4.58e-2 vs -3.95e-11), because of the beta_gamma*phi3 inside G;
explorer 1's patch restores it (-3.97e-11, checked). Class corrected to parameters and formula. LDA choice: explorer 1's
patch keeps VWN5 (XCFun's choice, consistent with the ADF context of the paper), leaving 2.4 mEh vs libxc, which uses
PW92; explorer 2's alternative patch switches to PW92 and matches libxc to 1e-15. Both are given in the report; the
choice is the maintainers'. Material: `verification/gemini-2026-11-26/`.

## Pending
Assess who uses XCFun's SPBEC (SSB-D in XCFun host programs). No test
reference data in the source.

