# OSTRA-019: simple-dftd3 gCP: effective atomic numbers beyond the 35/36-entry parameter tables

| Field | Value |
|---|---|
| library | simple-dftd3 |
| location | src/dftd3/gcp/param.f90 (effective_atomic_number; emiss_hf_* tables) |
| versions | v1.7.0 (84cb7cd, 2026-10-28) |
| detected | 2026-11-22 |
| detected_by | Explorer 2 (Gemini, Antigravity), batch 3; verified by explorer 1 (Claude) |
| oracle | Code reading; bounds-checked build (explorer) |
| libraries | simple-dftd3 |
| pattern | fixed-size-table |
| oracle_class | sanitizer-code-reading |
| class | Memory safety / fixed-size table, new |
| impact_level | III (heavy elements: Kr/Xe/Rn with MINIS, Fr/Ra/Ac, Z > 94) |
| status | Verified (code) |
| report | Pending |

effective_atomic_number maps 72:89 to Z - 50 (Fr/Ra/Ac -> 37-39) and leaves Z > 94 unchanged, while the
emission/basis tables have 36 entries (emiss_hf_minis has 35, so Kr/Xe/Rn -> 36 is already out of range). No element
validation in the parameter setup. Explorer's bounds-checked build reports E = +Inf / NaN gradients for Fr/Ra/Ac-O
(r2SCAN-3c) and a silent read of hydrogen's value for Kr-O (HF/MINIS); verifier confirmed the indices by reading the
code and counting table lengths, not by running the Fortran build.

Explorer's original record: `agents/explorer2/exploracion3-snapshot-2026-11-22/01-simple-dftd3-gcp-zeff-oob/finding.md` (copied below as `explorer_finding.md`).
