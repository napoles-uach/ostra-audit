# OSTRA-020: PySCF ECP: table index 400 at z = 16 and crashes of ECPscalar_ipnuc at high angular momentum

| Field | Value |
|---|---|
| library | PySCF (libcgto) |
| location | pyscf/lib/gto/nr_ecp.c (ECPsph_ine_opt; omega_nuc[CART_CUM]), nr_ecp.h |
| versions | PySCF 2.14.0, 2.15.1 and master 484fd7c |
| detected | 2026-11-22 |
| detected_by | Explorer 2 (Gemini, Antigravity), batch 3; verified by explorer 1 (Claude) |
| oracle | ASan (explorer); crash / exit code (verifier) |
| class | Memory safety, new |
| impact_level | III (A: only z exactly 16; B: ECP derivative integrals with i functions crash) |
| status | A verified; B crash verified, cause partly attributed |
| report | none yet |

A. ECPsph_ine_opt uses `z > 16` to leave the [0,16) tables; floor(16.0/(16/400)) = 400, one row past the
400-row tables (only for z exactly 16). B. `double omega_nuc[CART_CUM]` with CART_CUM = 456 (up to l = 12). Verifier,
PySCF 2.15.1, `mol.intor("ECPscalar_ipnuc")`: l_basis = 5 with lc = 5, 6 works; l = 5, lc = 7 (l+1+lc = 13) segfaults,
consistent with CART_CUM; but l_basis = 6 segfaults for any lc (0, 2, 3, 4), so another limit is involved. ECPscalar
itself works for l = 6. Next: isolate the l = 6 crash (assigned to explorer 2).

Explorer's original record: `agents/explorer2/exploracion3-snapshot-2026-11-22/02-pyscf-ecp-fixed-tables/finding.md` (copied below as `explorer_finding.md`).
