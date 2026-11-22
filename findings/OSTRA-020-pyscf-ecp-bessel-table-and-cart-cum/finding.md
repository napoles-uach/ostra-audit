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
| impact_level | II-III (A: only z exactly 16; B: nuclear gradients crash with built-in i-function bases and ECPs with a local channel) |
| status | A verified; B attributed (local-channel omega_nuc) and impact measured; extra ASan claims pending |
| report | Pending |

A. ECPsph_ine_opt uses `z > 16` to leave the [0,16) tables; floor(16.0/(16/400)) = 400, one row past the
400-row tables (only for z exactly 16). B. `double omega_nuc[CART_CUM]` with CART_CUM = 456 (up to l = 12). Verifier,
PySCF 2.15.1, `mol.intor("ECPscalar_ipnuc")`: l_basis = 5 with lc = 5, 6 works; l = 5, lc = 7 (l+1+lc = 13) segfaults,
consistent with CART_CUM; but l_basis = 6 segfaults for any lc (0, 2, 3, 4), so another limit is involved. ECPscalar
itself works for l = 6. Next: isolate the l = 6 crash (assigned to explorer 2).

Explorer's original record: `agents/explorer2/exploracion3-snapshot-2026-11-22/02-pyscf-ecp-fixed-tables/finding.md` (copied below as `explorer_finding.md`).

## Part B isolated (explorer 2, batch 4, 2026-11-22) and independently verified
**Cause (verified by reading master 484fd7c).** `double omega_nuc[CART_CUM]` (CART_CUM = 456, up to l = 12) is also
declared in `type1_rad_ang` (nr_ecp.c:5253), the local (l = -1) channel, called with lmax = li + lj (nr_ecp.c:6026),
independent of the ECP channels. A derivative raises li by one, so a diagonal block with l = 6 shells gives
lmax = 7 + 6 = 13 and a stack overflow. This explains why l_b = 6 crashes for any lc. With l_b = 5, lc = 7 the
overflow is the non-local `omega_nuc` in `type2_facs_ang` (nr_ecp.c:5078).

**Practical impact (verified, PySCF 2.15.1).** AgH with Ag cc-pV5Z-PP (contains i functions) and the def2-SVP ECP
(local channel present): `ECPscalar` works, `ECPscalar_ipnuc` and `ECPscalar_iprinv` segfault (exit 139). These are
the integrals of the nuclear gradient, so `mf.nuc_grad_method()` cannot run with built-in i-function bases combined
with ECPs that have a local channel (def2, LANL2DZ, SBKJC, CRENBL). Control: def2-QZVPP (l <= 4) with the same ECP
works. Explorer: 27 built-in basis sets contain i functions; patched libcgto reproduces FD4 gradients to 5e-14
(synthetic NaH) and 1.2e-10 Eh/bohr (AgH cc-pV5Z-PP / def2-SVP).

**Additional claims, not independently confirmed.** Under ASan the explorer reports a heap underestimate in
`ECPscalar_cache_size` for second derivatives (SEGV in `ECPscalar_ipipnuc`, l_b = lc = 5) and a heap overflow from the
hard-coded ECP_LMAX = 5 (`ECPscalar_ipnucip`, l_b = 5, lc = 6). In the verifier's ordinary build, called through
`mol.intor(..., comp=9)`, these run without crashing (a first attempt without `comp` crashed only because PySCF
allocated one component instead of nine). A heap overflow need not crash, so this neither confirms nor refutes them;
it needs an ASan build of libcgto. The explorer's patch also covers `_offset_cart[15]`, `_cart_pow_y/z[120]` and
`rad_all` zeroing; for an upstream report the minimal fix for the confirmed defect is enlarging CART_CUM and
`_offset_cart` to ANG_MAX = 15.

Status: B attributed and impact measured (crash of nuclear gradients with real basis/ECP pairs). Ready to report.
