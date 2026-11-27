#!/usr/bin/env python3
"""
Phase 1 (blind) verification of compiled baseline vs patched libxcfun.dylib:
  1. Built-in xcfun_test() on baseline and patched libraries
  2. Patched 'tw' vs exact H 1s (0.5 Eh), H2O, polarized point, libxc GGA_K_VW, and 4-pt FD
  3. Patched 'spbec' vs H2O, polarized point, t->inf limit, libxc GGA_C_SPBE, and 4-pt FD
  4. Exhaustive bit-for-bit identity check (exact IEEE-754 64-bit representation)
     across ALL OTHER functionals in XCFun (excluding tw and spbec).
"""
import ctypes
import struct
import numpy as np
from pyscf import gto, dft

BASE = "/Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-027-028"

# xcfun enums from api/xcfun.h
XC_PARTIAL_DERIVATIVES = 1
XC_A_B_GAA_GAB_GBB = 6
XC_A_B_GAA_GAB_GBB_LAPA_LAPB_TAUA_TAUB = 16


def load_xcfun(path):
    lib = ctypes.CDLL(path)
    lib.xcfun_test.restype = ctypes.c_int
    lib.xcfun_new.restype = ctypes.c_void_p
    lib.xcfun_delete.argtypes = [ctypes.c_void_p]
    lib.xcfun_set.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_double]
    lib.xcfun_set.restype = ctypes.c_int
    lib.xcfun_eval_setup.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_int]
    lib.xcfun_eval_setup.restype = ctypes.c_int
    lib.xcfun_input_length.argtypes = [ctypes.c_void_p]
    lib.xcfun_input_length.restype = ctypes.c_int
    lib.xcfun_output_length.argtypes = [ctypes.c_void_p]
    lib.xcfun_output_length.restype = ctypes.c_int
    lib.xcfun_eval.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
    ]
    lib.xcfun_eval_vec.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
    ]
    lib.xcfun_enumerate_parameters.argtypes = [ctypes.c_int]
    lib.xcfun_enumerate_parameters.restype = ctypes.c_char_p
    return lib


def eval_gga_vec(lib, name, d_arr, order=1):
    """Evaluate GGA functional on array d_arr of shape (npts, 5): [a, b, gaa, gab, gbb]."""
    fun = lib.xcfun_new()
    rc = lib.xcfun_set(fun, name.encode("ascii"), 1.0)
    assert rc == 0, f"xcfun_set({name}) failed: {rc}"
    rc = lib.xcfun_eval_setup(fun, XC_A_B_GAA_GAB_GBB, XC_PARTIAL_DERIVATIVES, order)
    assert rc == 0, f"xcfun_eval_setup({name}) failed: {rc}"
    nout = lib.xcfun_output_length(fun)
    d_c = np.ascontiguousarray(d_arr, dtype=np.float64)
    npts = d_c.shape[0]
    out = np.zeros((npts, nout), dtype=np.float64)
    lib.xcfun_eval_vec(
        fun,
        npts,
        d_c.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        5,
        out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        nout,
    )
    lib.xcfun_delete(fun)
    return out


lib_base = load_xcfun(f"{BASE}/build-base/libxcfun.dylib")
lib_pat = load_xcfun(f"{BASE}/build-patched/libxcfun.dylib")

print("=" * 80)
print("1. BUILT-IN XCFUN SELF-TESTS (xcfun_test)")
print("=" * 80)
nf_base = lib_base.xcfun_test()
nf_pat = lib_pat.xcfun_test()
print(f"Baseline xcfun_test() failed count: {nf_base}")
print(f"Patched  xcfun_test() failed count: {nf_pat}")

print("\n" + "=" * 80)
print("2. PATCHED 'tw' VERIFICATION (H 1s, H2O, polarized point, 4-pt FD)")
print("=" * 80)
# 2a. Exact H 1s on PySCF grid
mol_h = gto.M(atom="H 0 0 0", basis="cc-pvqz", spin=1, verbose=0)
grids_h = dft.gen_grid.Grids(mol_h)
grids_h.level = 5
grids_h.build()
r = np.linalg.norm(grids_h.coords, axis=1)
rho_a = (1.0 / np.pi) * np.exp(-2.0 * r)
gaa = 4.0 * rho_a**2
d_h = np.column_stack([rho_a, np.zeros_like(rho_a), gaa, np.zeros_like(rho_a), np.zeros_like(rho_a)])
out_h_base = eval_gga_vec(lib_base, "tw", d_h, order=0)
out_h_pat = eval_gga_vec(lib_pat, "tw", d_h, order=0)
E_h_base = float(np.sum(grids_h.weights * out_h_base[:, 0]))
E_h_pat = float(np.sum(grids_h.weights * out_h_pat[:, 0]))
print(f"Exact H 1s T_W (exact = 0.5 Eh):")
print(f"  Baseline tw : {E_h_base:.15f} Eh (err = {E_h_base - 0.5:+.6e})")
print(f"  Patched  tw : {E_h_pat:.15f} Eh (err = {E_h_pat - 0.5:+.3e})")

# 2b. H2O molecular density
mol_w = gto.M(atom="O 0 0 0; H 0 0.757 0.587; H 0 -0.757 0.587", basis="cc-pvtz", verbose=0)
mf_w = dft.RKS(mol_w)
mf_w.xc = "PBE"
mf_w.kernel()
dm_w = mf_w.make_rdm1()
grids_w = mf_w.grids
ni = dft.numint.NumInt()
ao = ni.eval_ao(mol_w, grids_w.coords, deriv=1)
rho_w = ni.eval_rho(mol_w, ao, dm_w, xctype="GGA")
n_w = rho_w[0]
gnn_w = rho_w[1]**2 + rho_w[2]**2 + rho_w[3]**2
d_w = np.column_stack([0.5 * n_w, 0.5 * n_w, 0.25 * gnn_w, 0.25 * gnn_w, 0.25 * gnn_w])

E_w_tw_base = float(np.sum(grids_w.weights * eval_gga_vec(lib_base, "tw", d_w, order=0)[:, 0]))
E_w_tw_pat = float(np.sum(grids_w.weights * eval_gga_vec(lib_pat, "tw", d_w, order=0)[:, 0]))
_, E_w_vw_libxc, _ = ni.nr_rks(mol_w, grids_w, "GGA_K_VW,", dm_w)
print(f"\nH2O RKS PBE/cc-pVTZ:")
print(f"  Baseline tw       : {E_w_tw_base:.12e} Eh")
print(f"  Patched  tw       : {E_w_tw_pat:.12f} Eh")
print(f"  libxc GGA_K_VW    : {E_w_vw_libxc:.12f} Eh")
print(f"  Diff patched-libxc: {E_w_tw_pat - E_w_vw_libxc:+.3e} Eh")

# 2c. Polarized point + derivatives + 4-pt FD
na_p, nb_p = 0.35, 0.15
ga_p = np.array([0.4, 0.1, 0.0])
gb_p = np.array([-0.1, 0.3, 0.2])
gaa_p = float(np.dot(ga_p, ga_p))
gab_p = float(np.dot(ga_p, gb_p))
gbb_p = float(np.dot(gb_p, gb_p))
d_p = np.array([[na_p, nb_p, gaa_p, gab_p, gbb_p]])
rho_pol = (
    np.array([[na_p], [ga_p[0]], [ga_p[1]], [ga_p[2]]]),
    np.array([[nb_p], [gb_p[0]], [gb_p[1]], [gb_p[2]]]),
)

out_tw_pat = eval_gga_vec(lib_pat, "tw", d_p, order=1)[0]
exc_vw_lib, vxc_vw_lib, _, _ = dft.libxc.eval_xc("GGA_K_VW,", rho_pol, spin=1, deriv=1)
ref_tw_vec = np.array([
    (na_p + nb_p) * exc_vw_lib[0],
    vxc_vw_lib[0][0, 0], vxc_vw_lib[0][0, 1],
    vxc_vw_lib[1][0, 0], vxc_vw_lib[1][0, 1], vxc_vw_lib[1][0, 2],
])
print(f"\nPolarized point (na={na_p}, nb={nb_p}) patched tw vs libxc GGA_K_VW:")
print(f"  [e, va, vb, vgaa, vgab, vgbb] max|diff| = {np.max(np.abs(out_tw_pat - ref_tw_vec)):.3e}")

for h in [1e-2, 1e-3, 1e-4, 1e-5]:
    fd = np.zeros(5)
    for k in range(5):
        dp1 = d_p.copy(); dp1[0, k] += h
        dm1 = d_p.copy(); dm1[0, k] -= h
        dp2 = d_p.copy(); dp2[0, k] += 2 * h
        dm2 = d_p.copy(); dm2[0, k] -= 2 * h
        e_p1 = eval_gga_vec(lib_pat, "tw", dp1, order=0)[0, 0]
        e_m1 = eval_gga_vec(lib_pat, "tw", dm1, order=0)[0, 0]
        e_p2 = eval_gga_vec(lib_pat, "tw", dp2, order=0)[0, 0]
        e_m2 = eval_gga_vec(lib_pat, "tw", dm2, order=0)[0, 0]
        fd[k] = (-e_p2 + 8 * e_p1 - 8 * e_m1 + e_m2) / (12 * h)
    print(f"  4pt-FD tw h={h:.0e}: max|deriv - FD| = {np.max(np.abs(out_tw_pat[1:] - fd)):.3e}")

print("\n" + "=" * 80)
print("3. PATCHED 'spbec' VERIFICATION (H2O, polarized point, t->inf limit, 4-pt FD)")
print("=" * 80)
E_w_sp_base = float(np.sum(grids_w.weights * eval_gga_vec(lib_base, "spbec", d_w, order=0)[:, 0]))
E_w_sp_pat = float(np.sum(grids_w.weights * eval_gga_vec(lib_pat, "spbec", d_w, order=0)[:, 0]))
_, E_w_sp_libxc, _ = ni.nr_rks(mol_w, grids_w, ",GGA_C_SPBE", dm_w)
print(f"H2O RKS PBE/cc-pVTZ:")
print(f"  Baseline spbec    : {E_w_sp_base:+.12f} Eh (diff vs libxc = {E_w_sp_base - E_w_sp_libxc:+.6e})")
print(f"  Patched  spbec    : {E_w_sp_pat:+.12f} Eh")
print(f"  libxc GGA_C_SPBE  : {E_w_sp_libxc:+.12f} Eh")
print(f"  Diff patched-libxc: {E_w_sp_pat - E_w_sp_libxc:+.3e} Eh")

out_sp_pat = eval_gga_vec(lib_pat, "spbec", d_p, order=1)[0]
exc_sp_lib, vxc_sp_lib, _, _ = dft.libxc.eval_xc(",GGA_C_SPBE", rho_pol, spin=1, deriv=1)
ref_sp_vec = np.array([
    (na_p + nb_p) * exc_sp_lib[0],
    vxc_sp_lib[0][0, 0], vxc_sp_lib[0][0, 1],
    vxc_sp_lib[1][0, 0], vxc_sp_lib[1][0, 1], vxc_sp_lib[1][0, 2],
])
print(f"\nPolarized point (na={na_p}, nb={nb_p}) patched spbec vs libxc GGA_C_SPBE:")
print(f"  Patched [e, va, vb, vgaa, vgab, vgbb]: {out_sp_pat}")
print(f"  libxc   [e, va, vb, vgaa, vgab, vgbb]: {ref_sp_vec}")
print(f"  max|diff| = {np.max(np.abs(out_sp_pat - ref_sp_vec)):.3e}")

for h in [1e-2, 1e-3, 1e-4, 1e-5]:
    fd = np.zeros(5)
    for k in range(5):
        dp1 = d_p.copy(); dp1[0, k] += h
        dm1 = d_p.copy(); dm1[0, k] -= h
        dp2 = d_p.copy(); dp2[0, k] += 2 * h
        dm2 = d_p.copy(); dm2[0, k] -= 2 * h
        e_p1 = eval_gga_vec(lib_pat, "spbec", dp1, order=0)[0, 0]
        e_m1 = eval_gga_vec(lib_pat, "spbec", dm1, order=0)[0, 0]
        e_p2 = eval_gga_vec(lib_pat, "spbec", dp2, order=0)[0, 0]
        e_m2 = eval_gga_vec(lib_pat, "spbec", dm2, order=0)[0, 0]
        fd[k] = (-e_p2 + 8 * e_p1 - 8 * e_m1 + e_m2) / (12 * h)
    print(f"  4pt-FD spbec h={h:.0e}: max|deriv - FD| = {np.max(np.abs(out_sp_pat[1:] - fd)):.3e}")

print("\nRapidly varying limit (t -> inf) on patched spbec:")
for grad_mag in [1.0, 10.0, 100.0, 1e4]:
    gnn_lt = grad_mag**2
    d_lt = np.array([[0.05, 0.05, 0.25 * gnn_lt, 0.25 * gnn_lt, 0.25 * gnn_lt]])
    e_b = eval_gga_vec(lib_base, "spbec", d_lt, order=0)[0, 0] / 0.1
    e_p = eval_gga_vec(lib_pat, "spbec", d_lt, order=0)[0, 0] / 0.1
    print(f"  n=0.1, |grad n|={grad_mag:7.1f}: base eps_c={e_b:+.10e} | patched eps_c={e_p:+.10e}")

print("\n" + "=" * 80)
print("4. BIT-FOR-BIT IDENTITY CHECK ACROSS ALL OTHER XCFUN FUNCTIONALS (ORDER 2)")
print("=" * 80)
# Full meta-GGA polarized test point: [a, b, gaa, gab, gbb, lapa, lapb, taua, taub]
d_mgga = np.array([0.35, 0.15, 0.17, -0.01, 0.14, 0.22, -0.11, 0.45, 0.28], dtype=np.float64)

all_params = []
idx = 0
while True:
    p = lib_base.xcfun_enumerate_parameters(idx)
    if not p:
        break
    all_params.append(p.decode("ascii"))
    idx += 1

checked = 0
skipped_params = []
mismatches = []
total_doubles = 0

for name in all_params:
    if name in ("tw", "spbec"):
        continue
    # Skip non-functional scalar settings (cam_alpha, cam_beta, rangesep_mu, etc.)
    fun_b = lib_base.xcfun_new()
    fun_p = lib_pat.xcfun_new()
    lib_base.xcfun_set(fun_b, name.encode("ascii"), 1.0)
    lib_pat.xcfun_set(fun_p, name.encode("ascii"), 1.0)
    rc_b = lib_base.xcfun_eval_setup(fun_b, XC_A_B_GAA_GAB_GBB_LAPA_LAPB_TAUA_TAUB, XC_PARTIAL_DERIVATIVES, 2)
    rc_p = lib_pat.xcfun_eval_setup(fun_p, XC_A_B_GAA_GAB_GBB_LAPA_LAPB_TAUA_TAUB, XC_PARTIAL_DERIVATIVES, 2)
    if rc_b != 0 or rc_p != 0:
        skipped_params.append((name, rc_b, rc_p))
        lib_base.xcfun_delete(fun_b)
        lib_pat.xcfun_delete(fun_p)
        continue
    nout = lib_base.xcfun_output_length(fun_b)
    out_b = np.zeros(nout, dtype=np.float64)
    out_p = np.zeros(nout, dtype=np.float64)
    lib_base.xcfun_eval(
        fun_b,
        d_mgga.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        out_b.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    )
    lib_pat.xcfun_eval(
        fun_p,
        d_mgga.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        out_p.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    )
    lib_base.xcfun_delete(fun_b)
    lib_pat.xcfun_delete(fun_p)

    # Exact bit-for-bit IEEE-754 comparison
    bytes_b = out_b.tobytes()
    bytes_p = out_p.tobytes()
    checked += 1
    total_doubles += nout
    if bytes_b != bytes_p:
        mismatches.append((name, float(np.max(np.abs(out_b - out_p)))))

print(f"Total parameters enumerated   : {len(all_params)}")
print(f"Other functionals bit-checked : {checked} ({total_doubles} IEEE-754 doubles, order 0..2)")
print(f"Non-functional params skipped : {len(skipped_params)} -> {[x[0] for x in skipped_params]}")
print(f"Bit-for-bit mismatches found  : {len(mismatches)}")
if mismatches:
    for m in mismatches:
        print("  MISMATCH:", m)
else:
    print("RESULT: 100% BIT-FOR-BIT IDENTICAL (0 ULP, 0 byte diffs) across all other functionals!")
