#!/usr/bin/env python3
"""
Phase 1 Steps 4 & 5 (blind verification of OSTRA-008):
1. Boys F_m(T) and SR Boys F_m(T, lower) vs 60-digit mpmath oracle.
2. High-l SR 2e integrals (int2e_sph) vs independent exact identity (full - LR).
3. Binary dependence across 4 builds:
   - clang -O2 (v6.2.1, no quadmath)
   - clang -O0 (v6.2.1, no quadmath)
   - gcc-16 -O2 (v6.2.1, WITH libquadmath)
   - PySCF bundled libcint.6.dylib (v6.1.3)
"""
import ctypes
import json
import math
import subprocess
from pathlib import Path
import numpy as np
import mpmath as mp

mp.mp.dps = 60

ROOT = Path("/Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008")
PYSCF_LIB = Path(
    "/Users/napolesd/Documents/Antigravity/PhysAud/.venv/lib/python3.13/site-packages/"
    "pyscf/lib/deps/lib/libcint.6.dylib"
)

BUILDS = {
    "clang_O2_v621": ROOT / "build-o2/libcint.dylib",
    "clang_O0_v621": ROOT / "build-o0/libcint.dylib",
    "gcc16_O2_quad_v621": ROOT / "build-gcc/libcint.dylib",
    "pyscf_v613": PYSCF_LIB,
}

# libcint constants
ATM_SLOTS = 6
BAS_SLOTS = 8
CHARGE_OF = 0
PTR_COORD = 1
ATOM_OF = 0
ANG_OF = 1
NPRIM_OF = 2
NCTR_OF = 3
PTR_EXP = 5
PTR_COEFF = 6
PTR_RANGE_OMEGA = 8
PTR_ENV_START = 20


def mp_boys(m: int, T: float) -> mp.mpf:
    """Exact Boys function F_m(T) = int_0^1 t^(2m) exp(-T t^2) dt."""
    t_mp = mp.mpf(T)
    if t_mp == 0:
        return mp.mpf(1) / (2 * m + 1)
    s = mp.mpf(m) + mp.mpf("0.5")
    return mp.gammainc(s, 0, t_mp) / (2 * mp.power(t_mp, s))


def mp_sr_boys(m: int, T: float, lower: float) -> mp.mpf:
    """Exact SR Boys function int_lower^1 t^(2m) exp(-T t^2) dt."""
    l_mp = mp.mpf(lower)
    t_mp = mp.mpf(T)
    if t_mp == 0:
        return (1 - mp.power(l_mp, 2 * m + 1)) / (2 * m + 1)
    s = mp.mpf(m) + mp.mpf("0.5")
    # int_lower^1 t^(2m) exp(-T t^2) dt = (gamma(s, T*lower^2) - gamma(s, T)) / (2 T^s)
    # Wait: mp.gammainc(s, a, b) = int_a^b x^(s-1) e^(-x) dx
    return mp.gammainc(s, t_mp * l_mp * l_mp, t_mp) / (2 * mp.power(t_mp, s))


def load_lib(path: Path):
    lib = ctypes.CDLL(str(path))
    lib.gamma_inc_like.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_double, ctypes.c_int]
    lib.gamma_inc_like.restype = None
    lib.fmt_erfc_like.argtypes = [
        ctypes.POINTER(ctypes.c_double), ctypes.c_double, ctypes.c_double, ctypes.c_int
    ]
    lib.fmt_erfc_like.restype = None
    lib.CINTgto_norm.argtypes = [ctypes.c_int, ctypes.c_double]
    lib.CINTgto_norm.restype = ctypes.c_double
    lib.int2e_sph.argtypes = [
        ctypes.c_void_p, ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int), ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p, ctypes.c_void_p,
    ]
    lib.int2e_sph.restype = ctypes.c_int64
    return lib


def call_gamma(lib, T: float, m: int) -> np.ndarray:
    buf = np.zeros(m + 1, dtype=np.float64)
    lib.gamma_inc_like(buf.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), ctypes.c_double(T), ctypes.c_int(m))
    return buf


def call_fmt_erfc(lib, T: float, lower: float, m: int) -> np.ndarray:
    buf = np.zeros(m + 1, dtype=np.float64)
    lib.fmt_erfc_like(
        buf.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_double(T), ctypes.c_double(lower), ctypes.c_int(m)
    )
    return buf


def eval_int2e_sph(lib, ls, R: float, omega: float, exps=(1.2, 0.9, 1.1, 0.8)):
    natm = 2
    nbas = 4
    atm = np.zeros(natm * ATM_SLOTS, dtype=np.int32)
    bas = np.zeros(nbas * BAS_SLOTS, dtype=np.int32)
    env = np.zeros(100, dtype=np.float64)

    off = PTR_ENV_START
    env[PTR_RANGE_OMEGA] = omega
    atm[0 * ATM_SLOTS + CHARGE_OF] = 1
    atm[0 * ATM_SLOTS + PTR_COORD] = off
    env[off:off + 3] = [0.0, 0.0, 0.0]
    off += 3
    atm[1 * ATM_SLOTS + CHARGE_OF] = 1
    atm[1 * ATM_SLOTS + PTR_COORD] = off
    env[off:off + 3] = [0.0, 0.0, R]
    off += 3

    atms = [0, 0, 1, 1]
    for i in range(4):
        bas[i * BAS_SLOTS + ATOM_OF] = atms[i]
        bas[i * BAS_SLOTS + ANG_OF] = ls[i]
        bas[i * BAS_SLOTS + NPRIM_OF] = 1
        bas[i * BAS_SLOTS + NCTR_OF] = 1
        bas[i * BAS_SLOTS + PTR_EXP] = off
        env[off] = exps[i]
        off += 1
        bas[i * BAS_SLOTS + PTR_COEFF] = off
        env[off] = lib.CINTgto_norm(ls[i], exps[i])
        off += 1

    shls = np.array([0, 1, 2, 3], dtype=np.int32)
    shape = tuple(2 * l + 1 for l in ls)
    n = int(np.prod(shape))
    out = np.zeros(n, dtype=np.float64)

    c_atm = atm.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    c_bas = bas.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    c_env = env.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    c_shls = shls.ctypes.data_as(ctypes.POINTER(ctypes.c_int))

    cache_sz = lib.int2e_sph(None, None, c_shls, c_atm, natm, c_bas, nbas, c_env, None, None)
    cache = np.empty(max(int(cache_sz), 1), dtype=np.float64)
    lib.int2e_sph(
        out.ctypes.data_as(ctypes.c_void_p), None, c_shls,
        c_atm, natm, c_bas, nbas, c_env, None,
        cache.ctypes.data_as(ctypes.c_void_p)
    )
    return out


def read_oob_turnover_bytes(lib_path: Path, max_m: int = 52):
    """Inspect TURNOVER_POINT[38..max_m] directly from Mach-O symbol table if present."""
    nm = subprocess.run(["nm", "-n", str(lib_path)], capture_output=True, text=True)
    if nm.returncode != 0:
        return None
    lines = nm.stdout.splitlines()
    tp_addr = None
    next_syms = []
    for idx, line in enumerate(lines):
        if "TURNOVER_POINT" in line:
            parts = line.split()
            tp_addr = int(parts[0], 16)
            for nxt in lines[idx + 1: idx + 8]:
                next_syms.append(nxt)
            break
    if tp_addr is None:
        return {"found": False, "note": "symbol stripped or not found"}
    return {"found": True, "addr": hex(tp_addr), "next_symbols": next_syms}


def probe_effective_turnover(lib, m: int):
    """
    Empirically detect the effective TURNOVER_POINT[m] read by gamma_inc_like!
    How? When T < TURNOVER_POINT[m], gamma_inc_like uses downward power-series
    (fmt1_gamma_inc_like), which is super accurate for small T and overflows/loses
    digits for very large T. When T >= TURNOVER_POINT[m], it uses upward recurrence
    from F_0(T), which suffers catastrophic cancellation for small T (especially T <= 15)!
    More cleanly: for tiny T (e.g. T = 1.0), upward recurrence for m >= 40 gives HUGE error
    (> 1e+10!), whereas downward recurrence gives rel_err < 1e-14!
    So we can classify branch(T) = 'downward' vs 'upward' and binary-search the exact
    switchover value T_eff!
    """
    def is_upward(T):
        if T <= 0:
            return False
        vals = call_gamma(lib, T, m)
        ref_m = float(mp_boys(m, T))
        # At T <= 25, upward recurrence for m >= 39 has rel_err >> 1e-4 (often > 1e5 or NaN/neg!)
        # Wait: what if T_eff > 25? Let's check using F_m(T) monotonicity/exactness against mpmath!
        # Better: check if vals[m] matches upward recurrence computed in Python float64!
        # Upward recurrence in C:
        # b = 0.5 / T; e = 0.5 * exp(-T); f0 = sqrt(pi/4/T)*erf(sqrt(T)); f[i] = b*((2i-1)*f[i-1] - e)
        b = 0.5 / T
        e = 0.5 * math.exp(-T)
        f = math.sqrt(math.pi / (4.0 * T)) * math.erf(math.sqrt(T))
        for i in range(1, m + 1):
            f = b * ((2 * i - 1) * f - e)
        # If gamma_inc_like used upward recurrence, vals[m] == f almost to ULP!
        if math.isnan(vals[m]) or math.isinf(vals[m]):
            return True
        if abs(f) > 0 and abs(vals[m] - f) <= 1e-12 * max(abs(f), abs(vals[m]), 1e-300):
            # Could both match if both accurate (large T)? Check if downward differs slightly or test small T first
            return True
        return False

    # First check at T = 1e-6 (tiny positive T): if upward already -> T_eff <= 1e-6 (e.g., 0.0, negative, or NaN!)
    if is_upward(1e-6):
        # Check if even T = 1e-300 is upward
        if is_upward(1e-300):
            return "<= 1e-300 (0, neg, subnormal, or NaN)"
        lo, hi = 1e-300, 1e-6
    else:
        # Find upper bound where it switches to upward
        lo = 1e-6
        hi = 1.0
        while hi < 1e300 and not is_upward(hi):
            lo = hi
            hi *= 2.0
        if not is_upward(hi):
            return "> 1e300 (+inf or huge)"
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if mid == lo or mid == hi:
            break
        if is_upward(mid):
            hi = mid
        else:
            lo = mid
    return f"{hi:.12g}"


def main():
    libs = {name: load_lib(path) for name, path in BUILDS.items()}
    report = {"symbols": {}, "effective_turnover": {}, "boys_errors": {}, "sr_boys_errors": {}, "eri_sr": {}}

    print("=== 0. Symbol layout after TURNOVER_POINT ===")
    for name, path in BUILDS.items():
        info = read_oob_turnover_bytes(path)
        report["symbols"][name] = info
        print(f"\n[{name}] ({path.name}):")
        if info and info.get("found"):
            print(f"  TURNOVER_POINT @ {info['addr']}")
            for s in info["next_symbols"]:
                print(f"    next: {s}")
        else:
            print(f"  {info}")

    print("\n=== 1. Empirically probed effective TURNOVER_POINT[m] across builds ===")
    ms_probe = [39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 50, 60]
    header = f"{'m':>3} {'true_t0':>10} | " + " | ".join(f"{k:>24}" for k in BUILDS)
    print(header)
    print("-" * len(header))
    for m in ms_probe:
        # true formula from fmt.c:39
        dfact = mp.fac2(2 * m - 1)
        t0_true = float(mp.mpf("0.5") * mp.power(dfact / mp.sqrt(2 * m - 1), mp.mpf(1) / (m - 1)))
        row = {}
        cells = []
        for name, lib in libs.items():
            eff = probe_effective_turnover(lib, m)
            row[name] = eff
            cells.append(f"{eff:>24}")
        report["effective_turnover"][str(m)] = {"true_t0": t0_true, "builds": row}
        print(f"{m:>3} {t0_true:>10.5f} | " + " | ".join(cells))

    print("\n=== 2. Boys F_m(T) max relative error (across k=0..m) vs 60-digit mpmath ===")
    Ts = [0.5, 2.0, 5.0, 10.0, 14.0, 15.0, 16.0, 18.0, 20.0, 25.0, 40.0]
    for m in [39, 40, 41, 43, 45, 47, 50]:
        print(f"\n--- gamma_inc_like m = {m} ---")
        print(f"{'T':>6} | " + " | ".join(f"{k:>20}" for k in BUILDS))
        m_res = {}
        for T in Ts:
            refs = [float(mp_boys(k, T)) for k in range(m + 1)]
            row = {}
            cells = []
            for name, lib in libs.items():
                vals = call_gamma(lib, T, m)
                rel_errs = [abs(vals[k] - refs[k]) / abs(refs[k]) for k in range(m + 1)]
                max_rel = max(rel_errs)
                err_m = rel_errs[m]
                row[name] = {"max_rel": float(max_rel), "rel_m": float(err_m), "val_m": float(vals[m]), "ref_m": refs[m]}
                cells.append(f"{max_rel:>20.3e}")
            m_res[str(T)] = row
            print(f"{T:>6.1f} | " + " | ".join(cells))
        report["boys_errors"][str(m)] = m_res

    print("\n=== 3. SR Boys fmt_erfc_like(T, lower=0.35, m) rel error at order m vs mpmath ===")
    lower = 0.35
    for m in [39, 40, 41, 43, 45, 47]:
        print(f"\n--- fmt_erfc_like m = {m}, lower = {lower} ---")
        print(f"{'T':>6} | " + " | ".join(f"{k:>20}" for k in BUILDS))
        m_res = {}
        for T in [0.5, 2.0, 5.0, 10.0, 14.0, 15.0, 16.0, 18.0, 20.0, 25.0]:
            ref_m = float(mp_sr_boys(m, T, lower))
            row = {}
            cells = []
            for name, lib in libs.items():
                vals = call_fmt_erfc(lib, T, lower, m)
                rel = abs(vals[m] - ref_m) / abs(ref_m)
                row[name] = {"rel_m": float(rel), "val_m": float(vals[m]), "ref_m": ref_m}
                cells.append(f"{rel:>20.3e}")
            m_res[str(T)] = row
            print(f"{T:>6.1f} | " + " | ".join(cells))
        report["sr_boys_errors"][str(m)] = m_res

    print("\n=== 4. High-l SR ERI (int2e_sph) vs exact reference (full - LR) ===")
    # Note: aij = 1.2 + 0.9 = 2.1, akl = 1.1 + 0.8 = 1.9, a0 = 2.1*1.9/4.0 = 0.9975
    # So x = a0 * R^2 ~= R^2!
    # theta = omega^2 / (omega^2 + a0). For omega = -0.5, theta = 0.25 / 1.2475 = 0.2004, lower = sqrt(theta) = 0.4477 >= 0.25 -> ALWAYS hits jacobi!
    # For omega = -0.2, theta = 0.04 / 1.0375 = 0.03855, lower = 0.1964 < 0.25 -> hits jacobi for x <= 60!
    eri_cases = [
        ("ctrl_n20_m39", (10, 10, 10, 9), 0.5, [1.5, 3.0, 4.0]),
        ("oob_n21_m41",  (10, 10, 10, 10), 0.5, [1.0, 2.0, 3.0, 3.8, 4.2, 5.0]),
        ("oob_n21_m41_w02", (10, 10, 10, 10), 0.2, [1.5, 3.0, 4.0, 5.0]),
        ("oob_n22_m43",  (11, 11, 10, 10), 0.5, [1.5, 3.0, 4.0]),
        ("oob_n23_m45",  (12, 11, 11, 10), 0.5, [1.5, 3.0, 4.0]),
        ("oob_n24_m47",  (12, 12, 11, 11), 0.5, [1.5, 3.0, 4.0]),
    ]
    for cname, ls, om, Rs in eri_cases:
        nroots = sum(ls) // 2 + 1
        m_idx = 2 * nroots - 1
        print(f"\n--- {cname}: l={ls} (nroots={nroots}, m={m_idx}, omega=-{om}) ---")
        print(f"{'R':>4} {'x':>6} {'lower':>6} | " + " | ".join(f"{k:>22}" for k in BUILDS))
        c_res = {}
        for R in Rs:
            a0 = (2.1 * 1.9) / 4.0
            x = a0 * R * R
            lower_val = math.sqrt(om * om / (om * om + a0))
            row = {}
            cells = []
            # Use gcc16 full - LR as cross-build gold reference too!
            ref_full_gcc = eval_int2e_sph(libs["gcc16_O2_quad_v621"], ls, R, 0.0)
            ref_lr_gcc = eval_int2e_sph(libs["gcc16_O2_quad_v621"], ls, R, +om)
            ref_sr_gold = ref_full_gcc - ref_lr_gcc
            gold_norm = float(np.linalg.norm(ref_sr_gold))
            for name, lib in libs.items():
                sr = eval_int2e_sph(lib, ls, R, -om)
                full = eval_int2e_sph(lib, ls, R, 0.0)
                lr = eval_int2e_sph(lib, ls, R, +om)
                ref_self = full - lr
                # Check self-consistency of full-LR across builds
                diff_self = float(np.linalg.norm(sr - ref_self)) / max(float(np.linalg.norm(ref_self)), 1e-300)
                diff_gold = float(np.linalg.norm(sr - ref_sr_gold)) / max(gold_norm, 1e-300)
                max_abs_diff = float(np.max(np.abs(sr - ref_sr_gold)))
                max_abs_ref = float(np.max(np.abs(ref_sr_gold)))
                has_nan = bool(np.any(~np.isfinite(sr)))
                row[name] = {
                    "rel_l2_vs_self_full_minus_lr": diff_self,
                    "rel_l2_vs_gcc_full_minus_lr": diff_gold,
                    "max_abs_err": max_abs_diff,
                    "max_abs_ref": max_abs_ref,
                    "nonfinite": has_nan,
                    "full_minus_lr_vs_gcc": float(np.linalg.norm(ref_self - ref_sr_gold)) / max(gold_norm, 1e-300),
                }
                if has_nan:
                    cells.append(f"{'NaN/Inf':>22}")
                else:
                    cells.append(f"{diff_gold:>22.3e}")
            c_res[str(R)] = {"x": x, "lower": lower_val, "gold_norm": gold_norm, "builds": row}
            print(f"{R:>4.1f} {x:>6.2f} {lower_val:>6.3f} | " + " | ".join(cells))
        report["eri_sr"][cname] = c_res

    out_json = ROOT / "outputs/phase1_oracle_results.json"
    out_json.write_text(json.dumps(report, indent=2))
    print(f"\nSaved JSON report to {out_json}")


if __name__ == "__main__":
    main()
