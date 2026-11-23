#!/usr/bin/env python3
"""
Phase 2 evaluation of proposed patch for OSTRA-008:
Extend TURNOVER_POINT up to m = 71 using t0 = 0.5 * ((2m-1)!! / sqrt(2m-1))**(1/(m-1)).
Checks:
 (b) Exact comparison of formula vs existing 40 entries (m=0..39) in src/fmt.c.
 (c) Build isolated patched copies (clang ASan, clang -O2, gcc-16 -O2 quadmath),
     verify ASan clean, verify Boys F_m & SR ERIs vs references, verify low-l bit-identity.
"""
import json
import math
import re
import subprocess
import sys
from pathlib import Path
import numpy as np
import mpmath as mp

mp.mp.dps = 60

ROOT = Path("/Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008")
ORIG_FMT = ROOT / "libcint/src/fmt.c"
PATCHED_DIR = ROOT / "patched"
PATCHED_FMT = PATCHED_DIR / "fmt.c"

sys.path.insert(0, str(ROOT / "scripts"))
from build_libcint import build_target
from phase1_oracle_and_builds import load_lib, call_gamma, call_fmt_erfc, eval_int2e_sph, mp_boys, mp_sr_boys


def formula_t0(m: int) -> mp.mpf:
    if m <= 1:
        return mp.mpf(0)
    dfact = mp.fac2(2 * m - 1)
    return mp.mpf("0.5") * mp.power(dfact / mp.sqrt(2 * m - 1), mp.mpf(1) / (m - 1))


def parse_existing_table(fmt_text: str):
    m = re.search(r"static double TURNOVER_POINT\[\]\s*=\s*\{([^}]+)\};", fmt_text, re.S)
    assert m is not None
    entries = [x.strip() for x in m.group(1).split(",") if x.strip()]
    return entries, m.start(), m.end()


def main():
    fmt_text = ORIG_FMT.read_text()
    existing_strs, start_pos, end_pos = parse_existing_table(fmt_text)
    assert len(existing_strs) == 40, f"Expected 40 entries, got {len(existing_strs)}"

    print("=== Phase 2(b): Formula vs 40 existing TURNOVER_POINT entries (m=0..39) ===")
    max_abs_diff = mp.mpf(0)
    max_abs_m = -1
    max_rel_diff = mp.mpf(0)
    max_rel_m = -1
    b_rows = []
    for m_idx, s in enumerate(existing_strs):
        val_table = mp.mpf(s)
        val_formula = formula_t0(m_idx)
        abs_d = abs(val_table - val_formula)
        rel_d = abs_d / val_formula if val_formula != 0 else mp.mpf(0)
        b_rows.append({
            "m": m_idx,
            "table_str": s,
            "formula_15g": mp.nstr(val_formula, 15),
            "abs_diff": float(abs_d),
            "rel_diff": float(rel_d),
        })
        if m_idx >= 2:
            if abs_d > max_abs_diff:
                max_abs_diff = abs_d
                max_abs_m = m_idx
            if rel_d > max_rel_diff:
                max_rel_diff = rel_d
                max_rel_m = m_idx

    print(f"Existing entries count: {len(existing_strs)} (indices 0..39)")
    print(f"m=0,1: table={existing_strs[0]}, {existing_strs[1]} (formula singular at m=1, undefined at m=0; defined as 0.)")
    print(f"m=2..39 max abs diff: {float(max_abs_diff):.3e} (at m={max_abs_m}: table={existing_strs[max_abs_m]} vs formula={mp.nstr(formula_t0(max_abs_m), 15)})")
    print(f"m=2..39 max rel diff: {float(max_rel_diff):.3e} (at m={max_rel_m}: table={existing_strs[max_rel_m]} vs formula={mp.nstr(formula_t0(max_rel_m), 15)})")
    print("Note: existing table entries have 12-13 significant digits (12 decimal places for m<10, 11 for m>=10).")
    for r in b_rows:
        if r["abs_diff"] > 5e-12:
            print(f"  notable diff m={r['m']}: table={r['table_str']} formula={r['formula_15g']} abs={r['abs_diff']:.3e}")

    # Create isolated patched fmt.c extending m=40..71
    PATCHED_DIR.mkdir(parents=True, exist_ok=True)
    new_lines = [f"        {s}," for s in existing_strs]
    for m_idx in range(40, 72):
        val = formula_t0(m_idx)
        # Format to 13 significant digits (matching existing 11 decimal places for 10 <= t0 < 100)
        formatted = mp.nstr(val, 13)
        comma = "," if m_idx < 71 else ""
        new_lines.append(f"        {formatted}{comma}")
    new_table = "static double TURNOVER_POINT[] = {\n" + "\n".join(new_lines) + "\n};"
    patched_text = fmt_text[:start_pos] + new_table + fmt_text[end_pos:]
    PATCHED_FMT.write_text(patched_text)
    print(f"\nWrote isolated patched fmt.c to {PATCHED_FMT} (72 entries: m=0..71)")
    print(f"  m=40: {mp.nstr(formula_t0(40), 13)}")
    print(f"  m=47: {mp.nstr(formula_t0(47), 13)}")
    print(f"  m=63: {mp.nstr(formula_t0(63), 13)}")
    print(f"  m=71: {mp.nstr(formula_t0(71), 13)}")

    print("\n=== Phase 2(c): Building patched libraries ===")
    build_target(
        ROOT / "build-patched-asan",
        cc="clang",
        cflags=["-O1", "-g", "-fsanitize=address", "-fno-omit-frame-pointer"],
        ldflags=["-fsanitize=address"],
        fmt_override=PATCHED_FMT,
        libname="libcint_asan.dylib",
    )
    build_target(
        ROOT / "build-patched-o2",
        cc="clang",
        cflags=["-O2", "-g"],
        ldflags=[],
        fmt_override=PATCHED_FMT,
        libname="libcint.dylib",
    )
    build_target(
        ROOT / "build-patched-gcc",
        cc="/opt/homebrew/bin/gcc-16",
        cflags=["-O2", "-g"],
        ldflags=["-lquadmath"],
        have_quadmath=True,
        fmt_override=PATCHED_FMT,
        libname="libcint.dylib",
    )

    # Compile and run patched ASan test executable (extended with m=71 and m=72 checks)
    asan_c = ROOT / "scripts/test_asan_patched.c"
    asan_c.write_text(
        (ROOT / "scripts/test_asan.c").read_text().replace(
            '    } else if (strcmp(mode, "eri_sr_3c2e") == 0) {\n        run_eri_3c2e(14, 14, 12, -0.5);',
            '    } else if (strcmp(mode, "eri_sr_3c2e") == 0) {\n        run_eri_3c2e(14, 14, 12, -0.5);\n'
            '    } else if (strcmp(mode, "boys_gamma_71") == 0) {\n        gamma_inc_like(f, 10.0, 71); printf("m=71 ok: %.12e\\n", f[71]);\n'
            '    } else if (strcmp(mode, "boys_gamma_72") == 0) {\n        gamma_inc_like(f, 10.0, 72); printf("m=72 ok: %.12e\\n", f[72]);'
        )
    )
    asan_bin = ROOT / "build-patched-asan/test_asan_patched"
    subprocess.run([
        "clang", "-O1", "-g", "-fsanitize=address", "-fno-omit-frame-pointer",
        f"-I{ROOT / 'build-patched-asan/include'}",
        f"-I{ROOT / 'libcint/include'}",
        f"-I{ROOT / 'build-patched-asan/src'}",
        f"-I{ROOT / 'libcint/src'}",
        str(asan_c),
        str(ROOT / "build-patched-asan/libcint_asan.dylib"),
        f"-Wl,-rpath,{ROOT / 'build-patched-asan'}",
        "-lm", "-o", str(asan_bin),
    ], check=True)

    print("\n=== Phase 2(c.1): ASan verification on patched build ===")
    asan_results = {}
    asan_log_lines = []
    for case in [
        "boys_ok", "boys_gamma_40", "boys_lgamma_40", "boys_erfc_40",
        "boys_lerfc_41", "boys_gamma_71", "eri_sr_l9", "eri_sr_l10",
        "eri_full_l10", "eri_lr_l10", "eri_sr_3c2e", "boys_gamma_72",
    ]:
        r = subprocess.run([str(asan_bin), case], capture_output=True, text=True)
        asan_results[case] = r.returncode
        asan_log_lines.append(f"===== CASE: {case} (EXIT={r.returncode}) =====\n{r.stdout}{r.stderr}")
        print(f"  {case:16s} -> EXIT={r.returncode}")
    (ROOT / "outputs/asan_phase2_patched.log").write_text("\n".join(asan_log_lines))

    # Load unpatched & patched libs for numerical verification
    orig_clang = load_lib(ROOT / "build-o2/libcint.dylib")
    pat_clang = load_lib(ROOT / "build-patched-o2/libcint.dylib")
    orig_gcc = load_lib(ROOT / "build-gcc/libcint.dylib")
    pat_gcc = load_lib(ROOT / "build-patched-gcc/libcint.dylib")

    print("\n=== Phase 2(c.2): Boys F_m(T) patched vs 60-digit mpmath ===")
    boys_patched = {}
    for m_idx in [40, 41, 45, 47, 50, 60, 71]:
        t0 = float(formula_t0(m_idx))
        Ts = [0.5, 2.0, 5.0, 10.0, 14.0, round(t0 - 0.5, 2), round(t0 + 0.5, 2), 30.0, 40.0]
        row = {}
        print(f"\n--- gamma_inc_like m={m_idx} (t0={t0:.4f}) ---")
        for T in Ts:
            ref_m = float(mp_boys(m_idx, T))
            v_orig = call_gamma(orig_gcc, T, m_idx)[m_idx]
            v_pat_c = call_gamma(pat_clang, T, m_idx)[m_idx]
            v_pat_g = call_gamma(pat_gcc, T, m_idx)[m_idx]
            e_orig = abs(v_orig - ref_m) / abs(ref_m)
            e_pat_c = abs(v_pat_c - ref_m) / abs(ref_m)
            e_pat_g = abs(v_pat_g - ref_m) / abs(ref_m)
            row[str(T)] = {"orig_gcc": e_orig, "patched_clang": e_pat_c, "patched_gcc": e_pat_g}
            print(f"  T={T:6.2f}: orig_gcc={e_orig:10.3e}  patched_clang={e_pat_c:10.3e}  patched_gcc={e_pat_g:10.3e}")
        boys_patched[str(m_idx)] = row

    print("\n=== Phase 2(c.3): High-l SR ERI (int2e_sph) patched vs exact reference (full - LR) ===")
    eri_patched = {}
    cases = [
        ("oob_n21_m41", (10, 10, 10, 10), 0.5, [0.2, 0.5, 0.6, 0.8, 1.0, 2.0, 4.0]),
        ("oob_n22_m43", (11, 11, 10, 10), 0.5, [0.2, 0.6, 0.8, 1.0, 2.0, 4.0]),
        ("oob_n23_m45", (12, 11, 11, 10), 0.5, [0.2, 0.8, 1.0, 3.0]),
        ("oob_n24_m47", (12, 12, 11, 11), 0.5, [0.2, 0.8, 1.0, 3.0]),
    ]
    for cname, ls, om, Rs in cases:
        print(f"\n--- {cname}: l={ls} omega=-{om} ---")
        print(f"{'R':>4} {'x':>6} | {'orig_clang':>12} | {'pat_clang':>12} | {'orig_gcc_q':>12} | {'pat_gcc_q':>12}")
        c_res = {}
        for R in Rs:
            x = (2.1 * 1.9 / 4.0) * R * R
            ref = eval_int2e_sph(pat_gcc, ls, R, 0.0) - eval_int2e_sph(pat_gcc, ls, R, +om)
            ref_n = float(np.linalg.norm(ref))
            vals = {}
            for tag, lib in [("orig_clang", orig_clang), ("pat_clang", pat_clang), ("orig_gcc_q", orig_gcc), ("pat_gcc_q", pat_gcc)]:
                sr = eval_int2e_sph(lib, ls, R, -om)
                vals[tag] = float(np.linalg.norm(sr - ref)) / ref_n
            c_res[str(R)] = {"x": x, **vals}
            print(f"{R:4.1f} {x:6.3f} | {vals['orig_clang']:12.3e} | {vals['pat_clang']:12.3e} | {vals['orig_gcc_q']:12.3e} | {vals['pat_gcc_q']:12.3e}")
        eri_patched[cname] = c_res

    print("\n=== Phase 2(c.4): Bit-for-bit identity check for low-l and m <= 39 ===")
    bit_checks = {"boys_m_le_39": {}, "eri_low_l": {}}
    # Boys m = 0..39
    for label, u_lib, p_lib in [("clang_O2", orig_clang, pat_clang), ("gcc16_quad", orig_gcc, pat_gcc)]:
        max_diff_g = 0.0
        max_diff_e = 0.0
        bit_eq = True
        for m_idx in range(40):
            for T in [0.0, 1e-8, 0.5, 2.0, 5.0, 10.0, 14.9, 15.0, 25.0, 50.0]:
                g_u = call_gamma(u_lib, T, m_idx)
                g_p = call_gamma(p_lib, T, m_idx)
                e_u = call_fmt_erfc(u_lib, T, 0.35, m_idx)
                e_p = call_fmt_erfc(p_lib, T, 0.35, m_idx)
                max_diff_g = max(max_diff_g, float(np.max(np.abs(g_u - g_p))))
                max_diff_e = max(max_diff_e, float(np.max(np.abs(e_u - e_p))))
                if not (np.array_equal(g_u.view(np.uint64), g_p.view(np.uint64)) and
                        np.array_equal(e_u.view(np.uint64), e_p.view(np.uint64))):
                    bit_eq = False
        bit_checks["boys_m_le_39"][label] = {"bit_identical": bit_eq, "max_abs_gamma": max_diff_g, "max_abs_erfc": max_diff_e}
        print(f"  Boys m=0..39 [{label}]: bit_identical={bit_eq} max_abs_gamma={max_diff_g} max_abs_erfc={max_diff_e}")

    # ERIs for l = 0..5 + control (10,10|10,9)
    for label, u_lib, p_lib in [("clang_O2", orig_clang, pat_clang), ("gcc16_quad", orig_gcc, pat_gcc)]:
        all_bit_eq = True
        max_abs = 0.0
        details = []
        for ls in [(0,0,0,0), (1,1,1,1), (2,2,2,2), (3,3,3,3), (4,4,4,4), (5,5,5,5), (10,10,10,9)]:
            for om in [0.0, 0.5, -0.5]:
                for R in [0.5, 1.5, 3.5]:
                    a = eval_int2e_sph(u_lib, ls, R, om)
                    b = eval_int2e_sph(p_lib, ls, R, om)
                    d = float(np.max(np.abs(a - b)))
                    eq = bool(np.array_equal(a.view(np.uint64), b.view(np.uint64)))
                    max_abs = max(max_abs, d)
                    if not eq:
                        all_bit_eq = False
                        details.append({"ls": ls, "omega": om, "R": R, "max_abs": d})
        bit_checks["eri_low_l"][label] = {"bit_identical": all_bit_eq, "max_abs": max_abs, "mismatches": details}
        print(f"  ERIs low-l & ctrl(m=39) [{label}]: bit_identical={all_bit_eq} max_abs={max_abs} mismatches={len(details)}")

    out_json = ROOT / "outputs/phase2_results.json"
    out_json.write_text(json.dumps({
        "phase2b_formula_vs_table": {
            "max_abs_diff": float(max_abs_diff), "max_abs_m": max_abs_m,
            "max_rel_diff": float(max_rel_diff), "max_rel_m": max_rel_m,
            "rows": b_rows,
        },
        "phase2c_asan": asan_results,
        "phase2c_boys": boys_patched,
        "phase2c_eri": eri_patched,
        "phase2c_bit_identity": bit_checks,
    }, indent=2))
    print(f"\nSaved Phase 2 results to {out_json}")


if __name__ == "__main__":
    main()
