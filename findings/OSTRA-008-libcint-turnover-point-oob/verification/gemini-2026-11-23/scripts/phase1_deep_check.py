#!/usr/bin/env python3
"""
Phase 1 deep check (Steps 4 & 5):
1. Dump exact raw IEEE-754 double values at TURNOVER_POINT[38..52] from the Mach-O binary
   of all 4 builds (clang -O2, clang -O0, gcc-16 -O2 quadmath, PySCF v6.1.3).
2. Sweep small-to-medium R (x = 0.01 .. 16) for high-l SR int2e_sph across all 4 builds,
   specifically demonstrating:
   - Catastrophic silent error EVEN on gcc-16 with 128-bit __float128!
   - Direct binary-dependent divergence on int2e_sph at x in (0, 0.5238) where clang -O0
     reads JACOBI_ALPHA[1]=0.5238 while clang -O2 / pyscf read subnormals ~1e-314!
   - Control l=(10,10|10,9) (nroots=20, m=39) remaining accurate at the exact same small R!
"""
import ctypes
import json
import math
import struct
import subprocess
from pathlib import Path
import numpy as np

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


def read_macho_doubles_at_symbol(path: Path, sym_substr: str, count: int):
    """Find virtual address of symbol via nm, map VM addr -> file offset via otool -l, read doubles."""
    nm = subprocess.run(["nm", "-n", str(path)], capture_output=True, text=True, check=True)
    vmaddr = None
    for line in nm.stdout.splitlines():
        if sym_substr in line:
            vmaddr = int(line.split()[0], 16)
            break
    if vmaddr is None:
        return None

    ot = subprocess.run(["otool", "-l", str(path)], capture_output=True, text=True, check=True)
    # Parse LC_SEGMENT_64 sections: addr, size, offset
    cur_addr = cur_size = cur_off = None
    file_off = None
    for line in ot.stdout.splitlines():
        s = line.strip()
        if s.startswith("sectname"):
            cur_addr = cur_size = cur_off = None
        elif s.startswith("addr "):
            cur_addr = int(s.split()[1], 16)
        elif s.startswith("size "):
            cur_size = int(s.split()[1], 16)
        elif s.startswith("offset "):
            cur_off = int(s.split()[1])
            if cur_addr is not None and cur_size is not None:
                if cur_addr <= vmaddr < cur_addr + cur_size:
                    file_off = cur_off + (vmaddr - cur_addr)
                    break
    if file_off is None:
        return None
    data = path.read_bytes()[file_off : file_off + 8 * count]
    vals = struct.unpack(f"<{count}d", data)
    raw_u64 = struct.unpack(f"<{count}Q", data)
    return vmaddr, file_off, vals, raw_u64


# Import helpers from phase1_oracle_and_builds
import sys
sys.path.insert(0, str(ROOT / "scripts"))
from phase1_oracle_and_builds import load_lib, eval_int2e_sph


def main():
    print("=== A. Exact IEEE-754 doubles sitting at TURNOVER_POINT[38..51] in each binary ===")
    dumps = {}
    for name, path in BUILDS.items():
        res = read_macho_doubles_at_symbol(path, "TURNOVER_POINT", 52)
        assert res is not None, f"Could not map TURNOVER_POINT in {name}"
        vmaddr, foff, vals, raw_u64 = res
        dumps[name] = {
            "vmaddr": hex(vmaddr),
            "file_offset": hex(foff),
            "entries": {str(m): {"double": vals[m], "hex": f"0x{raw_u64[m]:016x}"} for m in range(38, 52)},
        }

    header = f"{'m':>3} | " + " | ".join(f"{k:>22}" for k in BUILDS)
    print(header)
    print("-" * len(header))
    for m in range(38, 52):
        cells = []
        for name in BUILDS:
            v = dumps[name]["entries"][str(m)]["double"]
            cells.append(f"{v:>22.6e}")
        tag = "  <-- IN BOUNDS" if m < 40 else ("  <-- FIRST OOB" if m == 40 else "")
        print(f"{m:>3} | " + " | ".join(cells) + tag)

    print("\n=== B. High-l SR ERI (int2e_sph) across small-to-medium R (demonstrating quadmath failure + binary dependence) ===")
    libs = {name: load_lib(path) for name, path in BUILDS.items()}
    cases = [
        ("ctrl_n20_m39", (10, 10, 10, 9), 0.5, [0.2, 0.5, 0.6, 0.8, 1.0, 2.0]),
        ("oob_n21_m41",  (10, 10, 10, 10), 0.5, [0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.5, 2.0]),
        ("oob_n22_m43",  (11, 11, 10, 10), 0.5, [0.2, 0.5, 0.6, 0.8, 1.0, 1.5]),
        ("oob_n24_m47",  (12, 12, 11, 11), 0.5, [0.2, 0.5, 0.6, 0.8, 1.0, 1.5]),
    ]
    eri_out = {}
    for cname, ls, om, Rs in cases:
        nroots = sum(ls) // 2 + 1
        m_idx = 2 * nroots - 1
        print(f"\n--- {cname}: l={ls} (nroots={nroots}, m={m_idx}, omega=-{om}) ---")
        print(f"{'R':>4} {'x':>6} | " + " | ".join(f"{k:>20}" for k in BUILDS))
        c_dict = {}
        for R in Rs:
            a0 = (2.1 * 1.9) / 4.0
            x = a0 * R * R
            ref_full = eval_int2e_sph(libs["gcc16_O2_quad_v621"], ls, R, 0.0)
            ref_lr = eval_int2e_sph(libs["gcc16_O2_quad_v621"], ls, R, +om)
            ref_sr = ref_full - ref_lr
            ref_norm = float(np.linalg.norm(ref_sr))
            row = {}
            cells = []
            for name, lib in libs.items():
                sr = eval_int2e_sph(lib, ls, R, -om)
                if np.any(~np.isfinite(sr)):
                    row[name] = {"rel_err": None, "nonfinite": True}
                    cells.append(f"{'NaN/Inf':>20}")
                else:
                    rel = float(np.linalg.norm(sr - ref_sr)) / max(ref_norm, 1e-300)
                    row[name] = {"rel_err": rel, "nonfinite": False}
                    cells.append(f"{rel:>20.3e}")
            c_dict[str(R)] = {"x": x, "ref_norm": ref_norm, "builds": row}
            print(f"{R:>4.1f} {x:>6.3f} | " + " | ".join(cells))
        eri_out[cname] = c_dict

    out = ROOT / "outputs/phase1_deep_check.json"
    out.write_text(json.dumps({"turnover_memory_dump": dumps, "eri_small_R": eri_out}, indent=2))
    print(f"\nSaved to {out}")


if __name__ == "__main__":
    main()
