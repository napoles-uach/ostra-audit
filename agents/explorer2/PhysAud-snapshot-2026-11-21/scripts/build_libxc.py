#!/usr/bin/env python3
"""Build libxc.dylib from libxc-ref matching CMakeLists.txt flags when cmake is absent."""
import concurrent.futures
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC_ROOT = ROOT / "libxc-ref"
BUILD_DIR = ROOT / "build-libxc"
OBJ_DIR = BUILD_DIR / "obj"
INC_SRC_DIR = BUILD_DIR / "src"


def main():
    OBJ_DIR.mkdir(parents=True, exist_ok=True)
    INC_SRC_DIR.mkdir(parents=True, exist_ok=True)

    conf_ac = (SRC_ROOT / "configure.ac").read_text()
    m = re.search(r"AC_INIT\(\[libxc\],\[(\d+)\.(\d+)\.(\d+)\]", conf_ac)
    maj, min_, mic = m.groups() if m else ("7", "1", "2")
    ver = f"{maj}.{min_}.{mic}"

    xc_ver = (SRC_ROOT / "xc_version.h.in").read_text()
    xc_ver = (
        xc_ver.replace("@VERSION@", ver)
        .replace("@XC_MAJOR_VERSION@", maj)
        .replace("@XC_MINOR_VERSION@", min_)
        .replace("@XC_MICRO_VERSION@", mic)
    )
    (BUILD_DIR / "xc_version.h").write_text(xc_ver)

    cfg = f'/* Generated config.h */\n#define PACKAGE_VERSION "{ver}"\n#include <stdio.h>\n'
    (BUILD_DIR / "config.h").write_text(cfg)

    xc_h = (SRC_ROOT / "src" / "xc.h").read_text().replace("<xc_version.h>", '"xc_version.h"')
    (INC_SRC_DIR / "xc.h").write_text(xc_h)

    sources_mk = (SRC_ROOT / "src" / "sources.mk").read_text()
    c_files = sorted(set(re.findall(r"[a-z_0-9][a-z_0-9]*\.c", sources_mk)))
    print(f"Found {len(c_files)} C source files to compile (version {ver})...", flush=True)

    cflags = [
        "/usr/bin/clang",
        "-std=c99",
        "-O2",
        "-fPIC",
        "-ffp-contract=off",
        "-DHAVE_CBRT",
        "-DXC_MAXORDER=2",
        "-DXC_ENFORCE_FERMI_HOLE_CURVATURE",
        f"-I{BUILD_DIR}",
        f"-I{INC_SRC_DIR}",
        f"-I{SRC_ROOT / 'src'}",
    ]

    objs = []

    def compile_one(cf):
        src_path = SRC_ROOT / "src" / cf
        obj_path = OBJ_DIR / (cf[:-2] + ".o")
        if obj_path.exists() and obj_path.stat().st_mtime >= src_path.stat().st_mtime:
            return str(obj_path)
        cmd = cflags + ["-c", str(src_path), "-o", str(obj_path)]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"Failed compiling {cf}:\n{res.stderr}")
        return str(obj_path)

    workers = os.cpu_count() or 8
    done = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(compile_one, cf): cf for cf in c_files}
        for fut in concurrent.futures.as_completed(futs):
            objs.append(fut.result())
            done += 1
            if done % 100 == 0 or done == len(c_files):
                print(f"Compiled {done}/{len(c_files)}", flush=True)

    dylib = BUILD_DIR / "libxc.dylib"
    obj_list_file = BUILD_DIR / "objects.txt"
    obj_list_file.write_text("\n".join(sorted(objs)) + "\n")
    link_cmd = [
        "/usr/bin/clang",
        "-dynamiclib",
        f"-filelist",
        str(obj_list_file),
        "-o",
        str(dylib),
        "-lm",
    ]
    res = subprocess.run(link_cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"Link failed:\n{res.stderr}")
    print(f"Successfully built {dylib} ({dylib.stat().st_size} bytes)", flush=True)


if __name__ == "__main__":
    main()
