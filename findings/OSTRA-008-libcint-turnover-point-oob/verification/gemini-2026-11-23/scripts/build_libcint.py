#!/usr/bin/env python3
"""Build untouched libcint v6.2.1 (or a patched copy) into an isolated build dir."""
import os
import sys
import subprocess
import concurrent.futures
from pathlib import Path

ROOT = Path("/Users/napolesd/Documents/Antigravity/PhysAud/verif-OSTRA-008")
DEFAULT_SRC = ROOT / "libcint"

CORE_SRCS = [
    "src/c2f.c", "src/cart2sph.c", "src/cint1e.c", "src/cint2e.c", "src/cint_bas.c",
    "src/fblas.c", "src/g1e.c", "src/g2e.c", "src/misc.c", "src/optimizer.c",
    "src/fmt.c", "src/rys_wheeler.c", "src/eigh.c", "src/rys_roots.c", "src/find_roots.c",
    "src/cint2c2e.c", "src/g2c2e.c", "src/cint3c2e.c", "src/g3c2e.c",
    "src/cint3c1e.c", "src/g3c1e.c", "src/breit.c",
    "src/cint1e_a.c", "src/cint3c1e_a.c",
    "src/cint1e_grids.c", "src/g1e_grids.c",
    "src/autocode/breit1.c", "src/autocode/dkb.c", "src/autocode/gaunt1.c",
    "src/autocode/grad1.c", "src/autocode/grad2.c", "src/autocode/hess.c",
    "src/autocode/int3c1e.c", "src/autocode/int3c2e.c",
    "src/autocode/intor1.c", "src/autocode/intor2.c", "src/autocode/intor3.c", "src/autocode/intor4.c",
    "src/autocode/deriv3.c", "src/autocode/int1e_grids1.c",
    "src/autocode/deriv4.c", "src/autocode/lresc.c",
]


def configure_headers(src_dir: Path, build_dir: Path, have_quadmath: bool = False):
    (build_dir / "include").mkdir(parents=True, exist_ok=True)
    (build_dir / "src").mkdir(parents=True, exist_ok=True)

    cint_h_in = (src_dir / "include" / "cint.h.in").read_text()
    cint_h = (
        cint_h_in
        .replace("@cint_VERSION@", "6.2.1")
        .replace("@cint_SOVERSION@", "6")
        .replace("#cmakedefine I8", "/* #undef I8 */")
        .replace("#cmakedefine CACHE_SIZE_I8", "#define CACHE_SIZE_I8")
    )
    (build_dir / "include" / "cint.h").write_text(cint_h)

    cfg_in = (src_dir / "src" / "cint_config.h.in").read_text()
    cfg = (
        cfg_in
        .replace("#cmakedefine HAVE_EXPL", "#define HAVE_EXPL")
        .replace("#cmakedefine HAVE_SQRTL", "#define HAVE_SQRTL")
        .replace("#cmakedefine HAVE_FABSL", "#define HAVE_FABSL")
        .replace(
            "#cmakedefine HAVE_QUADMATH_H",
            "#define HAVE_QUADMATH_H" if have_quadmath else "/* #undef HAVE_QUADMATH_H */",
        )
        .replace("#cmakedefine WITH_RANGE_COULOMB", "/* #undef WITH_RANGE_COULOMB */")
    )
    (build_dir / "src" / "cint_config.h").write_text(cfg)


def build_target(
    build_dir: Path,
    cc: str,
    cflags: list[str],
    ldflags: list[str],
    src_dir: Path = DEFAULT_SRC,
    fmt_override: Path | None = None,
    have_quadmath: bool = False,
    libname: str = "libcint.dylib",
):
    configure_headers(src_dir, build_dir, have_quadmath=have_quadmath)
    obj_dir = build_dir / "obj"
    obj_dir.mkdir(parents=True, exist_ok=True)

    inc_flags = [
        f"-I{build_dir / 'include'}",
        f"-I{src_dir / 'include'}",
        f"-I{build_dir / 'src'}",
        f"-I{src_dir / 'src'}",
    ]

    jobs = []
    objs = []
    for rel in CORE_SRCS:
        src_path = fmt_override if (fmt_override is not None and rel == "src/fmt.c") else (src_dir / rel)
        obj_path = obj_dir / (rel.replace("/", "_") + ".o")
        objs.append(str(obj_path))
        cmd = [cc, "-std=c99", "-fPIC", "-DWITH_CINT2_INTERFACE", *cflags, *inc_flags, "-c", str(src_path), "-o", str(obj_path)]
        jobs.append((rel, cmd))

    def run_one(item):
        rel, cmd = item
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"Failed compiling {rel}:\n{res.stderr}")
        return rel

    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
        list(ex.map(run_one, jobs))

    out_lib = build_dir / libname
    link_cmd = [cc, "-shared", "-o", str(out_lib), *objs, "-lm", *ldflags]
    res = subprocess.run(link_cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"Failed linking {out_lib}:\n{res.stderr}")
    print(f"Built {out_lib} ({len(objs)} objects, cc={cc}, quadmath={have_quadmath})")
    return out_lib


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    if mode in ("asan", "all"):
        build_target(
            ROOT / "build-asan",
            cc="clang",
            cflags=["-O1", "-g", "-fsanitize=address", "-fno-omit-frame-pointer"],
            ldflags=["-fsanitize=address"],
            libname="libcint_asan.dylib",
        )
    if mode in ("o2", "all"):
        build_target(
            ROOT / "build-o2",
            cc="clang",
            cflags=["-O2", "-g"],
            ldflags=[],
            libname="libcint.dylib",
        )
    if mode in ("o0", "all"):
        build_target(
            ROOT / "build-o0",
            cc="clang",
            cflags=["-O0", "-g"],
            ldflags=[],
            libname="libcint.dylib",
        )
    if mode in ("gcc", "all"):
        # Check if gcc-16 has quadmath.h on arm64 (usually x86_64 only, or maybe available)
        test_q = subprocess.run(
            ["/opt/homebrew/bin/gcc-16", "-x", "c", "-", "-lquadmath", "-o", "/tmp/test_q_cint"],
            input="#include <quadmath.h>\nint main(){ __float128 x = 1.0q; return (int)sqrtq(x) - 1; }\n",
            capture_output=True, text=True,
        )
        has_q = (test_q.returncode == 0)
        print(f"gcc-16 libquadmath available: {has_q}")
        build_target(
            ROOT / "build-gcc",
            cc="/opt/homebrew/bin/gcc-16",
            cflags=["-O2", "-g"],
            ldflags=["-lquadmath"] if has_q else [],
            have_quadmath=has_q,
            libname="libcint.dylib",
        )
