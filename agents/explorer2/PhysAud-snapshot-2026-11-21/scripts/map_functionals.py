#!/usr/bin/env python3
"""Map libxc source files (maple/*.mpl and src/*.c) to all 727 registered functionals."""
import ctypes
import json
import os
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC_ROOT = ROOT / "libxc-ref"
DYLIB = ROOT / "build-libxc" / "libxc.dylib"


def load_libxc_registry():
    L = ctypes.CDLL(str(DYLIB))
    L.xc_number_of_functionals.restype = ctypes.c_int
    L.xc_available_functional_numbers.argtypes = [ctypes.POINTER(ctypes.c_int)]
    L.xc_functional_get_name.argtypes = [ctypes.c_int]
    L.xc_functional_get_name.restype = ctypes.c_char_p
    L.xc_func_alloc.restype = ctypes.c_void_p
    L.xc_func_init.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
    L.xc_func_init.restype = ctypes.c_int
    L.xc_func_end.argtypes = [ctypes.c_void_p]
    L.xc_func_free.argtypes = [ctypes.c_void_p]
    L.xc_func_get_info.argtypes = [ctypes.c_void_p]
    L.xc_func_get_info.restype = ctypes.c_void_p
    L.xc_func_info_get_family.argtypes = [ctypes.c_void_p]
    L.xc_func_info_get_kind.argtypes = [ctypes.c_void_p]
    L.xc_func_info_get_flags.argtypes = [ctypes.c_void_p]
    L.xc_func_info_get_references.argtypes = [ctypes.c_void_p, ctypes.c_int]
    L.xc_func_info_get_references.restype = ctypes.c_void_p
    L.xc_func_reference_get_ref.argtypes = [ctypes.c_void_p]
    L.xc_func_reference_get_ref.restype = ctypes.c_char_p
    L.xc_func_reference_get_doi.argtypes = [ctypes.c_void_p]
    L.xc_func_reference_get_doi.restype = ctypes.c_char_p
    L.xc_num_aux_funcs.argtypes = [ctypes.c_void_p]
    L.xc_num_aux_funcs.restype = ctypes.c_int
    L.xc_aux_func_ids.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]

    n = L.xc_number_of_functionals()
    arr = (ctypes.c_int * n)()
    L.xc_available_functional_numbers(arr)
    fids = list(arr)

    funcs = {}
    for fid in fids:
        fname = L.xc_functional_get_name(fid).decode()
        p = L.xc_func_alloc()
        rc = L.xc_func_init(p, fid, 1)
        if rc != 0:
            rc = L.xc_func_init(p, fid, 2)
        if rc != 0:
            L.xc_func_free(p)
            continue
        info = L.xc_func_get_info(p)
        fam = L.xc_func_info_get_family(info)
        kind = L.xc_func_info_get_kind(info)
        flags = L.xc_func_info_get_flags(info)
        refs = []
        for r_idx in range(5):
            rptr = L.xc_func_info_get_references(info, r_idx)
            if not rptr:
                break
            rstr = L.xc_func_reference_get_ref(rptr)
            rdoi = L.xc_func_reference_get_doi(rptr)
            refs.append({
                "ref": rstr.decode("utf-8", errors="replace") if rstr else "",
                "doi": rdoi.decode("utf-8", errors="replace") if rdoi else "",
            })
        n_aux = L.xc_num_aux_funcs(p)
        aux_ids = []
        if n_aux > 0:
            aux_arr = (ctypes.c_int * n_aux)()
            L.xc_aux_func_ids(p, aux_arr)
            aux_ids = list(aux_arr)
        L.xc_func_end(p)
        L.xc_func_free(p)
        funcs[fid] = {
            "id": fid,
            "name": fname,
            "family": fam,
            "kind": kind,
            "flags": flags,
            "aux_ids": aux_ids,
            "refs": refs,
        }
    return funcs


def build_source_graph(funcs):
    # 1. Resolve Maple $include graph
    mpl_files = list(SRC_ROOT.glob("maple/**/*.mpl"))
    mpl_by_name = {}
    for mp in mpl_files:
        rel = str(mp.relative_to(SRC_ROOT))
        mpl_by_name[mp.name] = rel
        mpl_by_name[rel] = rel

    direct_includes = {}  # mpl_rel -> set of included mpl_rel
    for mp in mpl_files:
        rel = str(mp.relative_to(SRC_ROOT))
        txt = mp.read_text()
        incs = set()
        # Every .mpl implicitly has maple/util.mpl in scope via maple2c!
        if rel != "maple/util.mpl":
            incs.add("maple/util.mpl")
        for m in re.finditer(r'\$include\s+"([^"]+)"', txt):
            inc_name = m.group(1)
            if inc_name in mpl_by_name:
                incs.add(mpl_by_name[inc_name])
            else:
                # check same subdir
                cand = str((mp.parent / inc_name).relative_to(SRC_ROOT))
                if cand in mpl_by_name:
                    incs.add(cand)
        direct_includes[rel] = incs

    # Transitive closure: which leaf .mpl files include target_mpl?
    def all_ancestors(rel, seen=None):
        if seen is None:
            seen = set()
        for inc in direct_includes.get(rel, ()):
            if inc not in seen:
                seen.add(inc)
                all_ancestors(inc, seen)
        return seen

    mpl_ancestors = {rel: all_ancestors(rel) | {rel} for rel in direct_includes}

    # 2. Map src/*.c to functional IDs and to maple2c/*.c -> maple/*.mpl
    name_to_id = {v["name"].lower(): k for k, v in funcs.items()}
    c_to_fids = {}
    fid_to_c = {}
    c_to_mpl = {}

    for cf in sorted((SRC_ROOT / "src").glob("*.c")):
        if cf.name.startswith("funcs_"):
            continue
        rel_c = f"src/{cf.name}"
        txt = cf.read_text()
        fids = set()
        for m in re.finditer(r"xc_func_info_(\w+)\s*=", txt):
            fn = m.group(1).lower()
            if fn in name_to_id:
                fids.add(name_to_id[fn])
                fid_to_c[name_to_id[fn]] = rel_c
        c_to_fids[rel_c] = fids

        mpls = set()
        for m in re.finditer(r'#include\s+"maple2c/([^"]+)\.c"', txt):
            sub = m.group(1)  # e.g. gga_exc/gga_x_pbe
            cand = f"maple/{sub}.mpl"
            if cand in mpl_ancestors:
                mpls.add(cand)
        c_to_mpl[rel_c] = mpls

    # 3. Expand mixes transitively: if functional M mixes A, then any source affecting A affects M!
    def transitive_dependents(base_fids):
        out = set(base_fids)
        changed = True
        while changed:
            changed = False
            for fid, info in funcs.items():
                if fid not in out and any(a in out for a in info["aux_ids"]):
                    out.add(fid)
                    changed = True
        return out

    mpl_to_fids = {}
    for target_mpl in mpl_ancestors:
        base = set()
        for rel_c, mpls in c_to_mpl.items():
            for leaf_mpl in mpls:
                if target_mpl in mpl_ancestors.get(leaf_mpl, set()):
                    base |= c_to_fids.get(rel_c, set())
        mpl_to_fids[target_mpl] = sorted(transitive_dependents(base))

    c_to_all_fids = {
        rel_c: sorted(transitive_dependents(fids))
        for rel_c, fids in c_to_fids.items()
    }

    return {
        "funcs": funcs,
        "fid_to_c": fid_to_c,
        "mpl_to_fids": mpl_to_fids,
        "c_to_all_fids": c_to_all_fids,
    }


if __name__ == "__main__":
    funcs = load_libxc_registry()
    graph = build_source_graph(funcs)
    out_file = ROOT / "scripts" / "source_graph.json"
    out_file.write_text(json.dumps(graph, indent=2))
    print(f"Mapped {len(funcs)} functionals. Saved {out_file}.")
    # Print priority functionals check
    prio_targets = [
        "gga_x_pbe", "gga_c_pbe", "hyb_gga_xc_pbeh", "hyb_gga_xc_hse06",
        "mgga_x_scan", "mgga_c_scan", "mgga_x_r2scan", "mgga_c_r2scan",
        "mgga_x_tpss", "mgga_c_tpss", "hyb_gga_xc_b3lyp", "hyb_gga_xc_wb97x"
    ]
    by_name = {v["name"]: v for v in funcs.values()}
    for pt in prio_targets:
        v = by_name.get(pt)
        print(f"  {pt:22s} -> id={v['id']:3d}, aux={v['aux_ids']}, c={graph['fid_to_c'].get(v['id'])}")
