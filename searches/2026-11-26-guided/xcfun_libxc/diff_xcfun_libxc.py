# OSTRA guided search, 2026-11-26: XCFun vs libxc in PySCF, same xc string or explicit component pairs.
# Realistic densities: grids of H2O (RKS), OH radical (UKS) and stretched N2 (UKS triplet-like, low density tails).
# Metrics: integrated E_xc difference (mEh) and max relative difference of v, f, k (deriv 1-3) over points with rho > 1e-8.
import json, sys, traceback, numpy as np
from pyscf import gto, dft
from pyscf.dft import xcfun, libxc
PAIRS_EXTRA = {  # XCFun component name -> libxc name, where the same string does not parse in libxc
 "SLATERX": "LDA_X", "VWN5C": "LDA_C_VWN", "VWN3C": "LDA_C_VWN_3", "PZ81C": "LDA_C_PZ", "PW92C": "LDA_C_PW",
 "BECKEX": "GGA_X_B88", "PBEX": "GGA_X_PBE", "PBEC": "GGA_C_PBE", "RPBEX": "GGA_X_RPBE", "REVPBEX": "GGA_X_PBE_R",
 "PBESOLX": "GGA_X_PBE_SOL", "PBEINTX": "GGA_X_PBEINT", "PBEINTC": "GGA_C_PBEINT", "PW86X": "GGA_X_PW86",
 "PW91X": "GGA_X_PW91", "PW91C": "GGA_C_PW91", "OPTX": "GGA_X_OPTX", "LYPC": "GGA_C_LYP", "P86C": "GGA_C_P86",
 "APBEX": "GGA_X_APBE", "APBEC": "GGA_C_APBE", "SPBEC": "GGA_C_SPBE", "PBELOCC": "GGA_C_PBELOC",
 "ZVPBEINTC": "GGA_C_ZVPBEINT", "ZVPBESOLC": "GGA_C_ZVPBESOL", "KT2XC": "GGA_XC_KT2", "KT3XC": "GGA_XC_KT3",
 "TPSSX": "MGGA_X_TPSS", "TPSSC": "MGGA_C_TPSS", "TPSSLOCC": "MGGA_C_TPSSLOC", "REVTPSSX": "MGGA_X_REVTPSS",
 "REVTPSSC": "MGGA_C_REVTPSS", "SCANX": "MGGA_X_SCAN", "SCANC": "MGGA_C_SCAN", "RSCANX": "MGGA_X_RSCAN",
 "RSCANC": "MGGA_C_RSCAN", "R2SCANX": "MGGA_X_R2SCAN", "R2SCANC": "MGGA_C_R2SCAN", "R4SCANX": "MGGA_X_R4SCAN",
 "R4SCANC": "MGGA_C_R2SCAN", "RPPSCANX": "MGGA_X_RPPSCAN", "RPPSCANC": "MGGA_C_RPPSCAN",
 "M05X": "HYB_MGGA_X_M05", "M05C": "MGGA_C_M05", "M052XX": "HYB_MGGA_X_M05_2X", "M052XC": "MGGA_C_M05_2X",
 "M06X": "HYB_MGGA_X_M06", "M06C": "MGGA_C_M06", "M062XX": "HYB_MGGA_X_M06_2X", "M062XC": "MGGA_C_M06_2X",
 "M06LX": "MGGA_X_M06_L", "M06LC": "MGGA_C_M06_L", "M06HFX": "HYB_MGGA_X_M06_HF", "M06HFC": "MGGA_C_M06_HF",
 "CSC": "MGGA_C_CS", "BRX": "MGGA_X_BR89", "TF": "LDA_K_TF", "TFK": "LDA_K_TF", "TW": "GGA_K_VW",
 "PW91K": "GGA_K_PW91",
}
MOLS = {"H2O": ("O 0 0 0.117; H 0 0.757 -0.469; H 0 -0.757 -0.469", 0),
        "OH": ("O 0 0 0; H 0 0 0.970", 1),
        "N2_stretched": ("N 0 0 0; N 0 0 2.2", 2)}
def densities():
    out = {}
    for name, (atom, spin) in MOLS.items():
        mol = gto.M(atom=atom, spin=spin, basis="def2-svp", verbose=0)
        mf = (dft.UKS(mol) if spin else dft.RKS(mol)); mf.xc = "pbe"; mf.grids.level = 3; mf.kernel()
        ni = dft.numint.NumInt(); dm = mf.make_rdm1(); g = mf.grids
        ao = ni.eval_ao(mol, g.coords, deriv=2)
        if spin:
            rho = np.array([ni.eval_rho(mol, ao, d, xctype="MGGA", with_lapl=True) for d in dm])
        else:
            rho = ni.eval_rho(mol, ao, dm, xctype="MGGA", with_lapl=True)
        out[name] = (1 if spin else 0, rho, g.weights)
    return out
def ev(mod, xc, rho, spin, deriv):
    return mod.eval_xc(xc, rho, spin=spin, deriv=deriv)
def rel(a, b, mask):
    if isinstance(a, (tuple, list)):
        vals = []
        for ai, bi in zip(a, b):
            if ai is None or bi is None:
                if (ai is None) != (bi is None) and not (ai is None and np.allclose(bi, 0)) and not (bi is None and np.allclose(ai, 0)):
                    vals.append("component missing")
                continue
            vals.append(rel(ai, bi, mask))
        bad = [v for v in vals if isinstance(v, str)]
        return bad[0] if bad else (max(vals) if vals else 0.0)
    a = np.asarray(a); b = np.asarray(b)
    if a.shape != b.shape: return f"shape {a.shape} vs {b.shape}"
    a = a[mask] if a.shape[0] == mask.size else a[..., mask]; b = b[mask] if b.shape[0] == mask.size else b[..., mask]
    if not np.all(np.isfinite(a)): return "nonfinite_xcfun"
    if not np.all(np.isfinite(b)): return "nonfinite_libxc"
    s = np.abs(b).max()
    return float(np.abs(a - b).max() / s) if s > 0 else float(np.abs(a - b).max())
def compare(xf, lx, D):
    res = {"xcfun": xf, "libxc": lx}
    try:
        nx = dft.numint.NumInt(); nx.libxc = xcfun
        nl = dft.numint.NumInt(); nl.libxc = libxc
        res["rsh_hyb"] = [list(map(float, nx.rsh_and_hybrid_coeff(xf))), list(map(float, nl.rsh_and_hybrid_coeff(lx)))]
    except Exception as e:
        res["hyb_err"] = f"{type(e).__name__}: {str(e)[:80]}"
    for mname, (spin, rho, w) in D.items():
        r = rho if (xcfun.xc_type(xf) == "MGGA" or libxc.xc_type(lx) == "MGGA") else rho[..., :4, :]
        typ = max(xcfun.xc_type(xf), libxc.xc_type(lx), key=["HF", "LDA", "GGA", "MGGA"].index)
        r = rho[..., :1, :] if typ == "LDA" else (rho[..., :4, :] if typ == "GGA" else rho)
        if typ == "LDA": r = r[..., 0, :]
        dens = (rho[0, 0] + rho[1, 0]) if spin else rho[0]
        mask = dens > 1e-6
        m = {}
        for deriv in (1, 2, 3):
            try:
                X = ev(xcfun, xf, r, spin, deriv)
            except Exception as e:
                m[f"d{deriv}"] = f"xcfun error: {type(e).__name__}: {str(e)[:200]}"; break
            try:
                L = ev(libxc, lx, r, spin, deriv)
            except Exception as e:
                m[f"d{deriv}"] = f"libxc error: {type(e).__name__}: {str(e)[:200]}"; break
            if deriv == 1:
                ex, el = np.asarray(X[0]), np.asarray(L[0])
                m["dE_mEh"] = float(np.sum(w * dens * (ex - el)) * 1000) if np.all(np.isfinite(ex)) and np.all(np.isfinite(el)) else "nonfinite"
                m["E_mEh"] = float(np.sum(w * dens * el) * 1000)
                m["exc"] = rel(ex, el, mask)
            m[f"d{deriv}"] = rel(X[deriv], L[deriv], mask) if X[deriv] is not None and L[deriv] is not None else "none"
        res[mname] = m
    return res
if __name__ == "__main__":
    D = densities()
    names = sorted(set(xcfun.XC_CODES) | set(getattr(xcfun, "XC_ALIAS", {})))
    out = []
    for xf in names:
        if xf in ("VV10",) or "VV10" in xf or xf.endswith("MV") or xf.endswith("M_V") or xf.endswith("XV") or xf.endswith("X_V"):
            continue  # NLC: see OSTRA-017
        lx = PAIRS_EXTRA.get(xf, xf)
        try:
            libxc.parse_xc(lx); xcfun.parse_xc(xf)
        except Exception as e:
            out.append({"xcfun": xf, "libxc": lx, "skip": f"{type(e).__name__}: {str(e)[:80]}"}); continue
        try:
            out.append(compare(xf, lx, D))
        except Exception as e:
            out.append({"xcfun": xf, "libxc": lx, "error": traceback.format_exc()[-300:]})
        print(json.dumps(out[-1])[:200], flush=True)
    json.dump(out, open("diff_results.json", "w"), indent=1)
