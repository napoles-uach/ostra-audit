"""Oracle: spin-polarized evaluation at rho_a = rho_b must reproduce the unpolarized one.

For each functional, evaluate on a grid (typical + edge region) with nspin=1 and with
nspin=2 at equal spin densities, and compare
  zk             : zk_pol == zk_unp
  vrho           : vrho_a == vrho_b == vrho_unp
  vsigma         : (vs_aa + vs_ab + vs_bb)/4 == vs_unp ;  vs_aa == vs_bb
  vtau, vlapl    : (v_a + v_b)/2 == v_unp ;  v_a == v_b
Usage: python3 spinprobe.py ID [ID ...]   -> one JSON line per ID on stdout
"""
import sys, json, ctypes, warnings
import numpy as np
import lxc

warnings.filterwarnings("ignore")
C_F = 0.3 * (3 * np.pi**2) ** (2 / 3)


def grid():
    rho = np.logspace(-6, 4, 41)                     # total density
    s = np.array([0, 1e-6, 1e-3, 0.05, 0.3, 1, 2, 4, 10, 100])
    alpha = np.array([0, 1e-6, 0.3, 1, 2, 10])
    R, S, A = [a.ravel() for a in np.meshgrid(rho, s, alpha, indexing="ij")]
    kf = (3 * np.pi**2 * R) ** (1 / 3)
    sigma = (2 * kf * R * S) ** 2
    tauW = sigma / (8 * R)
    tau = tauW + A * C_F * R ** (5 / 3)
    q = np.tile(np.array([-2.0, 0.0, 1.5]), len(R) // 3 + 1)[: len(R)]
    lapl = q * R ** (5 / 3)
    return R, S, A, sigma, tau, lapl


def call(f, x):
    """exc (+vxc if available) using lxc.Func, falling back to exc-only."""
    if f.flags & lxc.FLAG_VXC:
        return f.exc_vxc(x)
    rho = np.ascontiguousarray(x["rho"]); n = rho.shape[0]
    zk = np.zeros(n)
    P = ctypes.c_void_p(f.p); N = ctypes.c_size_t(n); p = lxc._p
    if f.family == "LDA":
        lxc.L.xc_lda_exc(P, N, p(rho), p(zk))
    elif f.family == "GGA":
        lxc.L.xc_gga_exc(P, N, p(rho), p(np.ascontiguousarray(x["sigma"])), p(zk))
    else:
        lxc.L.xc_mgga_exc(P, N, p(rho), p(np.ascontiguousarray(x["sigma"])),
                          p(np.ascontiguousarray(x["lapl"])), p(np.ascontiguousarray(x["tau"])), p(zk))
    return {"zk": zk}


def rel(a, b, scale):
    with np.errstate(all="ignore"):
        d = np.abs(a - b) / (np.abs(scale) + 1e-300)
    return d


def lower_thresholds(f):
    D = ctypes.c_double
    for fn, v in (("xc_func_set_dens_threshold", 1e-25), ("xc_func_set_sigma_threshold", 1e-40),
                  ("xc_func_set_tau_threshold", 1e-30)):
        g = getattr(lxc.L, fn); g.argtypes = [ctypes.c_void_p, D]; g(ctypes.c_void_p(f.p), D(v))


def probe(fid, low=False):
    out = {"id": fid, "name": lxc.name(fid), "low_thr": low}
    try:
        fu = lxc.Func(fid, 1)
    except Exception:
        out["status"] = "init_unp_failed"; return out
    try:
        fp = lxc.Func(fid, 2)
    except Exception:
        fu.close(); out["status"] = "init_pol_failed"; return out
    out["family"] = fu.family
    if not (fu.flags & lxc.FLAG_EXC):
        fu.close(); fp.close(); out["status"] = "no_exc"; return out
    if low:
        lower_thresholds(fu); lower_thresholds(fp)
    R, S, A, sigma, tau, lapl = grid()
    n = len(R)
    xu = {"rho": R[:, None], "sigma": sigma[:, None], "tau": tau[:, None], "lapl": lapl[:, None]}
    h = 0.5
    xp = {"rho": np.stack([h * R, h * R], 1),
          "sigma": np.stack([sigma / 4, sigma / 4, sigma / 4], 1),
          "tau": np.stack([h * tau, h * tau], 1),
          "lapl": np.stack([h * lapl, h * lapl], 1)}
    ou, op = call(fu, xu), call(fp, xp)
    fu.close(); fp.close()

    checks = {}
    # reference scale per point: magnitude of the unpolarized quantity, or of e for energy
    def add(key, a, b, scale):
        d = rel(a, b, scale)
        nan_mis = np.isfinite(a) != np.isfinite(b)
        fin = np.isfinite(a) & np.isfinite(b)
        dd = np.where(fin, d, 0.0)
        i = int(np.argmax(dd))
        checks[key] = {"max": float(dd[i]), "n_gt_1e-6": int((dd > 1e-6).sum()),
                       "n_nan_mismatch": int(nan_mis.sum()),
                       "at": {"rho": float(R[i]), "s": float(S[i]), "alpha": float(A[i])},
                       "vals": [float(a[i]), float(b[i])]}
        if nan_mis.any():
            j = int(np.argmax(nan_mis))
            checks[key]["nan_at"] = {"rho": float(R[j]), "s": float(S[j]), "alpha": float(A[j]),
                                     "vals": [float(a[j]), float(b[j])]}

    add("zk", op["zk"], ou["zk"], ou["zk"])
    if "vrho" in ou:
        vr = ou["vrho"][:, 0]
        add("vrho", op["vrho"][:, 0], vr, vr)
        add("vrho_ab_sym", op["vrho"][:, 0], op["vrho"][:, 1], vr)
    if "vsigma" in ou:
        vs = ou["vsigma"][:, 0]
        # natural scale for vsigma: compare sigma*vsigma against e = rho*zk
        e = R * ou["zk"]
        sc = np.maximum(np.abs(vs) * sigma, np.abs(e)) / np.maximum(sigma, 1e-300)
        add("vsigma", op["vsigma"].sum(1) / 4, vs, sc)
        add("vsigma_ab_sym", op["vsigma"][:, 0], op["vsigma"][:, 2], sc)
    for k in ("vtau", "vlapl"):
        if k in ou and np.any(ou[k] != 0):
            v = ou[k][:, 0]
            q = tau if k == "vtau" else np.abs(lapl)
            e = R * ou["zk"]
            sc = np.maximum(np.abs(v) * q, np.abs(e)) / np.maximum(q, 1e-300)
            add(k, op[k].mean(1), v, sc)
            add(k + "_ab_sym", op[k][:, 0], op[k][:, 1], sc)
    out["checks"] = checks
    out["worst"] = max(c["max"] for c in checks.values())
    out["nan_mismatch"] = sum(c["n_nan_mismatch"] for c in checks.values())
    out["status"] = "ok"
    return out


if __name__ == "__main__":
    for a in sys.argv[1:]:
        r = probe(int(a))
        if r.get("status") == "ok":
            r["low"] = probe(int(a), low=True)
        print(json.dumps(r), flush=True)
