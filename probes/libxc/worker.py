"""Probe one libxc functional: (a) vxc vs finite differences of e=rho*zk, (b) NaN/Inf at edge inputs."""
import sys, json, warnings
import numpy as np
import lxc

warnings.filterwarnings("ignore")
C_F = 0.3 * (3 * np.pi**2) ** (2 / 3)
NPT = 300
TOL = 1e-4
KEYS = {"LDA": ["rho"], "GGA": ["rho", "sigma"], "MGGA": ["rho", "sigma", "lapl", "tau"]}


def sample(rng, ns, n, edge=False):
    if not edge:
        rho_t = 10 ** rng.uniform(-5, 2, n)
        zeta = rng.uniform(-0.8, 0.8, n) if ns == 2 else np.zeros(n)
        s = rng.uniform(0.0, 4.0, (n, ns))
        alpha = rng.uniform(0.2, 3.0, (n, ns))
        cosg = rng.uniform(-0.8, 0.8, n)
    else:
        rho_t = np.concatenate([10 ** rng.uniform(-16, -6, n - 8), [0, 1e-30, 1e-20, 1e-12, 1e3, 1e5, 1e-10, 1e-10]])
        zeta = rng.choice([-1, -0.999999, 0, 0.999999, 1], n) if ns == 2 else np.zeros(n)
        s = 10 ** rng.uniform(-6, 4, (n, ns)); s[: n // 5] = 0
        alpha = rng.choice([0, 1e-10, 1, 1e4], (n, ns))
        cosg = rng.choice([-1, 0, 1], n)
    if ns == 1:
        rho = rho_t[:, None]
    else:
        rho = np.stack([rho_t * (1 + zeta) / 2, rho_t * (1 - zeta) / 2], 1)
    # per-spin-channel quantities (for ns==1 the channel is the total density)
    rr = rho if ns == 2 else rho  # channel densities
    kf = (3 * np.pi**2 * (2 * rr if ns == 2 else rr)) ** (1 / 3)
    grad = 2 * kf * (2 * rr if ns == 2 else rr) * s / (2 if ns == 2 else 1)
    sig_ch = grad**2
    if ns == 1:
        sigma = sig_ch
    else:
        sigma = np.stack([sig_ch[:, 0], cosg * np.sqrt(sig_ch[:, 0] * sig_ch[:, 1]), sig_ch[:, 1]], 1)
    with np.errstate(all="ignore"):
        tauW = np.where(rr > 0, sig_ch / (8 * rr), 0)
        tauU = C_F * (2 ** (2 / 3) if ns == 2 else 1) * rr ** (5 / 3)
    tau = tauW + alpha * tauU
    lapl = rng.uniform(-3, 3, (n, ns)) * rr ** (5 / 3)
    return {"rho": rho, "sigma": sigma, "tau": tau, "lapl": lapl}


def edens(f, x):
    o = f.exc_vxc(x)
    return o["zk"] * x["rho"].sum(1), o


def fd_check(f, x):
    e0, o = edens(f, x)
    vars_ = [("rho", "vrho")]
    if f.family in ("GGA", "MGGA"):
        vars_.append(("sigma", "vsigma"))
    if f.family == "MGGA":
        vars_.append(("tau", "vtau"))
        if f.flags & lxc.FLAG_NEEDS_LAPL:
            vars_.append(("lapl", "vlapl"))
    res = {}
    for xv, vv in vars_:
        for j in range(x[xv].shape[1]):
            if xv == "sigma" and f.nspin == 2 and j == 1:
                hs = np.sqrt(np.abs(x["sigma"][:, 0] * x["sigma"][:, 2])) + 1e-30
            elif xv == "lapl":
                hs = np.abs(x[xv][:, j]) + x["rho"].sum(1) ** (5 / 3)
            else:
                hs = np.abs(x[xv][:, j]) + 1e-30
            ana = o[vv][:, j]

            def dev(hrel):
                h = hrel * hs
                vals = []
                for k in (-2, -1, 1, 2):
                    xx = {kk: vv_.copy() for kk, vv_ in x.items()}
                    xx[xv][:, j] += k * h
                    vals.append(edens(f, xx)[0])
                fd = (vals[0] - 8 * vals[1] + 8 * vals[2] - vals[3]) / (12 * h)
                A, F = ana * hs, fd * hs
                scale = np.maximum.reduce([np.abs(A), np.abs(F), 1e-6 * np.abs(e0), np.full_like(A, 1e-14)])
                d = np.abs(A - F) / scale
                return np.where(np.isfinite(d), d, np.inf), fd

            d, fd = dev(1e-3)
            # step-size scan: a true inconsistency survives every h; an FD artifact does not
            dmin = d.copy()
            for hr in (1e-2, 1e-4, 1e-5, 1e-6):
                dmin = np.minimum(dmin, dev(hr)[0])
            persist = dmin > TOL
            bad = d > TOL
            if bad.any():
                res.setdefault("_persist", {})[f"{vv}[{j}]"] = {
                    "n": int(persist.sum()),
                    "max_dmin": float(dmin[persist].max()) if persist.any() else 0.0,
                }
                if persist.any():
                    d = np.where(persist, dmin, 0)
                    fd = dev(1e-3)[1]
                i = int(np.argmax(d))
                res[f"{vv}[{j}]"] = {
                    "frac_bad": float(bad.mean()),
                    "max_rel": float(d[i]) if np.isfinite(d[i]) else "inf",
                    "worst_point": {k: x[k][i].tolist() for k in KEYS[f.family]},
                    "ana": float(ana[i]), "fd": float(fd[i]),
                }
    return res


def edge_check(f, x):
    _, o = edens(f, x)
    bad = {}
    for k, v in o.items():
        m = ~np.isfinite(v)
        if m.ndim > 1:
            m = m.any(1)
        if m.any():
            i = int(np.argmax(m))
            bad[k] = {"n_nonfinite": int(m.sum()), "example": {kk: x[kk][i].tolist() for kk in x}}
    return bad


def main():
    fid, ns = int(sys.argv[1]), int(sys.argv[2])
    rec = {"id": fid, "name": lxc.name(fid), "nspin": ns}
    try:
        f = lxc.Func(fid, ns)
    except Exception as e:
        rec["skip"] = "init"; print(json.dumps(rec)); return
    rec["family"] = f.family
    if f.family == "OTHER" or not (f.flags & lxc.FLAG_VXC):
        rec["skip"] = "family/novxc"; print(json.dumps(rec)); return
    rng = np.random.default_rng(1234 + fid)
    rec["fd"] = fd_check(f, sample(rng, ns, NPT))
    rec["edge"] = edge_check(f, sample(rng, ns, 400, edge=True))
    f.close()
    print(json.dumps(rec))


if __name__ == "__main__":
    main()
