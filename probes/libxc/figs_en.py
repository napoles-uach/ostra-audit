"""Figures for the libxc xc_integrate / lda_x_1d_soft finding (+ the other two classes)."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import k0, k1, modstruve
import lxc

# ---------- style (reference palette, light) ----------
SURF, INK, INK2, MUTED, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#c3c2b7"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
    "text.color": INK, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": "#ecebe8", "grid.linewidth": 0.8,
    "font.size": 10.5, "axes.titlesize": 11.5, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "lines.linewidth": 2, "legend.frameon": False,
})

# ---------- faithful Python port of libxc devel xc_integrate (tanh-sinh) ----------
EPS = np.finfo(float).eps
TMAX = np.arcsinh(-np.log(EPS) / np.pi)


def xc_integrate(f, a, b, fixed=False, log=None):
    tol = 2 * EPS; pi2 = np.pi / 2; mid = 0.5 * (a + b); rad = 0.5 * (b - a)
    stop = (lambda w, fv, T: w <= tol and w * fv <= tol * (1 + abs(T))) if fixed else \
           (lambda w, fv, T: w * fv <= tol * (1 + abs(T)))
    h = 1.0; T = pi2 * f(mid)
    if log is not None: log.append((0, mid))
    k = 1
    while True:
        t = k * h
        if t > TMAX: break
        w = pi2 * np.cosh(t) / np.cosh(pi2 * np.sinh(t)) ** 2; xt = rad * np.tanh(pi2 * np.sinh(t))
        fp, fm = f(mid + xt), f(mid - xt); T += w * (fp + fm)
        if log is not None: log += [(0, mid + xt), (0, mid - xt)]
        if stop(w, abs(fp) + abs(fm), T): break
        k += 1
    T *= h; Tp = T; dp = 0.0
    for level in range(1, 13):
        new = 0.0; h *= 0.5; k = 1
        while True:
            t = k * h
            if t > TMAX: break
            w = pi2 * np.cosh(t) / np.cosh(pi2 * np.sinh(t)) ** 2; xt = rad * np.tanh(pi2 * np.sinh(t))
            fp, fm = f(mid + xt), f(mid - xt); new += w * (fp + fm)
            if log is not None: log += [(level, mid + xt), (level, mid - xt)]
            if stop(w, abs(fp) + abs(fm), T): break
            k += 2
        T = 0.5 * T + h * new
        sc = 1 + abs(T) + abs(Tp); d = abs(T - Tp)
        if level >= 2 and d == 0.0: return rad * T
        if level >= 3 and dp > 0 and d * d <= tol * sc * dp: return rad * T
        dp = d; Tp = T
    return rad * T


f1 = lambda t: 2 * k0(t) if t > 0 else 0.0
f2 = lambda t: 2 * t * k0(t) if t > 0 else 0.0
I1 = np.vectorize(lambda R: np.pi * R * (k0(R) * modstruve(-1, R) + k1(R) * modstruve(0, R)) if R < 40 else np.pi)  # tail beyond 40 < 1e-17
I2 = lambda R: 2 * (1 - R * k1(R))

# sanity: port reproduces the C numbers
assert abs(xc_integrate(f1, 0, 77.5) - 0.006135923151881) < 1e-12
assert abs(xc_integrate(f1, 0, 77.5, fixed=True) - np.pi) < 1e-12

# ================= Fig 1: symptom =================
fig, ax = plt.subplots(figsize=(7.2, 4.2))
bkf = np.logspace(-1, 2.3, 600)
for beta, c in ((0.5, BLUE), (1.0, ORANGE), (2.0, AQUA)):
    rho = 2 * bkf / (np.pi * beta)
    f = lxc.Func(21, 1); lxc.L.xc_func_set_ext_params.argtypes = [lxc.ctypes.c_void_p, lxc.ctypes.POINTER(lxc.ctypes.c_double)]
    lxc.L.xc_func_set_ext_params(f.p, (lxc.ctypes.c_double * 1)(beta))
    zk = f.exc_vxc({"rho": rho[:, None]})["zk"]
    ax.plot(bkf, zk * beta, color=c, lw=2, label=f"libxc 7.1.2, β = {beta}")
R = 2 * bkf
exact = -2 * (I1(R) - I2(R) / R) / (4 * np.pi)          # beta * eps_x, universal in beta*kF
rd, zd = np.load("/mnt/user-data/uploads/libxc-devel/_claude_probe/soft1d_devel.npy")
ax.plot((np.pi*rd/2)[::12], zd[::12], "o", ms=4.5, mfc="none", mec=INK2, mew=1.1, label="libxc devel (caa72bd), β = 1")
ax.plot(bkf, exact, color=INK, lw=1.4, ls=(0, (4, 3)), label="exact reference (Bessel/Struve)")
ax.axvline(38.6, color=MUTED, lw=1, ls=":")
ax.text(40, -0.30, "β·k_F ≈ 38.6", color=INK2, fontsize=9.5)
ax.axhline(0, color=AXIS, lw=1)
ax.set_xscale("log"); ax.set_xlabel("β · k_F   (k_F = πρ/2)"); ax.set_ylabel("β · ε_x   (Ha·bohr)")
ax.set_title("1. lda_x_1d_soft: ε_x collapses and changes sign at high density\n    (all three β coincide: the failure depends only on β·k_F)")
ax.legend(loc="lower left", fontsize=9.5)
fig.tight_layout(); fig.savefig("fig1_sintoma_1d_soft_en.png", dpi=200); plt.close(fig)

# ================= Fig 2: mechanism =================
fig, axs = plt.subplots(2, 1, figsize=(7.2, 6.2), sharex=False)
for ax, Rv, tag in ((axs[0], 50.0, "R = 50: correct"), (axs[1], 100.0, "R = 100: fails")):
    log = []; val = xc_integrate(f1, 0, Rv, log=log)
    t = np.linspace(1e-4, Rv, 4000)
    ax.fill_between(t, 2 * k0(t), color=BLUE, alpha=0.18, lw=0)
    ax.plot(t, 2 * k0(t), color=BLUE, lw=2, label="integrand 2K₀(t)")
    xs = np.array([x for _, x in log]); lv = np.array([l for l, _ in log])
    ax.scatter(xs, np.full_like(xs, -0.35) - 0.12 * lv, s=14, color=ORANGE, zorder=3, label="evaluated abscissae (row = refinement level)")
    ax.set_xlim(-1, Rv + 1); ax.set_ylim(-0.35 - 0.12 * (lv.max() + 1), 4.2)
    ax.set_yticks([0, 1, 2, 3, 4])
    ax.set_title(f"{tag}   →   result {val:.6f}  (exact π ≈ 3.141593)")
    ax.set_ylabel("2K₀(t)")
    ax.text(Rv * 1.012, -0.35, "level 0", fontsize=8, color=INK2, va="center"); ax.text(Rv * 1.012, -0.35 - 0.12 * lv.max(), f"level {lv.max()}", fontsize=8, color=INK2, va="center")
axs[1].set_xlabel("t")
axs[1].annotate("at every refinement level (lower rows) the loop first\nevaluates the midpoint, where K₀ ≈ e^(−50) is 'negligible',\nand breaks: t ≈ 0 is never sampled again", xy=(51, -1.05), xytext=(56, 1.7), color=INK2, fontsize=9.5,
                arrowprops=dict(arrowstyle="->", color=MUTED, lw=1))
axs[0].legend(loc="upper right", fontsize=9)
fig.suptitle("2. Mechanism: the xc_integrate tail check assumes the integrand\n    is not negligible at the interval midpoint",
             x=0.02, ha="left", fontsize=11.5, fontweight="bold")
fig.tight_layout(); fig.savefig("fig2_mecanismo_integrador_en.png", dpi=200); plt.close(fig)

# ================= Fig 3: integral vs R, original vs fixed =================
Rs = np.logspace(0.5, 4, 220)
orig = np.array([xc_integrate(f1, 0, r) for r in Rs]); fix = np.array([xc_integrate(f1, 0, r, fixed=True) for r in Rs])
ex = np.array([I1(r) for r in Rs])
fig, axs = plt.subplots(1, 2, figsize=(10, 4))
ax = axs[0]
ax.plot(Rs, orig / np.pi, color=BLUE, label="xc_integrate, devel caa72bd")
ax.plot(Rs, fix / np.pi, color=ORANGE, ls=(0, (5, 2)), label="with the one-line fix")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(1e-6, 3); ax.axvline(77.2, color=MUTED, lw=1, ls=":")
ax.text(85, 2e-5, "R ≈ 77.2", color=INK2, fontsize=9)
ax.text(130, 2.5e-3, "plateaus π·2⁻⁹, π·2⁻¹⁰…:\neach level only halves T", color=INK2, fontsize=9)
ax.set_xlabel("R  (upper limit)"); ax.set_ylabel("(∫₀ᴿ 2K₀ dt) / π   (exact → 1)"); ax.set_title("3a. Integral value (normalized)")
ax.legend(fontsize=9, loc="lower left")
ax = axs[1]
eo = np.abs(orig - ex) / ex; ef = np.maximum(np.abs(fix - ex) / ex, 1e-17)
ax.plot(Rs, eo, color=BLUE, label="current"); ax.plot(Rs, ef, color=ORANGE, label="fixed")
ax.axhline(2 * EPS, color=MUTED, lw=1, ls=":"); ax.text(2000, 2 * EPS * 1.6, "2× machine ε", color=INK2, fontsize=9)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("R"); ax.set_ylabel("relative error")
ax.set_title("3b. Relative error vs. exact reference"); ax.legend(fontsize=9, loc="center left")
fig.tight_layout(); fig.savefig("fig3_correccion_en.png", dpi=200); plt.close(fig)

print("ok")
