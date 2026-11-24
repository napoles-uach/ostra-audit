# Figure: error in atomization energies introduced by two XCFun defects (released - fixed), kcal/mol.
import json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
S = json.load(open("res/impact_summary.json"))
panels = [("OSTRA-014_r4scan", "r4SCAN exchange: damping parameters\n(XCFun 0.232/0.232 vs published 0.802/0.178)"),
          ("OSTRA-013_revtpssc", "revTPSS correlation: z² instead of z³")]
order = sorted(S["OSTRA-014_r4scan"]["AE"], key=lambda m: S["OSTRA-014_r4scan"]["AE"][m])
ink, mute, blue, surf = "#0b0b0b", "#52514e", "#2a78d6", "#fcfcfb"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": mute, "text.color": ink,
                     "axes.labelcolor": ink, "xtick.color": mute, "ytick.color": ink})
fig, axs = plt.subplots(1, 2, figsize=(7.2, 5.6), sharey=True, facecolor=surf)
for ax, (key, title) in zip(axs, panels):
    d = S[key]["AE"]; ys = range(len(order)); xs = [d[m] for m in order]
    ax.set_facecolor(surf)
    ax.axvspan(-1, 1, color="#e9e8e4", zorder=0, lw=0)
    ax.axvline(0, color=mute, lw=0.8, zorder=1)
    ax.hlines(list(ys), 0, xs, color=blue, lw=2, zorder=2)
    ax.plot(xs, list(ys), "o", ms=5, color=blue, mec=surf, mew=1.2, zorder=3)
    ax.set_title(title, fontsize=9, loc="left", color=ink)
    ax.set_xlim(-2.3, 1.2); ax.grid(axis="x", color="#e3e1dc", lw=0.6, zorder=0)
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.text(-0.97, len(order) - 0.4, "±1 kcal/mol", fontsize=8, color=mute, va="bottom")
axs[0].set_yticks(list(range(len(order)))); axs[0].set_yticklabels(order)
fig.supxlabel("Error in atomization energy (released XCFun − fixed), kcal/mol", fontsize=9, color=ink, y=0.045)
fig.text(0.01, 0.005, "PySCF 2.15.1, XCFun 2.1.1, def2-TZVP, grid level 5. r4SCAN: self-consistent; revTPSSC: first order on the revTPSS density.",
         fontsize=7, color=mute)
fig.tight_layout(rect=(0, 0.06, 1, 1))
fig.savefig("res/fig_xcfun_impact.png", dpi=200); fig.savefig("res/fig_xcfun_impact.svg")
