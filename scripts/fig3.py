"""Figure 3: matched-level design law, phase-bias route, recipe chart."""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullFormatter
from figstyle import C, DBL, panel_label

base = os.path.join(os.path.dirname(__file__), "..", "data")
D = json.load(open(os.path.join(base, "design.json")))
mp = json.load(open(os.path.join(base, "matched_points.json")))

fig = plt.figure(figsize=(DBL, 4.2))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.82], hspace=0.62,
                      wspace=0.40, left=0.11, right=0.975, top=0.94,
                      bottom=0.10)

# (a) sigma_E vs Tc*/T scan --------------------------------------------
ax = fig.add_subplot(gs[0, 0])
Tcs = np.array(D["scan_Tcs"]) / D["scan_T0"]
for key, c, lab in (("1e-05", C[2], r"$\tau_{\rm A}=10\ \mu$s"),
                    ("1e-06", C[1], r"$\tau_{\rm A}=1\ \mu$s"),
                    ("1e-07", C[0], r"$\tau_{\rm A}=100$ ns")):
    ax.semilogy(Tcs, D["scans"][key], color=c, lw=1.1, label=lab)
iopt = int(np.argmin(D["scans"]["1e-06"]))
ax.plot(Tcs[iopt], D["scans"]["1e-06"][iopt], "o", color=C[1], ms=4)
ax.axvline(Tcs[iopt], color="0.6", lw=0.6, ls=":")
ax.annotate(r"$\Delta^*(T)=2.40\,k_{\rm B}T$",
            xy=(Tcs[iopt] + 0.05, D["scans"]["1e-06"][iopt] * 0.9),
            xytext=(2.3, 4.6), fontsize=7, ha="left", va="center",
            arrowprops=dict(arrowstyle="->", lw=0.6, shrinkB=2))
ax.axvline(5.7, color=C[3], lw=0.8, ls="--")
ax.text(5.55, 3500, "Ta/Ti/Au", color=C[3], fontsize=7,
        ha="right", va="center")
ax.legend(loc="upper left", bbox_to_anchor=(0.13, 1.0), fontsize=6.5,
          handlelength=1.2, labelspacing=0.2, borderaxespad=0.2)
ax.set_xlabel(r"$T_c^*/T$  ($T=100$ mK, Ta/Ti/Au geometry)")
ax.set_ylabel(r"$\sigma_E/h$ (GHz)")
ax.set_ylim(3, 2e4)
panel_label(ax, "(a)", dx=-0.27)

# (b) phase-bias route --------------------------------------------------
ax = fig.add_subplot(gs[0, 1])
for tau, c in (("0.78", C[0]), ("0.95", C[1]), ("0.99", C[2])):
    d = D["phi_scan"][tau]
    ax.semilogy(d["phi"], d["sigE"], color=c, lw=1.1,
                label=rf"$\tau={tau}$")
    E = np.array(d["EkT"])
    if E.min() < 2.4:
        j = int(np.argmin(np.abs(E - 2.4)))
        ax.plot(d["phi"][j], d["sigE"][j], "o", color=c, ms=4)
ax.yaxis.set_major_locator(FixedLocator([30, 100, 300]))
ax.yaxis.set_minor_formatter(NullFormatter())
ax.set_yticklabels(["30", "100", "300"])
ax.legend(loc="upper right", handlelength=1.2, fontsize=6.5,
          labelspacing=0.2, borderaxespad=0.2)
ax.text(0.04, 0.06, r"$\bullet$  $E(\varphi_0)=2.40\,k_{\rm B}T$",
        transform=ax.transAxes, fontsize=7)
ax.set_xlabel(r"phase bias $\varphi_0$ (rad)  ($T_c^*=0.57$ K)")
ax.set_ylabel(r"$\sigma_E/h$ (GHz)")
panel_label(ax, "(b)", dx=-0.27)

# (c) recipes and designs -----------------------------------------------
ax = fig.add_subplot(gs[1, :])
from materials import RECIPES
labels, vals = [], []
for r in RECIPES:                      # entry at T = Tc*/6 of each recipe
    b = [b for b in D["budgets"] if b["label"] == r.label and
         abs(b["T"] - round(r.Tc / 6.0, 4)) < 1e-9][0]
    labels.append(b["label"])
    vals.append(b["sigE_GHz"])
labels += ["matched\n100 mK", "matched\n50 mK"]
vals += [mp["T0.1_tauA1e-06"]["sigE_GHz"],
         mp["T0.05_tauA1e-06"]["sigE_GHz"]]
x = np.arange(len(labels))
cols = [C[0]] * (len(labels) - 2) + [C[1], C[1]]
ax.bar(x, vals, color=cols, width=0.62, log=True)
ax.set_xticks(x, labels, fontsize=7)
ax.axhline(26, color="k", lw=0.9, ls="--",
           label="energy of one 26 GHz photon")
ax.legend(loc="upper left", fontsize=7, handlelength=2.2,
          borderaxespad=0.4)


def fmt(v):
    if v < 1e3:
        return f"{v:.3g}"
    ex = int(np.floor(np.log10(v)))
    return rf"${v / 10**ex:.1f}\times10^{{{ex}}}$"


for xi, v in zip(x, vals):
    # label above each bar; the white box keeps it clear of the photon line
    # bars that end near the photon line get their label above the line
    ax.text(xi, max(v * 1.35, 26 * 1.5), fmt(v), fontsize=7,
            ha="center", va="bottom")
ax.set_ylabel(r"$\sigma_E/h$ (GHz)")
ax.set_ylim(0.5, 2e7)
ax.text(0.99, 0.95, r"$\tau_{\rm A}=1\ \mu$s;  recipes at $T=T_c^*/6$",
        transform=ax.transAxes, fontsize=7, va="top", ha="right")
panel_label(ax, "(c)", dx=-0.115, dy=1.08)

fig.savefig(os.path.join(os.path.dirname(__file__), "..", "figures",
                         "fig3.pdf"))
print("fig3 done")
