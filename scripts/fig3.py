"""Figure 3: matched-level design law, phase-bias route, recipe chart."""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import matplotlib.pyplot as plt
from figstyle import C, DBL, panel_label

base = os.path.join(os.path.dirname(__file__), "..", "data")
D = json.load(open(os.path.join(base, "design.json")))

fig, axs = plt.subplots(1, 3, figsize=(DBL, 2.35))
plt.subplots_adjust(wspace=0.44, left=0.077, right=0.985, top=0.90,
                    bottom=0.19)

# (a) sigma_E vs Tc*/T scan --------------------------------------------
ax = axs[0]
Tcs = np.array(D["scan_Tcs"]) / D["scan_T0"]
for key, c, lab in (("1e-07", C[0], "100 ns"), ("1e-06", C[1],
                    r"1 $\mu$s"), ("1e-05", C[2], r"10 $\mu$s")):
    ax.semilogy(Tcs, D["scans"][key], color=c, lw=1.1)
    ax.text(Tcs[-1] - 0.15, D["scans"][key][-1] * 1.2, lab, color=c,
            fontsize=6.3, ha="right")
iopt = int(np.argmin(D["scans"]["1e-06"]))
ax.plot(Tcs[iopt], D["scans"]["1e-06"][iopt], "o", color=C[1], ms=4)
ax.axvline(Tcs[iopt], color="0.6", lw=0.6, ls=":")
ax.annotate(r"$\Delta^*(T)=2.40\,k_{\rm B}T$", xy=(Tcs[iopt] + 0.05,
            42), xytext=(2.75, 700), fontsize=7,
            arrowprops=dict(arrowstyle="->", lw=0.6, shrinkB=2))
ax.axvline(5.7, color=C[3], lw=0.8, ls="--")
ax.text(5.55, 12, "Ta/Ti/Au\n(measured)", color=C[3], fontsize=6,
        ha="right", va="bottom")
ax.set_xlabel(r"$T_c^*/T$   ($T=100$ mK)")
ax.set_ylabel(r"$\sigma_E/h$ (GHz)")
ax.set_ylim(8, 1e4)
panel_label(ax, "(a)", dx=-0.23)

# (b) phase-bias route --------------------------------------------------
ax = axs[1]
for tau, c in (("0.78", C[0]), ("0.95", C[1]), ("0.99", C[2])):
    d = D["phi_scan"][tau]
    ax.semilogy(d["phi"], d["sigE"], color=c, lw=1.1,
                label=rf"$\tau={tau}$")
    E = np.array(d["EkT"])
    if E.min() < 2.4:
        j = int(np.argmin(np.abs(E - 2.4)))
        ax.plot(d["phi"][j], d["sigE"][j], "o", color=c, ms=4)
ax.legend(loc="upper right", handlelength=1.4)
ax.text(0.03, 0.05, r"$\circ$: $E(\varphi_0)=2.40\,k_{\rm B}T$",
        transform=ax.transAxes, fontsize=6.5)
ax.set_xlabel(r"phase bias $\varphi_0$ (rad)")
ax.set_ylabel(r"$\sigma_E/h$ (GHz)")
panel_label(ax, "(b)")

# (c) recipes and designs -----------------------------------------------
ax = axs[2]
# take the T = Tc*/6 entry of each recipe (DkT = 1.764*6 there)
seen = {}
labels, vals = [], []
for b in D["budgets"]:
    if b["label"] not in seen and abs(b["DkT"] - 10.5) < 3.0:
        seen[b["label"]] = True
        labels.append(b["label"])
        vals.append(b["sigE_GHz"])
mp = json.load(open(os.path.join(base, "matched_points.json")))
labels += ["matched\n100 mK", "matched\n50 mK"]
vals += [mp["T0.1_tauA1e-06"]["sigE_GHz"], mp["T0.05_tauA1e-06"]["sigE_GHz"]]
y = np.arange(len(labels))[::-1]
cols = [C[0]] * (len(labels) - 2) + [C[1], C[1]]
ax.barh(y, vals, color=cols, height=0.62, log=True)
ax.set_yticks(y, labels, fontsize=6)
ax.axvline(26, color="k", lw=0.9, ls="--")
ax.annotate("26 GHz\nphoton", xy=(26, 0.42), xytext=(300, 0.38),
            fontsize=6, va="center",
            arrowprops=dict(arrowstyle="->", lw=0.6, shrinkB=1))


def fmt(v):
    if v < 1e3:
        return f"{v:.3g}"
    ex = int(np.floor(np.log10(v)))
    return rf"${v / 10**ex:.1f}\times10^{{{ex}}}$"


for yi, v in zip(y, vals):
    ax.text(v * 1.25, yi, fmt(v), fontsize=5.8, va="center")
ax.set_xlabel(r"$\sigma_E/h$ (GHz) at $T=T_c^*/6$, $\tau_{\rm A}=1\,\mu$s")
ax.set_xlim(0.5, 3e6)
panel_label(ax, "(c)", dx=-0.32)

fig.savefig(os.path.join(os.path.dirname(__file__), "..", "figures",
                         "fig3.pdf"))
print("fig3 done")
