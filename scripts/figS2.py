"""Supplementary figure S2: self-consistent click Monte Carlo
(revision, 2026): d_KC versus exchange time for several designs, versus
resonator linewidth, and raw single-shot records."""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import matplotlib.pyplot as plt
from figstyle import C, DBL, panel_label

base = os.path.join(os.path.dirname(__file__), "..", "data")
R = json.load(open(os.path.join(base, "click_coupled.json")))
TR = np.load(os.path.join(base, "click_coupled_T_3e-08.npz"))

fig = plt.figure(figsize=(DBL, 4.3))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.9], hspace=0.50,
                      wspace=0.36, left=0.10, right=0.975, top=0.95,
                      bottom=0.09)

# (a) d_KC vs tauA ---------------------------------------------------------
ax = fig.add_subplot(gs[0, 0])
geo = [((0.05, 5.3, 0.5), C[0], "o", r"$5.3\times0.5$, 50 mK"),
       ((0.05, 2.0, 0.5), C[2], "s", r"$2\times0.5$, 50 mK"),
       ((0.05, 5.3, 1.5), C[3], "d", r"$5.3\times1.5$, 50 mK"),
       ((0.1, 5.3, 0.5), C[1], "v", r"$5.3\times0.5$, 100 mK")]
for (T0, W, L), c, mk, lab in geo:
    for scen, ls, mfc in (("T", "-", c), ("C", "--", "none")):
        pts = sorted([(v["tauA_ns"], v["dkc"]) for v in R["vs_tauA"]
                      if v["scen"] == scen and v["design"]["T0"] == T0
                      and v["design"]["W_um"] == W
                      and v["design"]["L_um"] == L])
        x, yv = zip(*pts)
        ax.semilogx(np.array(x) * 1e-9, yv, ls=ls, marker=mk, ms=3.5,
                    color=c, mfc=mfc, lw=1.0,
                    label=lab if scen == "T" else None)
ax.legend(loc="upper right", fontsize=6.3, handlelength=1.6,
          labelspacing=0.2, borderaxespad=0.2)
ax.text(0.03, 0.04, "solid: T,  dashed: C", transform=ax.transAxes,
        fontsize=7)
ax.set_xlabel(r"$\tau_{\rm A}$ (s)")
ax.set_ylabel(r"$d_{\rm KC}$ (26 GHz photon)")
ax.set_ylim(0, 6.2)
panel_label(ax, "(a)", dx=-0.25)

# (b) d_KC vs resonator linewidth -------------------------------------------
ax = fig.add_subplot(gs[0, 1])
for scen, c in (("T", C[0]), ("C", C[1])):
    for nadd, ls, mfc in ((0.0, "-", c), (14.0, "--", "none")):
        pts = sorted([(r["kappa_MHz"], r["dkc"]) for r in R["readout"]
                      if r["scen"] == scen and r["nadd"] == nadd])
        x, yv = zip(*pts)
        ax.semilogx(x, yv, ls=ls, marker="o", ms=3.5, color=c, mfc=mfc,
                    lw=1.0, label=rf"{scen}, $n_{{\rm add}}={nadd:g}$")
    ideal = [e["dkc"] for e in R["breakdown"]
             if e["step"] == "+ energy exchange (full)"
             and e["scen"] == scen and abs(e["tauA_ns"] - 30) < 1][0]
    ax.axhline(ideal, color=c, lw=0.7, ls=":")
ax.text(1.1, 4.25, "ideal readout (dotted)", fontsize=7, va="bottom")
ax.legend(loc="center right", fontsize=6.3, handlelength=1.6,
          labelspacing=0.2, borderaxespad=0.2)
ax.set_xlabel(r"resonator linewidth $\kappa/2\pi$ (MHz)")
ax.set_ylabel(r"$d_{\rm KC}$ ($5.3\times0.5\ \mu$m$^2$, 30 ns)")
ax.set_ylim(0, 5.0)
panel_label(ax, "(b)", dx=-0.25)

# (c) raw records ------------------------------------------------------------
ax = fig.add_subplot(gs[1, :])
t = TR["tgrid"] * 1e6
i0 = int(TR["i0"])
t = t - t[i0]
ker = np.ones(10) / 10.0
sm = lambda yy: np.convolve(yy, ker, mode="same")
sel = (t > -0.25) & (t < 1.6)
ax.plot(t[sel], sm(TR["y0"])[sel] / 1e6, color="0.6", lw=0.6,
        label="dark record")
ax.plot(t[sel], sm(TR["y1"])[sel] / 1e6, color=C[0], lw=0.6,
        label="photon record")
ax.plot(t[sel], TR["tpl"][sel] / 1e6, color=C[1], lw=1.4,
        label="noiseless click")
ax.legend(loc="lower right", ncol=3, fontsize=6.5, handlelength=1.4,
          borderaxespad=0.2, columnspacing=1.0)
ax.set_xlabel(r"time after absorption ($\mu$s)")
ax.set_ylabel(r"$\delta\nu_r$ (MHz), 5 ns average")
panel_label(ax, "(c)", dx=-0.09, dy=1.07)

fig.savefig(os.path.join(os.path.dirname(__file__), "..", "figures",
                         "figS2.pdf"))
print("figS2 done")
