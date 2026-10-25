"""Figure 2: predicted frequency-noise spectra, temperature-resolution
budgets, and the regime map."""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import matplotlib.pyplot as plt
from figstyle import C, DBL, panel_label
from constants import KB
from materials import RECIPES
from sensor_limits import SensorBudget

D = json.load(open(os.path.join(os.path.dirname(__file__), "..",
                                "data", "design.json")))
freqs = np.array(D["freqs"])

fig, axs = plt.subplots(1, 3, figsize=(DBL, 2.25))
plt.subplots_adjust(wspace=0.55, left=0.085, right=0.91, top=0.89,
                    bottom=0.20)

# (a) spectra -----------------------------------------------------------
ax = axs[0]
labels = ["Ta/Ti/Au", "Ti/Al/Au", "Ti/Al(thin)", "Ti/Al(thick)"]
order = [1, 0, 2, 3]      # legend sorted by plateau height
for lab, c in zip(labels, C):
    sp = np.array(D["spectra"][lab]["Snu"])
    ax.loglog(freqs, np.sqrt(sp), color=c, lw=1.1, label=lab)
ax.axhline(np.sqrt(D["Snu_floor"]), color="0.4", lw=0.8, ls=":")
ax.text(2.0, np.sqrt(D["Snu_floor"]) * 1.12, "ideal readout",
        fontsize=7, color="0.35", va="bottom")
ax.annotate(r"$1/2\pi\tau_{\rm A}$", xy=(1.6e5, 150),
            xytext=(3e5, 420), fontsize=7, va="center", ha="center",
            arrowprops=dict(arrowstyle="->", lw=0.6, shrinkA=1,
                            shrinkB=1))
hnd, lbl = ax.get_legend_handles_labels()
ax.legend([hnd[i] for i in order], [lbl[i] for i in order],
          loc="lower left", fontsize=6.5, handlelength=1.2,
          borderaxespad=0.2, labelspacing=0.2)
ax.set_xlabel(r"$f$ (Hz)")
ax.set_ylabel(r"$\sqrt{S_\nu}$ (Hz$/\sqrt{\rm Hz}$)")
ax.set_xlim(1, 1e7); ax.set_ylim(3, 900)
ax.set_xticks([1, 1e2, 1e4, 1e6])
panel_label(ax, "(a)", dx=-0.30)

# (b) temperature-resolution budgets vs T (Ta/Ti/Au) --------------------
ax = axs[1]
sb = SensorBudget(RECIPES[0])
Ts = np.linspace(0.03, 0.42, 60)
for tA, c, lab in ((1e-8, C[0], "10 ns"), (1e-6, C[1], r"1 $\mu$s"),
                   (1e-4, C[2], r"100 $\mu$s")):
    dT = [sb.dT_andreev(T, tA, 1.0)[0] * 1e6 for T in Ts]
    ax.semilogy(Ts * 1e3, dT, color=c, lw=1.1)
    xl = 370 if tA == 1e-8 else 310
    off = 0.30 if tA == 1e-8 else 1.7
    va = "top" if tA == 1e-8 else "bottom"
    ax.text(xl, dT[np.argmin(np.abs(Ts * 1e3 - xl))] * off, lab,
            color=c, fontsize=7, ha="center", va=va)
ph = [sb.dT_phonon(T, 1.0) * 1e6 for T in Ts]
ax.semilogy(Ts * 1e3, ph, color="k", lw=1.2, ls="--")
ax.text(60, 0.95, "electron-phonon", fontsize=7, va="top")
ax.text(0.96, 0.95, "Ta/Ti/Au", transform=ax.transAxes, fontsize=7,
        ha="right", va="top")
ax.set_xlabel(r"$T$ (mK)")
ax.set_ylabel(r"$\delta T_{\min}$ ($\mu$K, $t=1$ s)")
ax.set_ylim(0.03, 3e4)
ax.set_xticks([0, 200, 400])
panel_label(ax, "(b)", dx=-0.30)

# (c) regime map --------------------------------------------------------
ax = axs[2]
M = np.array(D["map"]["log10ratio"])
Ts_m = np.array(D["map"]["T"]) * 1e3
tAs = np.array(D["map"]["tauA"])
pc = ax.pcolormesh(Ts_m, tAs, M, cmap="RdBu_r", vmin=-3, vmax=3,
                   shading="auto", rasterized=True)
cs = ax.contour(Ts_m, tAs, M, levels=[0], colors="k", linewidths=1.0)
ax.set_yscale("log")
ax.set_xlabel(r"$T$ (mK)")
ax.set_ylabel(r"$\tau_{\rm A}$ (s)")
ax.text(40, 1e-5, "Andreev\nlimited", fontsize=7.5, color="w",
        ha="left")
ax.text(330, 4e-9, "phonon\nlimited", fontsize=7.5, color="k",
        ha="center")
ax.set_xticks([0, 200, 400])
cb = fig.colorbar(pc, ax=ax, pad=0.03, aspect=25)
cb.set_label(r"$\log_{10}(\delta T_{\rm A}/\delta T_{\rm ph})$",
             fontsize=7.5)
cb.ax.tick_params(labelsize=7)
panel_label(ax, "(c)", dx=-0.42)

fig.savefig(os.path.join(os.path.dirname(__file__), "..", "figures",
                         "fig2.pdf"))
print("fig2 done")
