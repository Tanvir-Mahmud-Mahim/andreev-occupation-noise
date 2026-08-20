"""Figure 2: predicted frequency-noise spectra, dT budgets, regime map."""
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

fig, axs = plt.subplots(1, 3, figsize=(DBL, 2.35))
plt.subplots_adjust(wspace=0.46, left=0.06, right=0.94, top=0.90,
                    bottom=0.19)

# (a) spectra -----------------------------------------------------------
ax = axs[0]
labels = ["Ta/Ti/Au", "Ti/Al/Au", "Ti/Al(thin)", "Ti/Al(thick)"]
for lab, c in zip(labels, C):
    sp = np.array(D["spectra"][lab]["Snu"])
    ax.loglog(freqs, np.sqrt(sp), color=c, lw=1.1)
ax.axhline(np.sqrt(D["Snu_floor"]), color="0.4", lw=0.8, ls=":")
ax.text(2.5, np.sqrt(D["Snu_floor"]) * 1.25, "quantum-limited readout",
        fontsize=6, color="0.35")
ax.annotate(r"knee $=1/2\pi\tau_{\rm A}$", xy=(1.6e5, 100),
            xytext=(1.1e3, 12), fontsize=6.5,
            arrowprops=dict(arrowstyle="->", lw=0.6))
ax.text(2.5, 230, "Ta/Ti/Au", color=C[0], fontsize=6.5)
ax.text(2.5, 350, "Ti/Al/Au", color=C[1], fontsize=6.5)
ax.text(2.5, 132, "Ti/Al(thin)", color=C[2], fontsize=6.5)
ax.text(2.5, 60, "Ti/Al(thick)", color=C[3], fontsize=6.5)
ax.set_xlabel(r"$f$ (Hz)")
ax.set_ylabel(r"$\sqrt{S_\nu}$ (Hz$/\sqrt{\rm Hz}$)")
ax.set_xlim(1, 1e7); ax.set_ylim(3, 700)
panel_label(ax, "(a)")

# (b) temperature-resolution budgets vs T (Ta/Ti/Au) --------------------
ax = axs[1]
sb = SensorBudget(RECIPES[0])
Ts = np.linspace(0.03, 0.42, 60)
for tA, c, lab in ((1e-8, C[0], "10 ns"), (1e-6, C[1], r"1 $\mu$s"),
                   (1e-4, C[2], r"100 $\mu$s")):
    dT = [sb.dT_andreev(T, tA, 1.0)[0] * 1e6 for T in Ts]
    ax.semilogy(Ts * 1e3, dT, color=c, lw=1.1)
    ax.text(390, dT[-1] * 1.15, lab, color=c, fontsize=6.2, ha="right")
ph = [sb.dT_phonon(T, 1.0) * 1e6 for T in Ts]
ax.semilogy(Ts * 1e3, ph, color="k", lw=1.2, ls="--")
ax.text(200, 1.1, "phonon TFN", fontsize=6.5, rotation=6)
ax.text(0.05, 0.93, "Ta/Ti/Au", transform=ax.transAxes, fontsize=7)
ax.set_xlabel(r"$T$ (mK)")
ax.set_ylabel(r"$\delta T_{\min}$ ($\mu$K, $t=1$ s)")
ax.set_ylim(0.03, 3e4)
panel_label(ax, "(b)")

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
ax.text(60, 3e-5, "Andreev\nlimited", fontsize=7, color="w",
        ha="center")
ax.text(340, 3e-9, "phonon\nlimited", fontsize=7, color="k",
        ha="center")
cb = fig.colorbar(pc, ax=ax, pad=0.02, aspect=28)
cb.set_label(r"$\log_{10}(\delta T_{\rm A}/\delta T_{\rm ph})$",
             fontsize=7)
cb.ax.tick_params(labelsize=6)
panel_label(ax, "(c)")

fig.savefig(os.path.join(os.path.dirname(__file__), "..", "figures",
                         "fig2.pdf"))
print("fig2 done")
