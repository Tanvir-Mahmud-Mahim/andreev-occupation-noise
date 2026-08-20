"""Figure 4: speed-sensitivity trade-off and photon discrimination."""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import matplotlib.pyplot as plt
from figstyle import C, DBL, panel_label

base = os.path.join(os.path.dirname(__file__), "..", "data")
D = json.load(open(os.path.join(base, "calorimetry.json")))
mp = json.load(open(os.path.join(base, "matched_points.json")))

fig, axs = plt.subplots(1, 2, figsize=(DBL * 0.72, 2.35))
plt.subplots_adjust(wspace=0.34, left=0.085, right=0.985, top=0.90,
                    bottom=0.19)

# (a) sigma_E vs tauA ---------------------------------------------------
ax = axs[0]
tA = np.array(D["tauAs"])
ax.loglog(tA, D["sigE_ana"], color="0.5", lw=0.9, ls="--",
          label="Andreev only (analytic)")
ax.loglog(tA, D["sigE_num_nofloor"], color=C[0], lw=1.1,
          label="+ phonon TFN")
ax.loglog(tA, D["sigE_num"], color=C[1], lw=1.1,
          label="+ readout floor")
ax.axhline(26, color="k", lw=0.8, ls=":")
ax.text(1.1e-6, 29, "26 GHz photon", fontsize=6)
ax.text(2.5e-6, 3.0, r"$\propto\sqrt{\tau_{\rm A}}$", fontsize=7,
        color="0.4", rotation=22)
ax.legend(loc="upper left", handlelength=1.5, fontsize=5.6)
ax.set_xlabel(r"$\tau_{\rm A}$ (s)")
ax.set_ylabel(r"$\sigma_E/h$ (GHz)")
ax.set_ylim(0.8, 200)
ax.text(0.96, 0.05, "matched design,\n100 mK",
        transform=ax.transAxes, fontsize=6.2, ha="right")
panel_label(ax, "(a)")

# (b) discrimination histograms ----------------------------------------
ax = axs[1]
x = np.linspace(-14, 44, 800)
for key, c, lab in (("T0.1_tauA1e-06", C[1],
                     r"100 mK: $\sigma_E=h\times7.9$ GHz"),
                    ("T0.05_tauA1e-06", C[0],
                     r"50 mK: $\sigma_E=h\times1.5$ GHz")):
    s = mp[key]["sigE_GHz"]
    g0 = np.exp(-0.5 * (x / s)**2)
    g1 = np.exp(-0.5 * ((x - 26.0) / s)**2)
    ax.plot(x, g0, color=c, lw=1.1)
    ax.plot(x, g1, color=c, lw=1.1, ls="--")
    ax.fill_between(x, g0, color=c, alpha=0.12, lw=0)
    ax.fill_between(x, g1, color=c, alpha=0.12, lw=0)
ax.axvline(13, color="0.4", lw=0.7, ls=":")
ax.text(13.5, 1.02, "threshold", fontsize=6, color="0.35")
ax.text(0, 1.06, "no photon", fontsize=6.5, ha="center")
ax.text(26, 1.06, r"$h\times26$ GHz photon", fontsize=6.5, ha="center")
ax.text(0.02, 0.83, "50 mK", color=C[0], transform=ax.transAxes,
        fontsize=7)
ax.text(0.02, 0.60, "100 mK", color=C[1], transform=ax.transAxes,
        fontsize=7)
ax.text(0.98, 0.45,
        "50 mK: dark rate\n$3\\times10^{-11}$ s$^{-1}$\n"
        "$\\mathrm{SNR}=17$",
        transform=ax.transAxes, fontsize=6.2, ha="right")
ax.set_xlabel(r"matched-filter output $/h$ (GHz)")
ax.set_ylabel("probability density (norm.)")
ax.set_ylim(0, 1.18)
panel_label(ax, "(b)", dx=-0.11)

fig.savefig(os.path.join(os.path.dirname(__file__), "..", "figures",
                         "fig4.pdf"))
print("fig4 done")
