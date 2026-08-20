"""Supplemental figure 2: nonlinear click Monte Carlo, SNR landscape
and raw single-shot records."""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import matplotlib.pyplot as plt
from figstyle import C, DBL, panel_label

base = os.path.join(os.path.dirname(__file__), "..", "data")
NL = json.load(open(os.path.join(base, "nonlinear_click.json")))
TR = np.load(os.path.join(base, "click_traces_3e-08.npz"))
tAs = [3e-8, 1e-7, 3e-7, 1e-6]

fig, axs = plt.subplots(1, 2, figsize=(DBL * 0.78, 2.5))
plt.subplots_adjust(wspace=0.30, left=0.075, right=0.985, top=0.87,
                    bottom=0.18)

# (a) SNR vs tauA ------------------------------------------------------
ax = axs[0]


def series(key, scen):
    return [NL[key][f"snr_mc_{scen}_{t:.0e}"] for t in tAs]


ax.loglog(tAs, series("T0.05_W5.3_L0.5", "T"), "o-", color=C[0],
          lw=1.1, ms=3, label=r"$5.3\times0.5\ \mu$m$^2$, T")
ax.loglog(tAs, series("T0.05_W5.3_L0.5", "C"), "o--", color=C[0],
          lw=1.0, ms=3, mfc="none", label=r"$5.3\times0.5\ \mu$m$^2$, C")
ax.loglog(tAs, series("T0.05_W1.0_L0.1", "T"), "s-", color=C[2],
          lw=1.0, ms=3, label=r"$1\times0.1\ \mu$m$^2$, T (invalid regime)")
ax.loglog(tAs, series("T0.05_W5.3_L1.5", "T"), "d-", color=C[3],
          lw=1.0, ms=3, label=r"$5.3\times1.5\ \mu$m$^2$, T")
ax.loglog(tAs, series("T0.1_W5.3_L0.5", "T"), "v-", color=C[1],
          lw=1.0, ms=3, label=r"$5.3\times0.5\ \mu$m$^2$, T, 100 mK")
ax.axhline(5.0, color="0.4", lw=0.7, ls=":")
ax.text(1.1e-6, 5.3, "SNR = 5", fontsize=6, color="0.35", ha="right")
ax.legend(loc="lower left", fontsize=5.6, handlelength=1.6,
          labelspacing=0.3)
ax.set_xlabel(r"$\tau_{\rm A}$ (s)")
ax.set_ylabel("Monte Carlo SNR (26 GHz photon)")
ax.set_ylim(0.1, 30)
panel_label(ax, "(a)", dx=-0.13)

# (b) raw records ------------------------------------------------------
ax = axs[1]
t_us = TR["tgrid"] * 1e6
n = len(t_us)
i0 = n // 4
tpl = np.zeros_like(t_us)
tpl[i0:] = TR["tpl"][:n - i0]
ker = np.ones(10) / 10.0
sm = lambda y: np.convolve(y, ker, mode="same")
ax.plot(t_us, sm(TR["y_no"]) / 1e6, color="0.65", lw=0.6,
        label="dark record")
ax.plot(t_us, sm(TR["y_photon"]) / 1e6, color=C[0], lw=0.6,
        label="photon record")
ax.plot(t_us, tpl / 1e6, color=C[1], lw=1.3, label="click template")
ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=3,
          fontsize=5.8, handlelength=1.4, frameon=False,
          borderaxespad=0.0, columnspacing=1.0)
ax.set_xlabel(r"$t$ ($\mu$s)")
ax.set_ylabel(r"$\delta\nu_r$ (MHz), 5 ns average")
panel_label(ax, "(b)", dx=-0.16)

fig.savefig(os.path.join(os.path.dirname(__file__), "..", "figures",
                         "figS2.pdf"))
print("figS2 done")
