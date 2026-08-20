"""Figure 4: speed-sensitivity trade-off and nonlinear click Monte
Carlo (single-shot trace and matched-filter discrimination)."""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import matplotlib.pyplot as plt
from figstyle import C, DBL, panel_label

base = os.path.join(os.path.dirname(__file__), "..", "data")
D = json.load(open(os.path.join(base, "calorimetry.json")))
NL = json.load(open(os.path.join(base, "nonlinear_click.json")))
TR = np.load(os.path.join(base, "click_traces_3e-08.npz"))

fig, axs = plt.subplots(1, 3, figsize=(DBL, 2.35))
plt.subplots_adjust(wspace=0.42, left=0.065, right=0.985, top=0.90,
                    bottom=0.19)

# (a) sigma_E vs tauA (linear response) ---------------------------------
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
ax.text(0.96, 0.05, "matched design,\n100 mK, linear response",
        transform=ax.transAxes, fontsize=6.2, ha="right")
panel_label(ax, "(a)")

# (b) running matched-filter output (nonlinear MC, 50 mK) ---------------
ax = axs[1]
n = len(TR["tgrid"])
dt = float(TR["tgrid"][1] - TR["tgrid"][0])
lag = (np.arange(n) - n // 2) * dt * 1e6
rp = np.roll(TR["run_photon"], n // 2)
rn = np.roll(TR["run_no"], n // 2)
ax.plot(lag, rn, color="0.6", lw=0.8)
ax.plot(lag, rp, color=C[0], lw=0.9)
ax.axhline(0.5, color="k", lw=0.7, ls=":")
ax.text(-0.046, 0.53, "threshold", fontsize=6, color="0.25",
        va="bottom")
ax.text(0.012, 0.99, "photon", color=C[0], fontsize=6.5)
ax.text(0.085, -0.33, "no photon", color="0.45", fontsize=6.5)
ax.set_xlabel(r"arrival-time offset ($\mu$s)")
ax.set_ylabel("running filter output")
ax.set_xlim(-0.05, 0.15)      # window around the arrival time,
                              # clear of the circular-wrap alias
ax.set_ylim(-0.45, 1.28)
ax.text(0.97, 0.95, "50 mK, $C_e=10\\,k_{\\rm B}$, "
        r"$\tau_{\rm A}=30$ ns", transform=ax.transAxes,
        fontsize=6.2, ha="right", va="top")
panel_label(ax, "(b)", dx=-0.20)

# (c) matched-filter score histograms (nonlinear MC) --------------------
ax = axs[2]
key = "T0.05_W5.3_L0.5"
snr = NL[key]["snr_mc_T_3e-08"]
eff = NL[key]["eff_T_3e-08"]
bins = np.linspace(-0.55, 1.55, 57)
ax.hist(TR["out0"], bins=bins, color="0.6", alpha=0.75,
        label="dark (1000 trials)")
ax.hist(TR["out1"], bins=bins, color=C[0], alpha=0.75,
        label="photon (1000 trials)")
ax.axvline(0.5, color="k", lw=0.8, ls=":")
ax.text(0.475, 62, "threshold", fontsize=6, ha="right", va="center",
        color="0.25", rotation=90)
ax.text(0.03, 0.93, f"MC SNR = {snr:.1f}\n"
        f"efficiency = {eff:.3f}\n"
        "dark fraction = 0/1000",
        transform=ax.transAxes, fontsize=6.2, va="top")
ax.legend(loc="upper right", fontsize=5.6, handlelength=1.2,
          borderaxespad=0.3)
ax.set_xlabel("matched-filter score (photon units)")
ax.set_ylabel("trials per bin")
ax.set_ylim(0, 150)
panel_label(ax, "(c)", dx=-0.20)

fig.savefig(os.path.join(os.path.dirname(__file__), "..", "figures",
                         "fig4.pdf"))
print("fig4 done")
