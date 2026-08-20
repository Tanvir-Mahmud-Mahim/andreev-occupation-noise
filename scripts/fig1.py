"""Figure: universal responsivity family and channel split.
(The device schematic is its own figure, scripts/fig_device.py.)"""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import matplotlib.pyplot as plt
from figstyle import C, DBL, panel_label

D = json.load(open(os.path.join(os.path.dirname(__file__), "..",
                                "data", "universal.json")))
ts = np.array(D["t"])

fig, axs = plt.subplots(1, 2, figsize=(DBL * 0.72, 2.35))
plt.subplots_adjust(wspace=0.36, left=0.09, right=0.985, top=0.90,
                    bottom=0.19)

# (a) universal responsivity family -------------------------------------
ax = axs[0]
show = ["0.3", "0.58", "0.78", "0.95", "1.0"]
cols = [C[0], C[2], C[3], C[4], C[1]]
names = {"1.0": r"$\tau=1$", "0.95": r"$\tau=0.95$",
         "0.78": r"$\tau=0.78$", "0.58": r"$\tau=0.58$",
         "0.3": r"$\tau=0.3$"}
for k, c in zip(show, cols):
    h = np.abs(np.array(D["hI"][k]))
    ax.semilogy(ts, h, color=c, lw=1.1, label=names[k])
hnd, lbl = ax.get_legend_handles_labels()
ax.legend(hnd[::-1], lbl[::-1], loc="lower right", fontsize=5.8,
          handlelength=1.2, labelspacing=0.25, borderaxespad=0.3)
ax.set_xlabel(r"$T/T_c^*$")
ax.set_ylabel(r"$-T_c^*\,\mathrm{d}\ln I_c/\mathrm{d}T$")
ax.set_xlim(0.06, 0.94); ax.set_ylim(1e-6, 30)
panel_label(ax, "(a)")

# (b) channel decomposition for tau = 0.3 -------------------------------
ax = axs[1]
tot = np.abs(np.array(D["hL"]["0.3"]))
occ = np.abs(np.array(D["hL_occ"]["0.3"]))
gap = np.abs(np.array(D["hL_gap"]["0.3"]))
ax.semilogy(ts, tot, color="k", lw=1.2, label="total")
ax.semilogy(ts, occ, color=C[0], lw=1.1, label="occupation (noisy)")
ax.semilogy(ts, gap, color=C[1], lw=1.1, ls="--",
            label="gap (mean field)")
ax.legend(loc="lower right", handlelength=1.6)
ax.set_xlabel(r"$T/T_c^*$")
ax.set_ylabel(r"$-T_c^*\,\mathrm{d}\ln I'(0)/\mathrm{d}T$")
ax.set_xlim(0.06, 0.94); ax.set_ylim(1e-7, 30)
panel_label(ax, "(b)")

fig.savefig(os.path.join(os.path.dirname(__file__), "..", "figures",
                         "fig1.pdf"))
print("fig1 done")
