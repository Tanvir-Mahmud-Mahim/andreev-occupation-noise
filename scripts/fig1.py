"""Figure 1: concept, universal responsivity family, channel split."""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle, Circle
from figstyle import C, DBL, panel_label

D = json.load(open(os.path.join(os.path.dirname(__file__), "..",
                                "data", "universal.json")))
ts = np.array(D["t"])

fig, axs = plt.subplots(1, 3, figsize=(DBL, 2.35))
plt.subplots_adjust(wspace=0.42, left=0.055, right=0.985, top=0.90,
                    bottom=0.19)

# (a) schematic ---------------------------------------------------------
ax = axs[0]
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
# superconductors
ax.add_patch(Rectangle((0.2, 2.6), 2.2, 4.4, fc="#c8d6e5", ec="k", lw=0.6))
ax.add_patch(Rectangle((7.6, 2.6), 2.2, 4.4, fc="#c8d6e5", ec="k", lw=0.6))
ax.add_patch(Rectangle((2.4, 3.6), 5.2, 2.4, fc="#f3e2c7", ec="k", lw=0.6))
ax.text(1.3, 7.4, r"S ($\Delta^*$)", ha="center", fontsize=7)
ax.text(8.7, 7.4, r"S ($\Delta^*$)", ha="center", fontsize=7)
ax.text(5.0, 6.35, "graphene", ha="center", fontsize=6.5)
# gap edges and ABS ladder
ax.plot([2.6, 7.4], [5.55, 5.55], color="0.45", lw=0.7, ls=":")
ax.text(7.3, 5.72, r"$\Delta^*$", fontsize=6.5, color="0.3",
        ha="right")
for y, c in zip((5.15, 4.95, 4.72), (C[0], C[0], C[0])):
    ax.plot([3.0, 7.0], [y, y], color=c, lw=1.0)
ax.text(5.0, 4.15, r"Andreev levels $E_m(\varphi)$", ha="center",
        fontsize=6.5, color=C[0])
# occupation exchange arrows
ar = FancyArrowPatch((3.35, 5.20), (3.35, 5.92), arrowstyle="<->",
                     mutation_scale=6, color=C[1], lw=1.0)
ax.add_patch(ar)
ax.text(3.60, 5.35, r"$\tau_{\rm A}$", color=C[1], fontsize=7.5,
        ha="left", va="center")
# readout resonator
ax.plot([5.0, 5.0], [3.6, 2.2], color="k", lw=0.7)
th = np.linspace(0, 4 * np.pi, 200)
ax.plot(5.0 + 0.55 * np.sin(th), 1.75 - 0.14 * th / np.pi, color="k",
        lw=0.7)
ax.text(6.0, 1.3, r"$\nu_r \propto 1/\sqrt{L_J}$", fontsize=6.5)
ax.text(5.0, 9.3, r"$\sigma_m = 1-n_\uparrow-n_\downarrow$:"
        "  signal and noise", ha="center", fontsize=7)
panel_label(ax, "(a)", dx=0.0)

# (b) universal responsivity family -------------------------------------
ax = axs[1]
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
panel_label(ax, "(b)")

# (c) channel decomposition for tau = 0.3 -------------------------------
ax = axs[2]
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
panel_label(ax, "(c)")

fig.savefig(os.path.join(os.path.dirname(__file__), "..", "figures",
                         "fig1.pdf"))
print("fig1 done")
