"""Figure 4: linear-response energy resolution and the self-consistent
single-photon click Monte Carlo (revision, 2026)."""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import matplotlib.pyplot as plt
from figstyle import C, DBL, panel_label

base = os.path.join(os.path.dirname(__file__), "..", "data")
D = json.load(open(os.path.join(base, "calorimetry.json")))
R = json.load(open(os.path.join(base, "click_coupled.json")))
TR = np.load(os.path.join(base, "click_coupled_T_3e-08.npz"))

fig = plt.figure(figsize=(DBL, 4.45))
gs = fig.add_gridspec(1, 2, wspace=0.38, left=0.105, right=0.975,
                      top=0.945, bottom=0.57)

# (a) sigma_E vs tauA (linear response, coupled model) -------------------
ax = fig.add_subplot(gs[0, 0])
tA = np.array(D["tauAs"])
ax.loglog(tA, D["closed"], color="k", lw=1.1,
          label="Andreev + phonon (closed form)")
ax.loglog(tA, D["floor"], color=C[0], lw=1.1, label="+ ideal readout")
ax.loglog(tA, D["floor_k10"], color=C[2], lw=1.1,
          label=r"+ resonator, $\kappa/2\pi=10$ MHz")
ax.loglog(tA, D["floor_k1"], color=C[1], lw=1.1,
          label=r"+ resonator, $\kappa/2\pi=1$ MHz")
ax.axhline(26, color="k", lw=0.8, ls=":")
ax.text(4e-9, 30, "26 GHz photon", fontsize=7, va="bottom")
ax.legend(loc="lower right", handlelength=1.3, fontsize=6.3,
          labelspacing=0.2, borderaxespad=0.2)
ax.set_xlabel(r"$\tau_{\rm A}$ (s)")
ax.set_ylabel(r"$\sigma_E/h$ (GHz)")
ax.set_ylim(0.3, 300)
ax.text(0.03, 0.96, r"$5.3\times0.5\ \mu$m$^2$, 50 mK",
        transform=ax.transAxes, fontsize=7, va="top")
panel_label(ax, "(a)", dx=-0.27)

# (b) score histograms, full model, scenario T, 30 ns ---------------------
ax = fig.add_subplot(gs[0, 1])
ex = R["example"]["T"]
bins = np.linspace(-1.0, 2.2, 49)
ax.hist(TR["out0"], bins=bins, color="0.6", alpha=0.75,
        label="dark")
ax.hist(TR["out1"], bins=bins, color=C[0], alpha=0.75,
        label="photon")
ax.axvline(0.5, color="k", lw=0.8, ls=":")
ax.text(0.03, 0.96, rf"$d_{{\rm KC}}={ex['dkc']:.1f}$" + "\n"
        + f"efficiency {ex['eff']:.2f}\n"
        + f"dark fraction {ex['dark']:.3f}",
        transform=ax.transAxes, fontsize=7, va="top")
ax.legend(loc="upper right", fontsize=6.5, handlelength=1.0,
          borderaxespad=0.2)
ax.set_xlabel("filter score (photon units)")
ax.set_ylabel(r"trials per bin (of $10^3$)")
ax.set_ylim(0, 1.25 * max(np.histogram(TR["out0"], bins)[0].max(),
                          np.histogram(TR["out1"], bins)[0].max()))
panel_label(ax, "(b)", dx=-0.27)

# (c) breakdown of d_KC from the earlier to the self-consistent model ----
gs2 = fig.add_gridspec(1, 1, left=0.44, right=0.975, top=0.40,
                       bottom=0.075)
ax = fig.add_subplot(gs2[0, 0])
rows = [("earlier model (6 GHz, ideal readout)", "earlier model, 6 GHz"),
        ("earlier model, 3 GHz readout", "earlier model, 3 GHz"),
        ("+ noise follows the hot electrons", "+ nonstationary noise"),
        ("+ energy exchange with Andreev levels",
         "+ energy exchange (full)")]
vals = {s: [] for s in ("T", "C")}
for _, key in rows:
    for s in ("T", "C"):
        v = [b["dkc"] for b in R["breakdown"] if b["step"] == key and
             b["scen"] == s and abs(b["tauA_ns"] - 30) < 1]
        vals[s].append(v[0])
labels = [r[0] for r in rows]
for nadd, lab in ((0.0, r"+ resonator, $\kappa/2\pi=50$ MHz, 30 photons"),
                  (14.0, r"+ amplifier with $n_{\rm add}=14$")):
    labels.append(lab)
    for s in ("T", "C"):
        v = [r["dkc"] for r in R["readout"] if r["scen"] == s and
             abs(r["kappa_MHz"] - 50) < 1 and r["nadd"] == nadd]
        vals[s].append(v[0])
y = np.arange(len(labels))[::-1]
h = 0.36
ax.barh(y + h / 2, vals["T"], height=h, color=C[0],
        label="scenario T (exchange speeds up)")
ax.barh(y - h / 2, vals["C"], height=h, color=C[1],
        label="scenario C (constant exchange)")
for yi, vt, vc in zip(y, vals["T"], vals["C"]):
    ax.text(vt + 0.08, yi + h / 2, f"{vt:.1f}", fontsize=7,
            va="center")
    ax.text(vc + 0.08, yi - h / 2, f"{vc:.1f}", fontsize=7,
            va="center")
ax.set_yticks(y, labels, fontsize=7)
ax.set_xlabel(r"$d_{\rm KC}$, 26 GHz photon ($5.3\times0.5\ \mu$m$^2$,"
              r" 50 mK, $\tau_{\rm A}=30$ ns)")
ax.set_xlim(0, 9.5)
ax.legend(loc="lower right", fontsize=6.5, handlelength=1.0,
          borderaxespad=0.2)
panel_label(ax, "(c)", dx=-0.80, dy=1.08)

fig.savefig(os.path.join(os.path.dirname(__file__), "..", "figures",
                         "fig4.pdf"))
print("fig4 done")
