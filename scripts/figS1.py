"""Supplemental figure: Monte Carlo validation of the telegraph PSD."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import matplotlib.pyplot as plt
from figstyle import C, SGL
from montecarlo import telegraph_traces, psd_single_sided

rng = np.random.default_rng(7)
f0, tauA = 0.3, 1.0
dt, n, M = 0.05, 400_000, 24
tr = telegraph_traces(f0, tauA, dt, n, M, rng)
x = np.sum(1.0 - 2.0 * tr, axis=1)
fr, S = psd_single_sided(x, dt)
bins = np.logspace(np.log10(fr[1]), np.log10(fr.max()), 40)
idx = np.digitize(fr, bins)
fb, Sb = [], []
for i in range(1, len(bins)):
    m = idx == i
    if m.sum() > 2:
        fb.append(np.mean(fr[m])); Sb.append(np.mean(S[m]))
S_an = M * 16 * f0 * (1 - f0) * tauA / (1 + (2 * np.pi *
                                             np.array(fb) * tauA)**2)
fig, ax = plt.subplots(figsize=(SGL, 2.3))
plt.subplots_adjust(left=0.16, right=0.97, top=0.95, bottom=0.2)
ax.loglog(fb, Sb, 'o', ms=3, color=C[0], label='Monte Carlo')
ax.loglog(fb, S_an, '-', color=C[1], lw=1.2, label='analytic Lorentzian')
ax.set_xlabel('frequency (units of $1/\\tau_{\\rm A}$)')
ax.set_ylabel('$S_O$ (arb. units)')
ax.legend()
fig.savefig(os.path.join(os.path.dirname(__file__), "..", "figures",
                         "figS1.pdf"))
print("figS1 done")
