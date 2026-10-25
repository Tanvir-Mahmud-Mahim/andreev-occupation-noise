"""Matched-gap design points: at 100 mK the Ta/Ti/Au geometry
(5.3 x 0.2 um^2, C_e = 8.1 kB) and at 50 mK the 5.3 x 0.5 um^2 design of
the photon simulation (C_e = 10.1 kB). Both keep the relative
temperature fluctuations sqrt(kB/C_e) near one third; a 1 x 0.1 um^2
channel has C_e < 1 kB and is outside linear response.

Energy resolution of the coupled two-node model (src/coupled.py) with
Andreev and phonon noise and a quantum-limited readout floor (30
photons; the floor is that of a 1 MHz linewidth). Two readout cases:
an ideal readout that follows the signal, and a resonator of 1 MHz
linewidth through which the signal must pass. The readout tone is 3 GHz
at 50 mK and 6 GHz at 100 mK,
below the pair-transition frequency 2 Delta*/h (10 and 5 GHz).
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from constants import KB, H_PLANCK
from materials import Recipe, RECIPES
from sensor_limits import SensorBudget
from scipy.optimize import brentq
from abs_model import gap_bcs
import coupled as cp

out = {}
kappa = 2 * np.pi * 1e6
nbar = 30.0
S_omega = (kappa / 4.0) ** 2 * 2.0 / (nbar * kappa)
for T0, nu, W, L in ((0.05, 3e9, 5.3e-6, 0.5e-6),
                     (0.1, 6e9, 5.3e-6, 0.2e-6)):
    g = lambda Tc: gap_bcs(T0, Tc, 1.7639 * KB * Tc) - 2.3994 * KB * T0
    Tc = brentq(g, T0 * 1.01, T0 * 6)
    r = Recipe("m", "m", Tc, 1e-6, L, W, 30.0, 0.3,
               RECIPES[0].Ic20, RECIPES[0].Rn)
    sb = SensorBudget(r, nu_r=nu)
    s = sb.sj.andreev_sums(1e-4, T0, "L")
    conv = (sb.participation(T0) / 2) / sb.sj.dIdphi0(T0)
    S_TA = S_omega / (2 * np.pi * nu) ** 2 / (s["R_occ"] * conv) ** 2
    Ce, CA, G = sb.Ce(T0), s["C_A"], sb.Gep(T0)
    for tA in (1e-7, 1e-6):
        sig = cp.sigmaE_numeric(T0, Ce, CA, G, tA, S_TA_ro=S_TA)
        sigk = cp.sigmaE_numeric(T0, Ce, CA, G, tA, S_TA_ro=S_TA,
                                 kappa=kappa)
        out[f"T{T0}_tauA{tA:.0e}"] = dict(
            Tc=Tc, nu_GHz=nu / 1e9, W_um=W * 1e6, L_um=L * 1e6,
            Ce_kB=Ce / 1.380649e-23,
            sigT_over_T=float(np.sqrt(1.380649e-23 / Ce)), sigE_GHz=sig / H_PLANCK / 1e9,
            sigE_k1_GHz=sigk / H_PLANCK / 1e9,
            sigE_closed_GHz=cp.sigmaE_closed(T0, Ce, CA, G, tA) /
            H_PLANCK / 1e9,
            pair_GHz=2 * sb.sj.Delta(T0) / H_PLANCK / 1e9)
        print(f"T={T0} tauA={tA:.0e}: Tc*={Tc:.3f} readout {nu/1e9:.0f} "
              f"GHz: sigE = h x {sig/H_PLANCK/1e9:.2f} GHz (ideal "
              f"readout), {sigk/H_PLANCK/1e9:.2f} GHz (1 MHz linewidth)")
base = os.path.join(os.path.dirname(__file__), "..", "data")
json.dump(out, open(os.path.join(base, "matched_points.json"), "w"),
          indent=1)
