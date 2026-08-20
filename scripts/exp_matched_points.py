"""Matched-gap design operating points at 50 and 100 mK (with readout floor)."""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from constants import KB, H_PLANCK
from materials import Recipe, RECIPES
from sensor_limits import SensorBudget
from scipy.optimize import brentq
from abs_model import gap_bcs

out = {}
kappa = 2 * np.pi * 1e6
nbar = 30.0
S_omega = (kappa / 4.0)**2 * 2.0 / (nbar * kappa)
for T0 in (0.05, 0.1):
    g = lambda Tc: gap_bcs(T0, Tc, 1.7639 * KB * Tc) - 2.3994 * KB * T0
    Tc = brentq(g, T0 * 1.01, T0 * 6)
    r = Recipe("m", "m", Tc, 1e-6, 0.1e-6, 1.0e-6, 30.0, 0.3,
               RECIPES[0].Ic20, RECIPES[0].Rn)
    sb = SensorBudget(r)
    S_ro_y = S_omega / (2 * np.pi * sb.nu_r)**2
    for tA in (1e-7, 1e-6):
        sig = sb.energy_resolution(T0, tA, S_ro_y=S_ro_y)
        snr = H_PLANCK * 26e9 / sig
        nu_eff = 1.0 / (2 * np.pi * max(tA, sb.tau_th(T0)))
        dark = nu_eff * np.exp(-0.5 * (0.5 * snr)**2)
        out[f"T{T0}_tauA{tA:.0e}"] = dict(
            Tc=Tc, sigE_GHz=sig / H_PLANCK / 1e9, snr26=snr, dark=dark)
        print(f"T={T0} tauA={tA:.0e}: Tc*={Tc:.3f} "
              f"sigE=h x {sig/H_PLANCK/1e9:.2f} GHz SNR26={snr:.1f} "
              f"dark={dark:.2e}/s")
base = os.path.join(os.path.dirname(__file__), "..", "data")
json.dump(out, open(os.path.join(base, "matched_points.json"), "w"))
