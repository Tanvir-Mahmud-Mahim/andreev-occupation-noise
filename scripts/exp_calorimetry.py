"""Fig. 4(a) data: linear-response energy resolution of the matched
5.3 x 0.5 um^2 design at 50 mK (the design of the photon simulation,
readout at 3 GHz) versus tau_A, from the coupled two-node
model (src/coupled.py): closed form (Andreev + phonon noise), with the
quantum-limited readout floor of an ideal fast readout, and with the
signal passing through a resonator of 1 MHz or 10 MHz linewidth.
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

OUT = os.path.join(os.path.dirname(__file__), "..", "data")
GHz = lambda J: J / H_PLANCK / 1e9
T0 = 0.05
g = lambda Tc: gap_bcs(T0, Tc, 1.7639 * KB * Tc) - 2.3994 * KB * T0
Tc_opt = brentq(g, T0 * 1.01, T0 * 6)
matched = Recipe("matched", "matched", Tc_opt, 1e-6, 0.5e-6, 5.3e-6,
                 30.0, 0.3, RECIPES[0].Ic20, RECIPES[0].Rn)
sb = SensorBudget(matched, nu_r=3e9)          # 3 GHz readout tone
s = sb.sj.andreev_sums(1e-4, T0, "L")
Ce, CA, G = sb.Ce(T0), s["C_A"], sb.Gep(T0)
kappa = 2 * np.pi * 1e6
nbar = 30.0
S_omega = (kappa / 4.0) ** 2 * 2.0 / (nbar * kappa)
conv = (sb.participation(T0) / 2) / sb.sj.dIdphi0(T0)
S_TA = S_omega / (2 * np.pi * sb.nu_r) ** 2 / (s["R_occ"] * conv) ** 2

tAs = np.logspace(-8.5, -4.5, 33)
res = dict(T0=T0, Tc_opt=Tc_opt, Ce_kB=Ce / KB, CA_kB=CA / KB,
           tauA_star_ns=CA / G * 1e9, tauAs=tAs.tolist())
res["closed"] = [GHz(cp.sigmaE_closed(T0, Ce, CA, G, t)) for t in tAs]
res["oneway"] = [GHz(cp.sigmaE_oneway_closed(T0, Ce, CA, G, t))
                 for t in tAs]
res["floor"] = [GHz(cp.sigmaE_numeric(T0, Ce, CA, G, t, S_TA_ro=S_TA))
                for t in tAs]
res["floor_k1"] = [GHz(cp.sigmaE_numeric(T0, Ce, CA, G, t, S_TA_ro=S_TA,
                                         kappa=2 * np.pi * 1e6))
                   for t in tAs]
res["floor_k10"] = [GHz(cp.sigmaE_numeric(T0, Ce, CA, G, t,
                                          S_TA_ro=S_TA,
                                          kappa=2 * np.pi * 1e7))
                    for t in tAs]
# check: numerical filter without readout floor equals the closed form
res["closed_check"] = max(abs(GHz(cp.sigmaE_numeric(T0, Ce, CA, G, t)) /
                              c - 1) for t, c in zip(tAs, res["closed"]))
i = int(np.argmin(np.abs(tAs - 1e-6)))
print(f"matched design at 50 mK: Tc*={Tc_opt:.3f} K, C_e={Ce/KB:.2f} kB,"
      f" C_A={CA/KB:.1f} kB, tauA*={CA/G*1e9:.0f} ns")
print(f"tauA=1us: closed {res['closed'][i]:.2f}, +floor "
      f"{res['floor'][i]:.2f}, +1 MHz {res['floor_k1'][i]:.2f}, "
      f"+10 MHz {res['floor_k10'][i]:.2f} GHz; closed-form check "
      f"{res['closed_check']:.1e}")
json.dump(res, open(os.path.join(OUT, "calorimetry.json"), "w"),
          indent=1)
