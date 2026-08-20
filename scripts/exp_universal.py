"""Fig. 1 data: universal responsivity of short-junction ABS ensembles.

Computes the dimensionless responsivity  -Tc dln I'(0)/dT  and
-Tc dln Ic/dT as functions of t = T/Tc for a family of transparencies,
the occupation/gap channel decomposition, and collapse metrics.
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from constants import KB, E_CHARGE, HBAR
from materials import Recipe, RECIPES
from short_junction import ShortJunction

OUT = os.path.join(os.path.dirname(__file__), "..", "data")

taus = [0.1, 0.3, 0.42, 0.53, 0.58, 0.78, 0.95, 1.0]
ts = np.linspace(0.06, 0.94, 89)


def make_sj(tau):
    r = Recipe("scan", "scan", 1.0, 1e-6, 0.2e-6, 2e-6, 30.0, tau,
               1e-6, 50.0)
    return ShortJunction(r)


res = {"t": ts.tolist(), "taus": taus, "hL": {}, "hI": {},
       "hL_occ": {}, "hL_gap": {}}
for tau in taus:
    sj = make_sj(tau)
    Tc = sj.Tc
    hL, hI, hLo, hLg = [], [], [], []
    for t in ts:
        T = t * Tc
        hL.append(-Tc * sj.dlnI1_dT(T, rel=1e-3))
        hI.append(-Tc * sj.dlnIc_dT(T, rel=1e-3))
        s = sj.andreev_sums(1e-4, T, "L")
        I1 = sj.dIdphi0(T)
        hLo.append(-Tc * s["R_occ"] / I1)
        hLg.append(hL[-1] - hLo[-1])
    res["hL"][str(tau)] = hL
    res["hI"][str(tau)] = hI
    res["hL_occ"][str(tau)] = hLo
    res["hL_gap"][str(tau)] = hLg

# transparency-splitting metrics (tau <= 0.78 band and tau -> 1 excess)
sub = [str(x) for x in taus if x <= 0.78]
H = np.array([res["hI"][k] for k in sub])
i35 = int(np.argmin(np.abs(ts - 0.35)))
i20 = int(np.argmin(np.abs(ts - 0.2)))
res["spread_t035"] = float((H[:, i35].max() - H[:, i35].min()) /
                           H[:, i35].mean())
res["tauSplit_t035"] = float(np.abs(np.array(res["hI"]["0.78"][i35]) /
                                    np.array(res["hI"]["0.3"][i35])))
H1 = np.array(res["hI"]["1.0"])
res["tau1_excess_at_0p2"] = float(H1[i20] / H[:, i20].mean())
res["matched_level_y"] = 2.3994

# activation energies E_A(tau) = Delta0 sqrt(1 - tau sin^2(phi_max/2))
act = {}
for tau in taus:
    sj = make_sj(tau)
    T = 0.1 * sj.Tc
    pm = sj.phi_max(T)
    act[str(tau)] = float(np.sqrt(1 - tau * np.sin(pm / 2)**2))
res["EA_over_Delta"] = act

with open(os.path.join(OUT, "universal.json"), "w") as f:
    json.dump(res, f)
print(f"tau split at t=0.35: {res['tauSplit_t035']:.2f}")
print(f"tau=1 excess factor at t=0.2: {res['tau1_excess_at_0p2']:.2f}")
print("E_A/Delta0:", {k: round(v, 3) for k, v in act.items()})
