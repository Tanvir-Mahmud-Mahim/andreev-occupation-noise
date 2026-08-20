"""Fig. 2 and Fig. 3 data: noise spectra, regime map, matched-level
design, recipe budgets, and the roadmap point.
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from constants import KB, H_PLANCK
from materials import RECIPES, Recipe
from sensor_limits import SensorBudget
from short_junction import ShortJunction

OUT = os.path.join(os.path.dirname(__file__), "..", "data")
res = {}

# ---- per-recipe budgets at T = Tc*/6 and at 0.1 K --------------------
tauA = 1e-6
budgets = []
for r in RECIPES:
    for T in sorted({round(r.Tc / 6.0, 4), 0.1}):
        sb = SensorBudget(r)
        s = sb.sj.andreev_sums(1e-4, T, "L")
        a, _ = sb.dT_andreev(T, tauA, 1.0)
        ph = sb.dT_phonon(T, 1.0)
        sig = sb.energy_resolution_analytic_A(T, tauA)
        Ic1 = sb.sj.Ic(T)
        fom = -sb.sj.dlnIc_dT(T)          # |dIc/dT| / Ic  (1/K)
        tstar = sb.tau_th(T) * s["C_A"] / sb.Ce(T)
        budgets.append(dict(
            label=r.label, T=T, DkT=r.Delta / (KB * T),
            CA_kB=s["C_A"] / KB, Ce_kB=sb.Ce(T) / KB,
            dT_A_uK=a * 1e6, dT_ph_uK=ph * 1e6,
            sigE_GHz=sig / H_PLANCK / 1e9, fom_perK=fom,
            tauA_star_ns=tstar * 1e9, LJ_nH=sb.LJ(T) * 1e9,
            tau_th_ns=sb.tau_th(T) * 1e9))
res["budgets"] = budgets

# ---- frequency-noise spectra (measurable prediction) ------------------
freqs = np.logspace(0, 7, 200)
spectra = {}
for r in RECIPES[:4]:
    sb = SensorBudget(r)
    T = r.Tc / 6.0
    Sy, Snu = sb.freq_noise_spectrum(T, tauA, freqs)
    spectra[r.label] = dict(T=T, Snu=Snu.tolist())
res["freqs"] = freqs.tolist()
res["spectra"] = spectra
# quantum-limited readout floor for reference (n=30, kappa/2pi=1 MHz)
kappa = 2 * np.pi * 1e6
nbar = 30.0
S_omega = (kappa / 4.0)**2 * 2.0 / (nbar * kappa)
res["Snu_floor"] = S_omega / (2 * np.pi)**2

# ---- regime map: log10(dT_A / dT_ph) over (T, tauA), Ta/Ti/Au ---------
r = RECIPES[0]
sb = SensorBudget(r)
Ts = np.linspace(0.02, 0.45, 60)
tAs = np.logspace(-9, -3, 61)
M = np.zeros((len(tAs), len(Ts)))
for j, T in enumerate(Ts):
    s = sb.sj.andreev_sums(1e-4, T, "L")
    ratio2 = (s and (2 * KB * T**2 / s["C_A"]))
    ph2 = 2 * KB * T**2 * sb.tau_th(T) / sb.Ce(T)
    for i, tA in enumerate(tAs):
        M[i, j] = 0.5 * np.log10(ratio2 * tA / ph2)
res["map"] = dict(T=Ts.tolist(), tauA=tAs.tolist(), log10ratio=M.tolist())

# ---- matched-level design scan: sigma_E vs Tc*/T at several tauA ------
T0 = 0.1
scans = {}
Tcs = np.linspace(0.105, 0.7, 120)
for tA in (1e-7, 1e-6, 1e-5):
    sig = []
    for Tc in Tcs:
        rr = Recipe("s", "s", Tc, 1e-6, RECIPES[0].L, RECIPES[0].W,
                    RECIPES[0].Vbg, RECIPES[0].tau, RECIPES[0].Ic20,
                    RECIPES[0].Rn)
        sbb = SensorBudget(rr)
        sig.append(sbb.energy_resolution_analytic_A(T0, tA) /
                   H_PLANCK / 1e9)
    scans[f"{tA:.0e}"] = sig
res["scan_Tcs"] = Tcs.tolist()
res["scan_T0"] = T0
res["scans"] = scans
best_i = int(np.argmin(scans["1e-06"]))
rr = Recipe("s", "s", Tcs[best_i], 1e-6, RECIPES[0].L, RECIPES[0].W,
            RECIPES[0].Vbg, RECIPES[0].tau, RECIPES[0].Ic20, RECIPES[0].Rn)
res["opt_Tc"] = Tcs[best_i]
res["opt_DkT"] = SensorBudget(rr).sj.Delta(T0) / (KB * T0)
res["opt_sigE_GHz"] = scans["1e-06"][best_i]

# ---- phase-bias route: sigma_E vs phi0 for high-tau junction ----------
phi_scan = {}
for tau in (0.78, 0.95, 0.99):
    rr = Recipe("p", "p", 0.57, 1e-6, RECIPES[0].L, RECIPES[0].W,
                RECIPES[0].Vbg, tau, RECIPES[0].Ic20, RECIPES[0].Rn)
    sbb = SensorBudget(rr)
    sj = sbb.sj
    phis = np.linspace(0.05, np.pi - 0.02, 80)
    sig = []
    for p0 in phis:
        s = sj.andreev_sums(p0, T0, "I")
        ST0 = s["S0_over_tau"] * tauA / s["R_occ"]**2
        sig.append(sbb.Ce(T0) * np.sqrt(ST0 / sbb.tau_th(T0)) /
                   H_PLANCK / 1e9)
    E0 = sj.E(phis, T0) / (KB * T0)
    phi_scan[str(tau)] = dict(phi=phis.tolist(), sigE=sig,
                              EkT=E0.tolist())
res["phi_scan"] = phi_scan

# ---- roadmap: matched-gap small-area design at 50 and 100 mK ----------
roadmap = []
for T in (0.05, 0.1):
    # choose Tc* to satisfy Delta(T) = 2.40 kB T (numerically)
    from scipy.optimize import brentq
    from abs_model import gap_bcs
    g = lambda Tc: gap_bcs(T, Tc, 1.7639 * KB * Tc) - 2.3994 * KB * T
    Tc = brentq(g, T * 1.01, T * 6)
    for (W, L) in ((1.0e-6, 0.1e-6), (5.3e-6, 0.2e-6)):
        rr = Recipe("r", "r", Tc, 1e-6, L, W, 30.0, 0.3,
                    RECIPES[0].Ic20, RECIPES[0].Rn)
        sbb = SensorBudget(rr)
        for tA in (1e-7, 1e-6):
            sig = sbb.energy_resolution_analytic_A(T, tA)
            roadmap.append(dict(T=T, Tc=Tc, W_um=W * 1e6, L_um=L * 1e6,
                                tauA_ns=tA * 1e9,
                                sigE_GHz=sig / H_PLANCK / 1e9))
res["roadmap"] = roadmap

with open(os.path.join(OUT, "design.json"), "w") as f:
    json.dump(res, f)

print("optimum Tc* =", res["opt_Tc"], " Delta(T)/kT =", res["opt_DkT"],
      " sigE =", res["opt_sigE_GHz"], "GHz")
for d in roadmap:
    print(d)
print("floor sqrt(Snu) =", np.sqrt(res["Snu_floor"]), "Hz/rtHz")
for lab, sp in spectra.items():
    print(lab, "sqrt(Snu(1Hz)) =", np.sqrt(sp["Snu"][0]), "Hz/rtHz at T=",
          sp["T"])
