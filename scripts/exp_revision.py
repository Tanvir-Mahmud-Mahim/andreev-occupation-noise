"""Additional analyses requested in review (revision, 2026).

(a) Electron-phonon exponent: crossover time and thermal time for
    delta = 3, 4, 5, 6, with Sigma_delta chosen so that the thermal
    conductance equals the measured delta = 3 value of Lee et al. at
    their 0.19 K operating point.
(b) Level-dependent (activated) exchange: floor for current readout at
    the phase of maximum current, and the gap optimum when the exchange
    time itself scales as exp(Delta*/kB T).
(c) Coupled two-node model: energy resolutions for the recipes and the
    matched designs (closed form and numerical filter with the readout
    floor and a finite resonator linewidth), the gain of the matched gap
    at fixed geometry, and the corner time of the noise spectrum. (The
    Fig. 4(a) curves are computed in exp_calorimetry.py.)
(d) Comparison of the occupation-noise plateau with the 1/f critical-
    current noise reported by Haque et al. (Sci. Rep. 11, 19900 (2021)):
    S_Ic/Ic^2 = a/f with a ~ 4e-6.
(e) Andreev transition frequencies relevant to readout back-action.
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from constants import KB, H_PLANCK
from materials import RECIPES, Recipe, heat_capacity, carrier_density
from sensor_limits import SensorBudget
from abs_model import gap_bcs
import coupled as cp

OUT = os.path.join(os.path.dirname(__file__), "..", "data")
GHz = lambda J: J / H_PLANCK / 1e9
res = {}
TA = RECIPES[0]


def matched(T0, W, L, y=2.3994):
    g = lambda Tc: gap_bcs(T0, Tc, 1.7639 * KB * Tc) - y * KB * T0
    Tc = brentq(g, T0 * 1.01, T0 * 40)
    return Recipe("m", "m", Tc, 1e-6, L * 1e-6, W * 1e-6, 30.0, 0.3,
                  TA.Ic20, TA.Rn)


def ro_floor_TA(sb, T, nu, kappa=2 * np.pi * 1e6, nbar=30.0, nadd=0.0):
    """White readout floor referred to the Andreev 'temperature' (K^2/Hz).
    """
    S_omega = (kappa / 4) ** 2 * 2 * (1 + 2 * nadd) / (nbar * kappa)
    S_ro_y = S_omega / (2 * np.pi * nu) ** 2
    s = sb.sj.andreev_sums(1e-4, T, "L")
    conv = (sb.participation(T) / 2) / sb.sj.dIdphi0(T)
    return S_ro_y / (s["R_occ"] * conv) ** 2


# ---- (a) electron-phonon exponent ---------------------------------------
T_anchor = 0.19
sb3 = SensorBudget(TA)
G_anchor = sb3.Gep(T_anchor)                 # delta = 3, Sigma = 2
ep = []
for T in (0.1, TA.Tc / 6.0, 0.05):
    s = sb3.sj.andreev_sums(1e-4, T, "L")
    for d in (3, 4, 5, 6):
        Sig = G_anchor / (d * TA.area * T_anchor ** (d - 1))
        sb = SensorBudget(TA, sigma=Sig, delta=d)
        ep.append(dict(T=T, delta=d, Sigma=Sig, G=sb.Gep(T),
                       tau_th_ns=sb.tau_th(T) * 1e9,
                       tauA_star_ns=s["C_A"] / sb.Gep(T) * 1e9,
                       dT_A_uK=sb.dT_andreev(T, 1e-6, 1.0)[0] * 1e6,
                       dT_ph_uK=sb.dT_phonon(T, 1.0) * 1e6))
res["eph_exponent"] = ep

# ---- (b) activated exchange ---------------------------------------------
T = 0.1
sb = SensorBudget(TA)
phm = sb.sj.phi_max(T)
E_m = sb.sj.E(phm, T)
D = sb.sj.Delta(T)
act = dict(T=T, phi_max=phm, E_over_Delta=E_m / D,
           depth_over_kT=(D - E_m) / (KB * T),
           tau_factor=float(np.exp((D - E_m) / (KB * T))),
           floor_factor=float(np.exp((D - E_m) / (2 * KB * T))))
# gap optimum: tauA = tau0 (fixed) versus tauA = tau0' exp(Delta(T)/kB T)
Tcs = np.linspace(0.101, 0.7, 600)
fixed, thermal = [], []
for Tc in Tcs:
    rr = Recipe("s", "s", Tc, 1e-6, TA.L, TA.W, TA.Vbg, TA.tau, TA.Ic20,
                TA.Rn)
    sbb = SensorBudget(rr)
    s = sbb.sj.andreev_sums(1e-4, T, "L")
    y = sbb.sj.Delta(T) / (KB * T)
    fixed.append(1.0 / s["C_A"])                 # floor^2 ~ tau/C_A
    thermal.append(np.exp(y) / s["C_A"])
fixed, thermal = np.array(fixed), np.array(thermal)
i_f, i_t = int(np.argmin(fixed)), int(np.argmin(thermal))
yy = lambda Tc: gap_bcs(T, Tc, 1.7639 * KB * Tc) / (KB * T)
act.update(opt_Tc_fixed=float(Tcs[i_f]), opt_y_fixed=float(yy(Tcs[i_f])),
           opt_Tc_thermal=float(Tcs[i_t]),
           opt_y_thermal=float(yy(Tcs[i_t])),
           # floor at the thermal-exchange optimum relative to the
           # floor at y = 2.40, both with exchange time ~ exp(y)
           thermal_gain_vs_matched=float(np.sqrt(thermal[i_f] /
                                                 thermal[i_t])))
res["activated"] = act

# ---- (c) coupled two-node model -----------------------------------------
rec = []
for r in RECIPES:
    for T in (round(r.Tc / 6.0, 4), 0.1):
        sb = SensorBudget(r)
        s = sb.sj.andreev_sums(1e-4, T, "L")
        Ce, CA, G = sb.Ce(T), s["C_A"], sb.Gep(T)
        rec.append(dict(label=r.label, T=T, Ce_kB=Ce / KB, CA_kB=CA / KB,
                        tauA_star_ns=CA / G * 1e9,
                        x_1us=cp.x_ratio(G, 1e-6, CA),
                        sigE_closed_GHz=GHz(cp.sigmaE_closed(T, Ce, CA, G,
                                                             1e-6)),
                        sigE_oneway_GHz=GHz(cp.sigmaE_oneway_closed(
                            T, Ce, CA, G, 1e-6)),
                        corner_ns=cp.corner_time(T, Ce, CA, G, 1e-6) *
                        1e9))
res["recipes_coupled"] = rec

# matched gap at fixed Ta/Ti/Au geometry, 100 mK
gain = {}
for tA in (1e-7, 1e-6, 1e-5):
    out = []
    for r_ in (TA, matched(0.1, TA.W * 1e6, TA.L * 1e6)):
        sb = SensorBudget(r_)
        s = sb.sj.andreev_sums(1e-4, 0.1, "L")
        out.append(cp.sigmaE_closed(0.1, sb.Ce(0.1), s["C_A"],
                                    sb.Gep(0.1), tA))
    gain[f"{tA:.0e}"] = float(out[0] / out[1])
res["matched_gain_fixed_geometry"] = gain

# scan sigma_E vs Tc* at 100 mK (Fig. 3a) with the closed two-node form
scans = {}
Tcs2 = np.linspace(0.105, 0.7, 120)
for tA in (1e-7, 1e-6, 1e-5):
    sig = []
    for Tc in Tcs2:
        rr = Recipe("s", "s", Tc, 1e-6, TA.L, TA.W, TA.Vbg, TA.tau,
                    TA.Ic20, TA.Rn)
        sbb = SensorBudget(rr)
        s = sbb.sj.andreev_sums(1e-4, 0.1, "L")
        sig.append(GHz(cp.sigmaE_closed(0.1, sbb.Ce(0.1), s["C_A"],
                                        sbb.Gep(0.1), tA)))
    scans[f"{tA:.0e}"] = sig
i = int(np.argmin(scans["1e-06"]))
res["scan"] = dict(Tcs=Tcs2.tolist(), scans=scans,
                   opt_Tc=float(Tcs2[i]),
                   opt_y=float(yy(Tcs2[i])),
                   opt_sigE_GHz=float(scans["1e-06"][i]))

# matched designs: closed, + readout floor, + 1 MHz linewidth. Both keep
# C_e near 10 kB (relative temperature fluctuations about one third); a
# 1 x 0.1 um^2 channel has C_e < 1 kB and is outside linear response.
md = []
for (T0, W, L, nu) in ((0.1, 5.3, 0.2, 6e9), (0.05, 5.3, 0.5, 3e9)):
    rr = matched(T0, W, L)
    sb = SensorBudget(rr, nu_r=nu)
    s = sb.sj.andreev_sums(1e-4, T0, "L")
    Ce, CA, G = sb.Ce(T0), s["C_A"], sb.Gep(T0)
    STA = ro_floor_TA(sb, T0, nu)
    for tA in (3e-8, 1e-7, 1e-6):
        md.append(dict(T0=T0, W_um=W, L_um=L, nu_GHz=nu / 1e9, Tc=rr.Tc,
                       Ce_kB=Ce / KB, CA_kB=CA / KB, tauA_ns=tA * 1e9,
                       closed=GHz(cp.sigmaE_closed(T0, Ce, CA, G, tA)),
                       oneway=GHz(cp.sigmaE_oneway_closed(T0, Ce, CA, G,
                                                          tA)),
                       floor=GHz(cp.sigmaE_numeric(T0, Ce, CA, G, tA,
                                                   S_TA_ro=STA)),
                       floor_k1=GHz(cp.sigmaE_numeric(
                           T0, Ce, CA, G, tA, S_TA_ro=STA,
                           kappa=2 * np.pi * 1e6)),
                       corner_ns=cp.corner_time(T0, Ce, CA, G, tA) * 1e9,
                       pair_GHz=2 * sb.sj.Delta(T0) / H_PLANCK / 1e9))
res["matched_designs"] = md

# ---- (d) 1/f comparison --------------------------------------------------
A_HAQUE = 4e-6
one_f = []
for r in RECIPES:
    sb = SensorBudget(r)
    for T in (r.Tc / 6.0, r.Tc / 5.0, r.Tc / 4.0, r.Tc / 3.0):
        s = sb.sj.andreev_sums(1e-4, T, "L")
        I1 = sb.sj.dIdphi0(T)
        # plateau of S_I'/I'^2 per unit tauA; ratio to a/f at the corner
        # f_c = 1/(2 pi tauA) is independent of tauA
        ratio = s["S0_over_tau"] / I1 ** 2 / (A_HAQUE * 2 * np.pi)
        one_f.append(dict(label=r.label, T=T, Tc_over_T=r.Tc / T,
                          ratio_at_corner=float(ratio)))
res["one_over_f"] = one_f
# temperature at which the occupation-noise plateau of Ta/Ti/Au equals
# the a/f level at the corner frequency
sbT = SensorBudget(TA)


def lratio(T):
    s = sbT.sj.andreev_sums(1e-4, T, "L")
    return np.log(s["S0_over_tau"] / sbT.sj.dIdphi0(T) ** 2 /
                  (A_HAQUE * 2 * np.pi))


Tx = brentq(lratio, TA.Tc / 6.0, TA.Tc / 2.5)
res["one_over_f_cross"] = dict(T=Tx, Tc_over_T=TA.Tc / Tx)
mres = []
for (T0, W, L) in ((0.1, 5.3, 0.2), (0.05, 5.3, 0.5)):
    rr = matched(T0, W, L)
    sb = SensorBudget(rr)
    s = sb.sj.andreev_sums(1e-4, T0, "L")
    I1 = sb.sj.dIdphi0(T0)
    mres.append(dict(T0=T0, W_um=W, L_um=L, ratio_at_corner=float(
        s["S0_over_tau"] / I1 ** 2 / (A_HAQUE * 2 * np.pi))))
res["one_over_f_matched"] = mres

# ---- (e) Andreev transition frequencies -----------------------------------
tr = []
for r in RECIPES:
    sb = SensorBudget(r)
    for T in (r.Tc / 6.0, 0.1):
        D = sb.sj.Delta(T)
        tr.append(dict(label=r.label, T=T,
                       pair_phi0_GHz=2 * D / H_PLANCK / 1e9,
                       pair_phimax_GHz=2 * sb.sj.E(sb.sj.phi_max(T), T) /
                       H_PLANCK / 1e9))
res["transitions"] = tr

# ---- (f) energy balance of the mean click response (coupled model) -----
from click_sim import ClickModel
r5 = matched(0.05, 5.3, 0.5)
sb5 = SensorBudget(r5, nu_r=3e9)
cm = ClickModel(0.05, sb5.sj.Delta(0.05), sb5.sj.Nch, sb5.Ce(0.05) / 0.05,
                sb5.sigma * r5.area, sb5.delta, 3e9, sb5.L_r,
                sb5.sj.dIdphi0(0.05))
Eg = H_PLANCK * 26e9
tg = np.concatenate([np.linspace(0, 50e-9, 5001)[:-1],
                     np.linspace(50e-9, 2e-5, 40000)])
eb = {}
for scen in ("C", "T"):
    U, n = cm.mean_trajectory(tg, 0, Eg, scen, 3e-8, 1.0)
    Te = cm.Te_of_U(U)
    P = cm.SA * (Te ** cm.d - cm.T0 ** cm.d)
    Eph = np.trapezoid(P, tg)
    UA = cm.Ns * cm.E * (n - cm.f0)
    left = (U[-1] - cm.U0) + UA[-1]
    eb[scen] = dict(residual=float(abs(Eph + left - Eg) / Eg),
                    max_share_andreev=float(UA.max() / Eg),
                    t_max_share_ns=float(tg[int(np.argmax(UA))] * 1e9))
res["energy_balance"] = eb

json.dump(res, open(os.path.join(OUT, "revision.json"), "w"), indent=1)
print("(f) energy balance:", eb)

# ---- summary ----------------------------------------------------------------
print("(a) e-ph exponent (Ta/Ti/Au):")
for e in ep:
    print(f"   T={e['T']:.4f} delta={e['delta']} tau_th={e['tau_th_ns']:.2f}"
          f" ns tauA*={e['tauA_star_ns']:.2f} ns dT_ph={e['dT_ph_uK']:.2f}"
          f" uK dT_A={e['dT_A_uK']:.1f} uK")
print("(b) activated:", {k: round(v, 4) for k, v in act.items()})
print("(c) recipes:")
for d in rec:
    print("   ", {k: (round(v, 4) if isinstance(v, float) else v)
                  for k, v in d.items()})
print("    matched gain at fixed geometry:", gain)
print("    scan optimum:", res["scan"]["opt_Tc"], res["scan"]["opt_y"],
      res["scan"]["opt_sigE_GHz"])
for d in md:
    print("   ", {k: (round(v, 3) if isinstance(v, float) else v)
                  for k, v in d.items()})
print("(d) 1/f ratio at corner:")
for d in one_f:
    print(f"   {d['label']:12s} Tc*/T={d['Tc_over_T']:.0f}: "
          f"{d['ratio_at_corner']:.3g}")
print("   matched:", mres)
print("   crossing:", res["one_over_f_cross"])
print("(e) transitions:")
for d in tr:
    print("   ", d)
