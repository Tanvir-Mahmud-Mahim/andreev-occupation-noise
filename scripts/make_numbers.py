"""Regenerate every number quoted in the manuscript from the simulation
outputs, as LaTeX macros (data/numbers.tex) and JSON (data/numbers.json).
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from constants import KB, H_PLANCK
from materials import RECIPES
from sensor_limits import SensorBudget
from finiteL import FiniteLJunction

base = os.path.join(os.path.dirname(__file__), "..", "data")
U = json.load(open(os.path.join(base, "universal.json")))
D = json.load(open(os.path.join(base, "design.json")))
MP = json.load(open(os.path.join(base, "matched_points.json")))
LM = json.load(open(os.path.join(base, "limits.json")))
CC = json.load(open(os.path.join(base, "click_coupled.json")))
RV = json.load(open(os.path.join(base, "revision.json")))

N = {}

# universality / responsivity
ts = np.array(U["t"])
i35 = int(np.argmin(np.abs(ts - 0.35)))
h03 = np.abs(U["hI"]["0.3"][i35])
h078 = np.abs(U["hI"]["0.78"][i35])
N["tauSplit"] = round(h078 / h03, 2)
N["tauOneExcess"] = round(U["tau1_excess_at_0p2"], 1)
N["EAlow"] = round(U["EA_over_Delta"]["0.3"], 2)
N["EAhigh"] = round(U["EA_over_Delta"]["0.78"], 2)
N["EAone"] = round(U["EA_over_Delta"]["1.0"], 2)
N["ystar"] = 2.40
N["optTc"] = round(D["opt_Tc"], 3)

# budgets at 100 mK (Ta/Ti/Au) and general
b100 = [b for b in D["budgets"] if b["label"] == "Ta/Ti/Au"
        and abs(b["T"] - 0.1) < 1e-6][0]
N["dTA100"] = round(b100["dT_A_uK"], 0)
N["dTph100"] = round(b100["dT_ph_uK"], 1)
N["ratio100"] = round(b100["dT_A_uK"] / b100["dT_ph_uK"], 0)
N["tauAstar"] = round(b100["tauA_star_ns"], 1)
N["sigEbest100"] = round(b100["sigE_GHz"], 0)
N["CAbest"] = round(b100["CA_kB"], 1)
N["Cebest"] = round(b100["Ce_kB"], 1)
N["tauTh100"] = round(b100["tau_th_ns"], 1)
# worst tauA* among the four low-gap recipes at 0.1 K
worst = min(b["tauA_star_ns"] for b in D["budgets"]
            if abs(b["T"] - 0.1) < 1e-6 and b["label"].startswith("Ti/Al"))
N["tauAstarWorstPs"] = round(worst * 1e3, 1)

# spectra
N["SnuTa"] = round(np.sqrt(D["spectra"]["Ta/Ti/Au"]["Snu"][0]), 0)
N["SnuFloor"] = round(np.sqrt(D["Snu_floor"]), 0)
N["SnuKneekHz"] = round(1.0 / (2 * np.pi * 1e-6) / 1e3, 0)
# exchange time at which the Lorentzian plateau (computed at tauA = 1 us,
# plateau ~ tauA) falls to the quantum-limited readout floor
N["tauVisNs"] = round(1e3 * D["Snu_floor"] /
                      D["spectra"]["Ta/Ti/Au"]["Snu"][0], 0)
# crossover exchange time of the second-lowest-gap recipe at 100 mK
N["tauAstarTiAlAuNs"] = round([b for b in D["budgets"]
                               if b["label"] == "Ti/Al/Au"
                               and abs(b["T"] - 0.1) < 1e-6][0]
                              ["tauA_star_ns"], 2)
N["occSharePct"] = None   # filled below from gapChanCorr

# finite-L deficits
defs = {}
for r in RECIPES:
    fl = FiniteLJunction(r)
    defs[r.label] = fl.saturation_ratio(0.02, r.Tc / 6.0, "L")
N["deficitShortMaxPct"] = round(100 * (max(v for k, v in defs.items()
                                if k != "MoRe") - 1), 1)
N["deficitMoRePct"] = round(100 * (defs["MoRe"] - 1), 0)
mr = [r for r in RECIPES if r.label == "MoRe"][0]
N["LxiMoRe"] = round(mr.L / mr.xi, 2)       # L and xi as tabulated by Jung et al.

# matched design points (coupled model, ideal readout / 1 MHz resonator)
N["sigEmatched100"] = round(MP["T0.1_tauA1e-06"]["sigE_GHz"], 1)
N["sigEmatched50"] = round(MP["T0.05_tauA1e-06"]["sigE_GHz"], 1)
N["sigEmatched100k"] = int(round(MP["T0.1_tauA1e-06"]["sigE_k1_GHz"]))
N["sigEmatched50k"] = int(round(MP["T0.05_tauA1e-06"]["sigE_k1_GHz"]))
N["ceMatched100"] = round(MP["T0.1_tauA1e-06"]["Ce_kB"], 1)
N["TcMatched100"] = round(MP["T0.1_tauA1e-06"]["Tc"], 3)
N["TcMatched50"] = round(MP["T0.05_tauA1e-06"]["Tc"], 3)
N["gainTen"] = round(RV["matched_gain_fixed_geometry"]["1e-06"], 1)

# channel correction and anchors
sb = SensorBudget(RECIPES[0])
s = sb.sj.andreev_sums(1e-4, 0.1, "L")
Rtot = (sb.sj.dIdphi0(0.1 + 1e-5) - sb.sj.dIdphi0(0.1 - 1e-5)) / 2e-5
N["gapChanCorr"] = round(Rtot / s["R_occ"], 2)
N["occSharePct"] = int(round(100.0 / (Rtot / s["R_occ"]), 0))

# phase-bias route: best sigma_E for tau=0.99 at 100 mK
pb = D["phi_scan"]["0.99"]
N["phaseBiasBest"] = round(min(pb["sigE"]), 0)

# resolved-limitation numbers (limits.json)
N["contShareMaxPct"] = round(100 * max(
    abs(v["share_dIdT"]) for k, v in LM["continuum_shares"].items()
    if k != "MoRe"), 1)
N["phaseBiasPenalty"] = int(round(LM["phase_bias_penalty_sigmaE"]))
N["phaseBiasDkT"] = round(LM["phase_bias_DkT"], 1)
import math
f100 = LM["noneq"]["T0.1"]["f"]
f50 = LM["noneq"]["T0.05"]["f"]
pen = lambda f, q: math.sqrt(q * (1 - q) / (f * (1 - f)))
N["noneqHundred"] = round(pen(f100, 1e-4), 1)
N["noneqFifty"] = int(round(pen(f50, 1e-4)))
N["tauSpreadPct"] = round(100 * (LM["tau_inhomogeneity"]["0.1"]
                                 ["deficit_L"] - 1), 1)
N["dorokhovL"] = round(LM["dorokhov"]["deficit_L"], 2)
N["dorokhovI"] = round(LM["dorokhov"]["deficit_I_phi2"], 2)
N["kneeRecs"] = int(min(int(k) for k, v in
                        LM["knee_measurement"].items() if v <= 0.10))
N["kneeDegradedPct"] = round(100 * LM["knee_degraded"]["ratio1.5_N1000"],
                             0)
N["multiKneeSpread"] = int(round(LM["knee_multilevel"]["spread"], -1))
N["multiFloorPct"] = round(100 * (np.sqrt(
    1 + LM["knee_multilevel"]["err_plateau_1000"]) - 1), 0)

# electron-phonon exponent (revision.json): tauA* at 100 mK
for d in (4, 5, 6):
    e = [x for x in RV["eph_exponent"] if x["delta"] == d and
         abs(x["T"] - 0.1) < 1e-9][0]
    N[f"tauAstarD{d}"] = round(e["tauA_star_ns"], 1) if d < 6 else \
        int(round(e["tauA_star_ns"]))
# noise-equivalent temperature of the occupation floor, Ta/Ti/Au, 100 mK
s100 = sb3 = SensorBudget(RECIPES[0]).sj.andreev_sums(1e-4, 0.1, "L")
for key, tA in (("netA", 1e-6), ("netAfast", 1e-8)):
    N[key] = int(round(np.sqrt(4 * KB * 0.01 * tA / s100["C_A"]) * 1e6))
# activated exchange and gap optimum for thermal-quasiparticle exchange
A = RV["activated"]
N["actDepth"] = round(A["depth_over_kT"], 2)
N["actTau"] = round(A["tau_factor"], 1)
N["actFloor"] = round(A["floor_factor"], 1)
N["optTcThermal"] = round(A["opt_Tc_thermal"], 3)
N["optYThermal"] = round(A["opt_y_thermal"], 1)
N["thermalGain"] = round(A["thermal_gain_vs_matched"], 1)
# 1/f comparison
oneF = [x for x in RV["one_over_f"] if x["label"] == "Ta/Ti/Au"]
r6 = [x for x in oneF if abs(x["Tc_over_T"] - 6) < 1e-6][0]
r3 = [x for x in oneF if abs(x["Tc_over_T"] - 3) < 1e-6][0]
N["oneFSix"] = int(round(1.0 / r6["ratio_at_corner"], -1))
N["oneFThird"] = round(r3["ratio_at_corner"], 1)
N["oneFCross"] = round(RV["one_over_f_cross"]["Tc_over_T"], 1)
N["oneFMatched"] = int(round(RV["one_over_f_matched"][0]
                             ["ratio_at_corner"]))

# self-consistent click Monte Carlo (click_coupled.json)
def bd(step, scen, tA):
    return [b for b in CC["breakdown"] if b["step"] == step and
            b["scen"] == scen and abs(b["tauA_ns"] - tA) < 1][0]


full = "+ energy exchange (full)"
N["kcFullT"] = round(bd(full, "T", 30)["dkc"], 1)
N["kcFullC"] = round(bd(full, "C", 30)["dkc"], 1)
N["kcFullTh"] = round(bd(full, "T", 100)["dkc"], 1)
N["kcFullCh"] = round(bd(full, "C", 100)["dkc"], 1)
N["effFullT"] = round(bd(full, "T", 30)["eff"], 2)
N["kcEarlierT"] = round(bd("earlier model, 6 GHz", "T", 30)["dkc"], 1)
N["kcEarlierC"] = round(bd("earlier model, 6 GHz", "C", 30)["dkc"], 1)
N["kcEarlierThreeT"] = round(bd("earlier model, 3 GHz", "T", 30)["dkc"], 1)
N["kcEarlierThreeC"] = round(bd("earlier model, 3 GHz", "C", 30)["dkc"], 1)
N["kcNonstT"] = round(bd("+ nonstationary noise", "T", 30)["dkc"], 1)
N["shiftEarlier"] = int(round(bd("earlier model, 3 GHz", "T", 30)
                              ["peak_dnu_MHz"]))
N["shiftFull"] = round(bd(full, "T", 30)["peak_dnu_MHz"], 1)
N["corrTime"] = int(round(bd(full, "T", 30)["corr_time_ns"]))
ro = lambda sc, k, na: [r for r in CC["readout"] if r["scen"] == sc and
                        abs(r["kappa_MHz"] - k) < 1e-6 and
                        r["nadd"] == na][0]
N["kcResFifty"] = round(ro("T", 50, 0.0)["dkc"], 1)
N["kcResHEMT"] = round(ro("T", 50, 14.0)["dkc"], 1)
N["kcBest"] = round(max([v["dkc"] for v in CC["vs_tauA"]] +
                        [v["dkc"] for v in CC["gap_scan"]]), 1)
N["kcHundredmK"] = float(np.ceil(10 * max(
    v["dkc"] for v in CC["vs_tauA"] if v["design"]["T0"] == 0.1)) / 10)
e3 = [e for e in CC["eph"] if e["design"]["delta"] == 3 and
      e["scen"] == "T"][0]["dkc"]
eb = max(e["dkc"] for e in CC["eph"] if e["design"]["delta"] in (4, 5)
         and e["scen"] == "T")
N["kcDgain"] = round(eb / e3, 1)
N["kcDfive"] = round(max(e["dkc"] for e in CC["eph"]
                         if e["design"]["delta"] in (4, 5)
                         and e["scen"] == "T"), 1)
EB = RV["energy_balance"]
N["ebResidual"] = float(f'{max(v["residual"] for v in EB.values()):.0e}')
N["ebShareLo"] = int(round(100 * EB["C"]["max_share_andreev"]))
N["ebShareHi"] = int(round(100 * EB["T"]["max_share_andreev"]))
N["ebTime"] = int(np.ceil(max(v["t_max_share_ns"] for v in EB.values())))
dm = CC["design_main"]
N["ceMain"] = round(dm["Ce_kB"], 1)
N["caMain"] = int(round(dm["CA_kB"]))
N["pairMain"] = round(dm["pair_GHz"], 1)
N["tauSlowMain"] = int(round(dm["tau_slow_ns"], -1))
# peak electron temperature of the 26 GHz deposit (energy balance)
from materials import heat_capacity
gam = dm["Ce_kB"] * KB / dm["T0"]
N["tpeakMain"] = round(float(np.sqrt(dm["T0"] ** 2 + 2 * H_PLANCK * 26e9
                                     / gam)), 3)

# LaTeX macro names (letters only; these are the names used in main.tex)
MACROS = {
    "DeficitShortMaxPct": N["deficitShortMaxPct"],
    "DeficitMoRePct": int(N["deficitMoRePct"]),
    "TauAstar": N["tauAstar"],
    "SigEmatchedHundred": N["sigEmatched100"],
    "SigEmatchedFifty": N["sigEmatched50"],
    "SigEmatchedHundredK": N["sigEmatched100k"],
    "SigEmatchedFiftyK": N["sigEmatched50k"],
    "GainTen": N["gainTen"],
    "SnuTa": int(N["SnuTa"]),
    "SnuFloor": int(N["SnuFloor"]),
    "EAlow": N["EAlow"],
    "EAhigh": N["EAhigh"],
    "EAone": N["EAone"],
    "TauSplit": N["tauSplit"],
    "TauOneExcess": int(round(N["tauOneExcess"])),
    "CAbest": N["CAbest"],
    "Cebest": N["Cebest"],
    "DTAHundred": int(N["dTA100"]),
    "DTphHundred": N["dTph100"],
    "RatioHundred": int(N["ratio100"]),
    "OptTc": N["optTc"],
    "LxiMoRe": N["LxiMoRe"],
    "ContShareMaxPct": N["contShareMaxPct"],
    "NoneqHundred": N["noneqHundred"],
    "NoneqFifty": N["noneqFifty"],
    "TauSpreadPct": N["tauSpreadPct"],
    "DorokhovL": N["dorokhovL"],
    "DorokhovI": N["dorokhovI"],
    "KneeDegradedPct": int(N["kneeDegradedPct"]),
    "MultiKneeSpread": N["multiKneeSpread"],
    "MultiFloorPct": int(N["multiFloorPct"]),
    "TauVisNs": int(N["tauVisNs"]),
    "TauAstarTiAlAu": N["tauAstarTiAlAuNs"],
    "OccSharePct": N["occSharePct"],
    "TauAstarDfour": N["tauAstarD4"],
    "TauAstarDfive": N["tauAstarD5"],
    "TauAstarDsix": N["tauAstarD6"],
    "NetA": N["netA"],
    "NetAfast": N["netAfast"],
    "ActDepth": N["actDepth"],
    "ActTau": N["actTau"],
    "ActFloor": N["actFloor"],
    "OptTcThermal": N["optTcThermal"],
    "OptYThermal": N["optYThermal"],
    "ThermalGain": N["thermalGain"],
    "OneFSix": N["oneFSix"],
    "OneFThird": N["oneFThird"],
    "OneFCross": N["oneFCross"],
    "OneFMatched": N["oneFMatched"],
    "KcFullT": N["kcFullT"],
    "KcFullC": N["kcFullC"],
    "KcFullTh": N["kcFullTh"],
    "KcFullCh": N["kcFullCh"],
    "EffFullT": f"{N['effFullT']:.2f}",
    "KcEarlierT": N["kcEarlierT"],
    "KcEarlierC": N["kcEarlierC"],
    "KcEarlierThreeT": N["kcEarlierThreeT"],
    "KcEarlierThreeC": N["kcEarlierThreeC"],
    "KcNonstT": N["kcNonstT"],
    "ShiftEarlier": N["shiftEarlier"],
    "ShiftFull": N["shiftFull"],
    "CorrTime": N["corrTime"],
    "KcResFifty": N["kcResFifty"],
    "KcResHEMT": N["kcResHEMT"],
    "KcBest": N["kcBest"],
    "KcHundredmK": N["kcHundredmK"],
    "KcDgain": N["kcDgain"],
    "KcDfive": N["kcDfive"],
    "EbShareLo": N["ebShareLo"],
    "EbShareHi": N["ebShareHi"],
    "EbTime": N["ebTime"],
    "CeMatchedHundred": N["ceMatched100"],
    "CeMain": N["ceMain"],
    "CAMain": N["caMain"],
    "PairMain": N["pairMain"],
    "TauSlowMain": N["tauSlowMain"],
    "TpeakMain": N["tpeakMain"],
}
j = json.dumps(N, indent=1)
open(os.path.join(base, "numbers.json"), "w").write(j)
with open(os.path.join(base, "numbers.tex"), "w") as f:
    for k, v in MACROS.items():
        f.write(f"\\newcommand{{\\n{k}}}{{{v}}}\n")
print(j)
print("macros written:", len(MACROS))
