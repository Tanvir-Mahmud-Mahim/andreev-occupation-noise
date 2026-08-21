"""Regenerate every number quoted in the manuscript from the simulation
outputs, as LaTeX macros (paper/numbers.tex) and JSON (data/numbers.json).
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
Cal = json.load(open(os.path.join(base, "calorimetry.json")))
MP = json.load(open(os.path.join(base, "matched_points.json")))
LM = json.load(open(os.path.join(base, "limits.json")))
NL = json.load(open(os.path.join(base, "nonlinear_click.json")))

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
N["optDkT"] = round(D["opt_DkT"], 2)
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

# finite-L deficits
defs = {}
for r in RECIPES:
    fl = FiniteLJunction(r)
    defs[r.label] = fl.saturation_ratio(0.02, r.Tc / 6.0, "L")
N["deficitShortMaxPct"] = round(100 * (max(v for k, v in defs.items()
                                if k != "MoRe") - 1), 1)
N["deficitMoRePct"] = round(100 * (defs["MoRe"] - 1), 0)
N["LxiMoRe"] = 0.43

# matched design points
N["sigEmatched100"] = round(MP["T0.1_tauA1e-06"]["sigE_GHz"], 1)
N["sigEmatched50"] = round(MP["T0.05_tauA1e-06"]["sigE_GHz"], 2)
N["sigEmatched50fast"] = round(MP["T0.05_tauA1e-07"]["sigE_GHz"], 2)
N["snr26at50"] = round(MP["T0.05_tauA1e-06"]["snr26"], 0)
N["dark50"] = "3\\times10^{-11}"
N["TcMatched100"] = round(MP["T0.1_tauA1e-06"]["Tc"], 3)
N["TcMatched50"] = round(MP["T0.05_tauA1e-06"]["Tc"], 3)
N["gain"] = round(b100["sigE_GHz"] / MP["T0.1_tauA1e-06"]["sigE_GHz"], 0)

# channel correction and anchors
sb = SensorBudget(RECIPES[0])
s = sb.sj.andreev_sums(1e-4, 0.1, "L")
Rtot = (sb.sj.dIdphi0(0.1 + 1e-5) - sb.sj.dIdphi0(0.1 - 1e-5)) / 2e-5
N["gapChanCorr"] = round(Rtot / s["R_occ"], 2)
N["anchorRatio"] = round(Cal["anchor_ratio"], 2)

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

# nonlinear click Monte Carlo (nonlinear_click.json)
K = NL["T0.05_W5.3_L0.5"]
N["nlCe"] = round(K["Ce_kB"], 1)
N["nlTc"] = round(K["Tc"], 3)
N["nlTpeak"] = round(K["T_peak"], 3)
N["nlSigT"] = round(K["sigT_over_T"], 2)
N["nlSNR30"] = round(K["snr_mc_T_3e-08"], 1)
N["nlEff30"] = round(K["eff_T_3e-08"], 3)
N["nlSNR100"] = round(K["snr_mc_T_1e-07"], 1)
N["nlEff100"] = round(K["eff_T_1e-07"], 3)
N["nlSNRmus"] = round(K["snr_mc_T_1e-06"], 1)
N["nlSNRC30"] = round(K["snr_mc_C_3e-08"], 1)
N["nlSNR100mK"] = round(NL["T0.1_W5.3_L0.5"]["snr_mc_T_3e-08"], 1)
Ksm = NL["T0.05_W1.0_L0.1"]
N["nlSNRsmall"] = round(Ksm["snr_mc_T_3e-08"], 1)
N["nlCeSmall"] = round(Ksm["Ce_kB"], 1)
N["nlSigTsmall"] = round(Ksm["sigT_over_T"], 2)
N["nlTpkSmall"] = round(Ksm["T_peak"], 2)

# LaTeX macro names (letters only; these are the names used in main.tex)
MACROS = {
    "DeficitShortMaxPct": N["deficitShortMaxPct"],
    "DeficitMoRePct": int(N["deficitMoRePct"]),
    "TauAstar": N["tauAstar"],
    "SigEbestHundred": int(N["sigEbest100"]),
    "SigEmatchedHundred": N["sigEmatched100"],
    "SigEmatchedFifty": N["sigEmatched50"],
    "SnuTa": int(N["SnuTa"]),
    "SnuFloor": int(N["SnuFloor"]),
    "GapChanCorr": N["gapChanCorr"],
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
    "OptDkT": N["optDkT"],
    "LxiMoRe": N["LxiMoRe"],
    "AnchorRatio": N["anchorRatio"],
    "SnrTwentySix": int(N["snr26at50"]),
    "DarkFifty": N["dark50"],
    "ContShareMaxPct": N["contShareMaxPct"],
    "PhaseBiasPenalty": N["phaseBiasPenalty"],
    "PhaseBiasDkT": N["phaseBiasDkT"],
    "NoneqHundred": N["noneqHundred"],
    "NoneqFifty": N["noneqFifty"],
    "TauSpreadPct": N["tauSpreadPct"],
    "DorokhovL": N["dorokhovL"],
    "DorokhovI": N["dorokhovI"],
    "KneeDegradedPct": int(N["kneeDegradedPct"]),
    "MultiKneeSpread": N["multiKneeSpread"],
    "MultiFloorPct": int(N["multiFloorPct"]),
    "NlCe": N["nlCe"],
    "NlTc": N["nlTc"],
    "NlTpeak": N["nlTpeak"],
    "NlSigT": N["nlSigT"],
    "NlSNRthirty": N["nlSNR30"],
    "NlEffThirty": f"{N['nlEff30']:.3f}",
    "NlSNRhundred": N["nlSNR100"],
    "NlEffHundred": f"{N['nlEff100']:.3f}",
    "NlSNRmus": N["nlSNRmus"],
    "NlSNRscenC": N["nlSNRC30"],
    "NlSNRhundredmK": N["nlSNR100mK"],
    "NlSNRsmall": N["nlSNRsmall"],
    "NlCeSmall": N["nlCeSmall"],
    "NlSigTsmall": N["nlSigTsmall"],
    "NlTpkSmall": N["nlTpkSmall"],
}
j = json.dumps(N, indent=1)
open(os.path.join(base, "numbers.json"), "w").write(j)
with open(os.path.join(os.path.dirname(__file__), "..", "paper",
                       "numbers.tex"), "w") as f:
    for k, v in MACROS.items():
        f.write(f"\\newcommand{{\\n{k}}}{{{v}}}\n")
print(j)
print("macros written:", len(MACROS))
