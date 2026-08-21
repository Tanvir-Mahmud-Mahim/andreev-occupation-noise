"""Quantitative resolution of the remaining model limitations.

Produces data/limits.json with:
 1. pair-process worst case: tau_eff/tau_s versus pair rate (exact
    4-state master equation) at the operating occupation of the
    matched design; the single-quasiparticle value is the maximum.
 2. continuum share of I'(0) and of dI'(0)/dT for all six recipes at
    T = Tc*/6 (exact finite-length solver, bound/continuum split of
    the free energy).
 3. scenario-B (activated tau(E)) penalty of the phase-bias route at
    the matched-level point, and the statement that the recipe route
    is unaffected (levels at the gap edge).
 4. nonequilibrium occupation floor: resolution penalty factor versus
    steady-state trapped probability q for the Ta/Ti/Au recipe at 100
    and 50 mK.
 5. transparency-inhomogeneity saturation deficit: channels with
    Beta-distributed tau (mean 0.3), analytic short-junction sums,
    versus the closed form sqrt(1 + CV^2) at phi -> 0.
 6. calibration-scale independence: the bound-saturation ratio,
    matched-level optimum, and regime-map contour recomputed with
    scale = 1 (uncalibrated) versus calibrated.
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from constants import KB, H_PLANCK
from materials import RECIPES, Recipe
from short_junction import ShortJunction
from sensor_limits import SensorBudget
from noise_general import sigma_spectrum, tau_activated, noneq_penalty
from abs_model import continuum_share, gap_bcs
from scipy.optimize import brentq

OUT = os.path.join(os.path.dirname(__file__), "..", "data")
res = {}

# ---- 1. pair processes ------------------------------------------------
f_matched = 1.0 / (np.exp(2.3994) + 1.0)
ratios = {}
for gp_over_gs in (0.0, 0.3, 1.0, 3.0, 10.0):
    _, _, _, te = sigma_spectrum(f_matched, 1.0, gp_over_gs,
                                 np.array([0.0]))
    ratios[str(gp_over_gs)] = te
res["pair_tau_ratio"] = ratios
res["pair_worst_is_singles"] = bool(
    all(v <= ratios["0.0"] + 1e-12 for v in ratios.values()))

# ---- 2. continuum shares ---------------------------------------------
shares = {}
for r in RECIPES:
    T = r.Tc / 6.0
    s = continuum_share(r, T)
    shares[r.label] = dict(share_I=s["share_I"],
                           share_dIdT=s["share_dIdT"])
    print(f"{r.label:12s} continuum share at pi/2: I "
          f"{s['share_I']*100:6.2f}%  dI/dT {s['share_dIdT']*100:6.2f}%")
res["continuum_shares"] = shares

# ---- 3. scenario-B phase-bias penalty --------------------------------
T0 = 0.1
rr = Recipe("p", "p", 0.57, 1e-6, RECIPES[0].L, RECIPES[0].W,
            RECIPES[0].Vbg, 0.99, RECIPES[0].Ic20, RECIPES[0].Rn)
sj = ShortJunction(rr); sj.calibrate()
D = sj.Delta(T0)
E_star = 2.3994 * KB * T0
penalty = float(np.sqrt(tau_activated(E_star, D, T0, 1.0)))
res["phase_bias_penalty_sigmaE"] = penalty
res["phase_bias_DkT"] = float(D / (KB * T0))
print(f"scenario-B phase-bias penalty on sigma_E at matched level: "
      f"x{penalty:.0f}  (Delta/kBT = {D/(KB*T0):.1f})")

# ---- 4. nonequilibrium floor -----------------------------------------
noneq = {}
for T in (0.1, 0.05):
    f = 1.0 / (np.exp(RECIPES[0].Delta / (KB * T)) + 1.0)
    qs = np.logspace(-6, np.log10(0.5), 40)
    noneq[f"T{T}"] = dict(f=float(f), q=qs.tolist(),
                          penalty=noneq_penalty(f, qs).tolist())
res["noneq"] = noneq
for T in (0.1, 0.05):
    f = noneq[f"T{T}"]["f"]
    p4 = noneq_penalty(f, 1e-4)
    print(f"T={T}: f_th={f:.2e}, penalty at q=1e-4: x{p4:.1f}")
res["noneq_penalty_1e4_100mK"] = float(
    noneq_penalty(noneq["T0.1"]["f"], 1e-4))

# ---- 5. transparency inhomogeneity -----------------------------------
rng = np.random.default_rng(5)
inh = {}
for sd in (0.05, 0.1, 0.2):
    mean = 0.3
    var = sd ** 2
    a = mean * (mean * (1 - mean) / var - 1)
    b = (1 - mean) * (mean * (1 - mean) / var - 1)
    taus = rng.beta(a, b, 2000)
    # L-readout at phi -> 0: E = Delta (common), g ~ tau
    cv2 = np.var(taus) / np.mean(taus) ** 2
    deficit_L = np.sqrt(np.mean(taus ** 2) / np.mean(taus) ** 2)
    # I-readout at phi = 2: E and g both tau-dependent
    D0, T = RECIPES[0].Delta, 0.1
    phi = 2.0
    E = D0 * np.sqrt(1 - taus * np.sin(phi / 2) ** 2)
    g = D0 ** 2 * taus * np.sin(phi) / (4 * E)
    v = 1.0 / (2 * np.cosh(E / (2 * KB * T)) ** 2)
    achieved2 = np.sum(g ** 2 * v) * np.sum(v * E ** 2) / \
        np.sum(g * v * E) ** 2
    inh[str(sd)] = dict(deficit_L=float(deficit_L),
                        closed_form=float(np.sqrt(1 + cv2)),
                        deficit_I=float(np.sqrt(achieved2)))
    print(f"tau spread sd={sd}: deficit(L)={deficit_L:.4f} "
          f"[closed form {np.sqrt(1+cv2):.4f}], "
          f"deficit(I,phi=2)={np.sqrt(achieved2):.4f}")
res["tau_inhomogeneity"] = inh

# ---- 5b. fully diffusive (Dorokhov) transparency distribution --------
# tau = sech^2(x), x uniform in [0, X]: the Dorokhov distribution of a
# diffusive conductor. X chosen so <tau> matches the measured Ta/Ti/Au
# transparency 0.30. The bound itself is distribution independent; only
# the saturation deficit changes.
from scipy.optimize import brentq as _brentq
X = _brentq(lambda x: np.tanh(x) / x - 0.30, 1e-3, 50.0)
xs = np.linspace(0.0, X, 200001)[1:]
taus_d = 1.0 / np.cosh(xs) ** 2
m1, m2 = taus_d.mean(), (taus_d ** 2).mean()
dor = {}
dor["X"] = float(X)
dor["mean_tau"] = float(m1)
dor["deficit_L"] = float(np.sqrt(m2 / m1 ** 2))
D0, Td = RECIPES[0].Delta, 0.1
phi = 2.0
E_d = D0 * np.sqrt(1 - taus_d * np.sin(phi / 2) ** 2)
g_d = D0 ** 2 * taus_d * np.sin(phi) / (4 * E_d)
v_d = 1.0 / (2 * np.cosh(E_d / (2 * KB * Td)) ** 2)
ach2 = np.sum(g_d ** 2 * v_d) * np.sum(v_d * E_d ** 2) /     np.sum(g_d * v_d * E_d) ** 2
dor["deficit_I_phi2"] = float(np.sqrt(ach2))
res["dorokhov"] = dor
print(f"Dorokhov (X={X:.2f}, <tau>=0.30): deficit(L)="
      f"{dor['deficit_L']:.3f}, deficit(I,phi=2)="
      f"{dor['deficit_I_phi2']:.3f}")

# ---- 5c. time to measure tau_A from the spectral knee ----------------
# OU process (the occupation noise) plus white readout floor at the
# Ta/Ti/Au plateau-to-floor amplitude ratio; Welch periodogram, 3-para
# fit (plateau, knee, floor); fractional error on tau_A versus record
# length in units of tau_A.
from scipy.signal import lfilter, welch
from scipy.optimize import curve_fit
ratio_amp = 198.0 / 26.0                # sqrt(S) plateau / floor
rngk = np.random.default_rng(9)
tauA_u = 1.0                            # work in units of tau_A
dt = tauA_u / 20.0
a1 = np.exp(-dt / tauA_u)
sig_ou = 1.0
S_ou0 = 4 * sig_ou ** 2 * tauA_u        # S(0) of the OU line
S_fl = S_ou0 / ratio_amp ** 2
sig_w = np.sqrt(S_fl / (2 * dt))


def lor(f, A, tau, C):
    return A / (1 + (2 * np.pi * f * tau) ** 2) + C


knee = {}
for n_tau in (300, 1000, 3000, 10000):
    n = int(n_tau / dt)
    errs = []
    for _ in range(40):
        xi = rngk.standard_normal(n)
        y = lfilter([sig_ou * np.sqrt(1 - a1 ** 2)], [1, -a1], xi)
        y = y + sig_w * rngk.standard_normal(n)
        f_w, S_w = welch(y, fs=1 / dt, nperseg=min(n, 4096))
        f_w, S_w = f_w[1:], S_w[1:]
        try:
            popt, _ = curve_fit(lor, f_w, S_w,
                                p0=(S_ou0, 0.5 * tauA_u, S_fl),
                                maxfev=20000)
            errs.append(abs(popt[1] - tauA_u) / tauA_u)
        except Exception:
            errs.append(1.0)
    knee[str(n_tau)] = float(np.median(errs))
    print(f"record {n_tau:6d} tau_A: median |dtau/tau| = "
          f"{np.median(errs):.3f}")
res["knee_measurement"] = knee

# ---- 6. calibration-scale independence -------------------------------
sb1 = SensorBudget(RECIPES[0])
sj_cal = sb1.sj
sj_unc = ShortJunction(RECIPES[0])          # scale = 1
outs = {}
for name, sjx in (("calibrated", sj_cal), ("uncalibrated", sj_unc)):
    a, b, s = sjx.temperature_bound(0.1, 1e-6, 1.0, which="L")
    outs[name] = dict(sat_ratio=float(a / b),
                      CA_kB=float(s["C_A"] / KB))
res["scale_independence"] = outs
res["scale_independent_ok"] = bool(
    abs(outs["calibrated"]["sat_ratio"] -
        outs["uncalibrated"]["sat_ratio"]) < 1e-9 and
    abs(outs["calibrated"]["CA_kB"] -
        outs["uncalibrated"]["CA_kB"]) < 1e-9)
print("scale independence of saturation ratio and C_A:",
      res["scale_independent_ok"])

json.dump(res, open(os.path.join(OUT, "limits.json"), "w"))
print("limits.json written")
