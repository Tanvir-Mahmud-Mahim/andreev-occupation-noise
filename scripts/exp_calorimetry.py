"""Fig. 4 data: matched-filter energy resolution versus tauA (numeric
matched filter with readout floor against the analytic Andreev-only
form), a simulated single-photon detection trace, and detection
statistics for a 26 GHz photon.
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from constants import KB, H_PLANCK
from materials import Recipe, RECIPES
from sensor_limits import SensorBudget
from scipy.optimize import brentq
from abs_model import gap_bcs

OUT = os.path.join(os.path.dirname(__file__), "..", "data")
res = {}

T0 = 0.1
g = lambda Tc: gap_bcs(T0, Tc, 1.7639 * KB * Tc) - 2.3994 * KB * T0
Tc_opt = brentq(g, T0 * 1.01, T0 * 6)
matched = Recipe("matched", "matched", Tc_opt, 1e-6, 0.1e-6, 1.0e-6,
                 30.0, 0.3, RECIPES[0].Ic20, RECIPES[0].Rn)
sb = SensorBudget(matched)

# readout floor: quantum-limited phase readout, n = 30, kappa/2pi = 1 MHz
kappa = 2 * np.pi * 1e6
nbar = 30.0
S_omega = (kappa / 4.0)**2 * 2.0 / (nbar * kappa)
nu_r = sb.nu_r
S_ro_y = S_omega / (2 * np.pi * nu_r)**2

tauAs = np.logspace(-8.5, -4.5, 33)
num, ana, num_nofloor = [], [], []
for tA in tauAs:
    num.append(sb.energy_resolution(T0, tA, S_ro_y=S_ro_y) /
               H_PLANCK / 1e9)
    num_nofloor.append(sb.energy_resolution(T0, tA, S_ro_y=0.0) /
                       H_PLANCK / 1e9)
    ana.append(sb.energy_resolution_analytic_A(T0, tA) / H_PLANCK / 1e9)
res["tauAs"] = tauAs.tolist()
res["sigE_num"] = num
res["sigE_num_nofloor"] = num_nofloor
res["sigE_ana"] = ana
res["S_ro_y"] = S_ro_y
res["Tc_opt"] = Tc_opt
res["T0"] = T0

# anchor: numeric (no floor) vs analytic in Andreev-dominated regime
i = int(np.argmin(np.abs(tauAs - 1e-6)))
res["anchor_ratio"] = num_nofloor[i] / ana[i]
print(f"matched design at 100 mK: Tc*={Tc_opt:.3f} K")
print(f"anchor tauA=1us: numeric(no floor)/analytic = "
      f"{res['anchor_ratio']:.3f}")
print(f"sigE (floor incl.): {num[i]:.2f} GHz; analytic {ana[i]:.2f} GHz")

# ---- simulated detection trace ---------------------------------------
rng = np.random.default_rng(3)
tauA = 1e-6
tth = sb.tau_th(T0)
Ce = sb.Ce(T0)
fs = 50e6                      # 50 MS/s
n = 60000
t = np.arange(n) / fs
Ephot = H_PLANCK * 26e9
dT_pulse = np.zeros(n)
i0 = n // 3
dT_pulse[i0:] = (Ephot / Ce) * np.exp(-(t[i0:] - t[i0]) / tth)
# frequency signal: through occupation lag
s = sb.sj.andreev_sums(1e-4, T0, "L")
I1 = sb.sj.dIdphi0(T0)
conv = (sb.participation(T0) / 2.0) / I1
R = s["R_occ"] * conv
HA = np.exp(-np.arange(0, 12 * tauA, 1 / fs) / tauA)
HA /= HA.sum()
y_sig = R * np.convolve(dT_pulse, HA, mode="full")[:n]
# noises: Andreev Lorentzian (Gaussian approx of many channels),
# phonon TFN (Lorentzian, tau_th), readout white
S_A0 = conv**2 * s["S0_over_tau"] * tauA
xw = rng.standard_normal(n) * np.sqrt(S_A0 * fs / 2)
a = np.exp(-1.0 / (fs * tauA))
yA = np.zeros(n)
for k in range(1, n):
    yA[k] = a * yA[k - 1] + np.sqrt(1 - a**2) * xw[k]
S_ph0 = R**2 * 4 * KB * T0**2 / sb.Gep(T0)
xw2 = rng.standard_normal(n) * np.sqrt(S_ph0 * fs / 2)
b = np.exp(-1.0 / (fs * tth))
yP = np.zeros(n)
for k in range(1, n):
    yP[k] = b * yP[k - 1] + np.sqrt(1 - b**2) * xw2[k]
yR = rng.standard_normal(n) * np.sqrt(S_ro_y * fs / 2)
y_tot = y_sig + yA + yP + yR
res["trace"] = dict(t_us=(t * 1e6).tolist(),
                    y_nu=(y_tot * nu_r).tolist(),
                    y_sig_nu=(y_sig * nu_r).tolist())

# detection statistics: matched-filter SNR for the 26 GHz photon
sigE = sb.energy_resolution(T0, tauA, S_ro_y=S_ro_y)
snr = Ephot / sigE
# dark counts: Gaussian threshold at Ephot/2, effective rate ~ bandwidth
nu_eff = 1.0 / (2 * np.pi * max(tauA, tth))
dark = nu_eff * np.exp(-0.5 * (0.5 * snr)**2)
res["snr26"] = snr
res["dark_rate_hz"] = dark
print(f"26 GHz photon: SNR = {snr:.1f}, dark rate at E/2 threshold = "
      f"{dark:.2e} /s")

with open(os.path.join(OUT, "calorimetry.json"), "w") as f:
    json.dump(res, f)
