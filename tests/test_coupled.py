"""Testbench for the coupled two-node model and the click simulation
(revision, 2026). Run: python3 tests/test_coupled.py"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from scipy.optimize import brentq
from scipy.linalg import solve_continuous_lyapunov
from constants import KB, H_PLANCK
from materials import Recipe, RECIPES
from sensor_limits import SensorBudget
from abs_model import gap_bcs
import coupled as cp
from click_sim import ClickModel

OUT = []


def check(name, val, tol):
    ok = abs(val) <= tol
    OUT.append(f"{'PASS' if ok else 'FAIL'}  {name}: {val:.3e} (tol {tol})")
    assert ok, name


def matched(T0, W, L):
    g = lambda Tc: gap_bcs(T0, Tc, 1.7639 * KB * Tc) - 2.3994 * KB * T0
    Tc = brentq(g, T0 * 1.01, T0 * 6)
    return Recipe("m", "m", Tc, 1e-6, L * 1e-6, W * 1e-6, 30.0, 0.3,
                  RECIPES[0].Ic20, RECIPES[0].Rn)


# 1. closed form = numerical optimal filter (no readout floor)
for r_, T in ((RECIPES[0], 0.1), (matched(0.1, 1.0, 0.1), 0.1),
              (matched(0.05, 5.3, 0.5), 0.05)):
    sb = SensorBudget(r_)
    s = sb.sj.andreev_sums(1e-4, T, "L")
    for tA in (1e-8, 1e-7, 1e-6):
        a = cp.sigmaE_closed(T, sb.Ce(T), s["C_A"], sb.Gep(T), tA)
        b = cp.sigmaE_numeric(T, sb.Ce(T), s["C_A"], sb.Gep(T), tA)
        check(f"closed vs numeric sigma_E, {r_.label} T={T} tauA={tA:.0e}",
              b / a - 1, 2e-3)

# 2. earlier one-way numerical filter (fixed integration range) agrees
sb = SensorBudget(matched(0.1, 1.0, 0.1))
s = sb.sj.andreev_sums(1e-4, 0.1, "L")
for tA in (1e-7, 1e-6):
    a = cp.sigmaE_closed(0.1, sb.Ce(0.1), s["C_A"], sb.Gep(0.1), tA)
    b = sb.energy_resolution(0.1, tA)
    check(f"one-way numeric = two-node closed form, tauA={tA:.0e}",
          b / a - 1, 5e-3)

# 3. linear-noise model of the click simulation: equilibrium covariance
#    and the linear-response energy resolution with and without the
#    energy coupling equal the closed form
r_ = matched(0.05, 5.3, 0.5)
sb = SensorBudget(r_, nu_r=3e9)
sj = sb.sj
cm = ClickModel(0.05, sj.Delta(0.05), sj.Nch, sb.Ce(0.05) / 0.05,
                sb.sigma * r_.area, sb.delta, 3e9, sb.L_r,
                sj.dIdphi0(0.05))
s = sj.andreev_sums(1e-4, 0.05, "L")
for tA in (3e-8, 1e-6):
    for cE in (0.0, 1.0):
        J = cm.jac(cm.U0, cm.f0, "C", tA, cE)
        Q = cm.Qmat(cm.U0, cm.f0, "C", tA, cE)
        S = solve_continuous_lyapunov(J, -Q)
        check(f"var(U_e) = kB T^2 C_e, cE={cE}, tauA={tA:.0e}",
              S[0, 0] * cm.U0 ** 2 / (KB * 0.05 ** 2 * sb.Ce(0.05)) - 1,
              1e-6)
        if cE == 1.0:
            check(f"var(n) = f(1-f)/N_s (coupled), tauA={tA:.0e}",
                  S[1, 1] / (cm.f0 * (1 - cm.f0) / cm.Ns) - 1, 1e-6)
        f = np.logspace(3, 11, 20000)
        w = 2 * np.pi * f
        I2 = np.eye(2)
        H = np.array([np.linalg.solve(1j * x * I2 - J,
                                      np.array([1.0 / cm.U0, 0.0]))[1]
                      for x in w])
        Sn = np.array([2 * (np.linalg.solve(1j * x * I2 - J, I2) @ Q @
                            np.linalg.solve(1j * x * I2 - J, I2).conj().T
                            )[1, 1].real for x in w])
        sig = 1 / np.sqrt(4 * np.trapezoid(np.abs(H) ** 2 / Sn, f))
        ref = cp.sigmaE_closed(0.05, sb.Ce(0.05), s["C_A"], sb.Gep(0.05),
                               tA)
        check(f"LNA linear sigma_E = closed form, cE={cE}, tauA={tA:.0e}",
              sig / ref - 1, 3e-3)

# 4. corner time of the coupled noise spectrum ~ tauA + C_A/G when
#    C_e << C_A, and ~ tauA for the measured recipe
sb = SensorBudget(RECIPES[0])
s = sb.sj.andreev_sums(1e-4, 0.1, "L")
tc = cp.corner_time(0.1, sb.Ce(0.1), s["C_A"], sb.Gep(0.1), 1e-6)
check("corner time = tauA for Ta/Ti/Au (tauA = 1 us)", tc / 1e-6 - 1,
      5e-3)

for line in OUT:
    print(line)
print("ALL COUPLED-MODEL TESTS PASSED")
