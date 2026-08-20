"""Testbench: occupation-noise model against Monte Carlo and analytics."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from constants import KB, HBAR, E_CHARGE
from materials import RECIPES
from short_junction import ShortJunction
from montecarlo import telegraph_traces, psd_single_sided


def test_saturation_uniform_tau():
    """Uniform-tau short junction saturates the Cauchy-Schwarz bound
    for both current and inductance readout."""
    for r in RECIPES[:3]:
        sj = ShortJunction(r)
        sj.calibrate()
        T = 0.4 * r.Tc
        for which in ("I", "L"):
            a, b, _ = sj.temperature_bound(T, tauA=1e-6, t_int=1e-3,
                                           which=which)
            ratio = a / b
            print(f"{r.label:12s} {which}: achieved/bound = {ratio:.6f}")
            assert abs(ratio - 1.0) < 1e-9


def test_responsivity_occupation_channel():
    """Analytic R_occ equals numerical dO/dT with the gap frozen."""
    r = RECIPES[1]
    sj = ShortJunction(r)
    sj.calibrate()
    T = 0.35 * r.Tc
    phi = sj.phi_max(T)
    s = sj.andreev_sums(phi, T, "I")
    # numerical: vary occupation temperature only, keep Delta(T) fixed
    D = sj.Delta(T)
    E = sj.E(phi, T)
    g = sj.scale * (2 * E_CHARGE / HBAR) * (-sj.dEdphi(phi, T))
    dT = 1e-6
    occ = lambda TT: np.tanh(E / (2 * KB * TT))
    R_num = sj.Nch * g * (occ(T + dT) - occ(T - dT)) / (2 * dT)
    err = abs(R_num - s["R_occ"]) / abs(R_num)
    print(f"R_occ analytic vs numeric: rel err = {err:.2e}")
    assert err < 1e-5


def test_telegraph_psd_and_variance():
    """MC telegraph noise reproduces the Lorentzian PSD and the
    var = S(0)/(2t) estimator convention."""
    rng = np.random.default_rng(7)
    f, tauA = 0.3, 1.0
    dt, n_steps, M = 0.05, 400_000, 24
    tr = telegraph_traces(f, tauA, dt, n_steps, M, rng)
    x = np.sum(1.0 - 2.0 * tr, axis=1)          # ~ sum of (1-2n)
    freqs, S = psd_single_sided(x, dt)
    S0_an = M * 4.0 * (2**2) * f * (1 - f) * tauA   # (1-2n): g=2? no:
    # var(1-2n) = 4 f(1-f); telegraph S(0) = 4 var tau = 16 f(1-f) tau
    S0_an = M * 4.0 * 4.0 * f * (1 - f) * tauA
    m = (freqs > 0.005) & (freqs < 0.03)        # below the knee
    S0_mc = np.mean(S[m])
    err = abs(S0_mc - S0_an) / S0_an
    print(f"telegraph S(0): MC={S0_mc:.1f} analytic={S0_an:.1f} "
          f"err={err:.2%}")
    assert err < 0.10
    # knee: S(1/(2 pi tauA)) should be S0/2
    fk = 1.0 / (2 * np.pi * tauA)
    mk = (freqs > 0.8 * fk) & (freqs < 1.25 * fk)
    ratio = np.mean(S[mk]) / S0_an
    print(f"knee ratio (expect ~0.5): {ratio:.3f}")
    assert abs(ratio - 0.5) < 0.12
    # variance of time-averages over windows of length t_int
    t_int = 50.0
    k = int(t_int / dt)
    nwin = (n_steps // k)
    means = x[:nwin * k].reshape(nwin, k).mean(axis=1)
    var_mc = np.var(means)
    var_an = S0_an / (2 * t_int)
    err2 = abs(var_mc - var_an) / var_an
    print(f"window-average variance: MC={var_mc:.3f} "
          f"analytic={var_an:.3f} err={err2:.2%}")
    assert err2 < 0.25


if __name__ == "__main__":
    test_saturation_uniform_tau()
    test_responsivity_occupation_channel()
    test_telegraph_psd_and_variance()
    print("ALL NOISE TESTS PASSED")
