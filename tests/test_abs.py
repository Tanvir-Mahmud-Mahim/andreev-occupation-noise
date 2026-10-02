"""Testbench: ABS solver analytic anchors."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from constants import HBAR, KB, E_CHARGE
from abs_model import abs_energies, JunctionModel
from materials import RECIPES

Delta = 1.5e-23  # J


def test_short_junction_formula():
    """L->0: E = Delta sqrt(1 - tau sin^2(phi/2))."""
    worst = 0.0
    for tau in (0.1, 0.3, 0.5, 0.78, 0.99, 1.0):
        for phi in (0.3, 1.0, 2.0, 3.0):
            E = abs_energies(phi, tau, 0.0, Delta)
            Ean = Delta * np.sqrt(1 - tau * np.sin(phi / 2) ** 2)
            assert len(E) == 1
            worst = max(worst, abs(E[0] - Ean) / Delta)
    print(f"short-junction anchor: max |dE|/Delta = {worst:.2e}")
    assert worst < 1e-12


def test_kulik_levels():
    """tau=1, finite L: 2 arccos(E/Delta) - cE = +-phi mod 2pi."""
    c = 3.0 / Delta   # eta(Delta) = 3 rad, a moderately long junction
    worst = 0.0
    for phi in (0.5, 1.5, 2.5):
        for E in abs_energies(phi, 1.0, c, Delta):
            th = 2 * np.arccos(E / Delta) - c * E
            r = min(abs(((th - phi) + np.pi) % (2 * np.pi) - np.pi),
                    abs(((th + phi) + np.pi) % (2 * np.pi) - np.pi))
            worst = max(worst, r)
    print(f"Kulik anchor: max residual = {worst:.2e} rad")
    assert worst < 1e-10


def test_ballistic_IcRn():
    """Single ballistic channel, short junction, T->0:
    e Ic Rn = pi Delta (Rn = h/2e^2 spin-degenerate)."""
    from materials import Recipe
    # fabricate an artificial 1-mode recipe: narrow W, tau=1, tiny L
    r = Recipe("test", "test", 1.0, 1e-5, 1e-9, 12e-9, 30.0, 1.0, 1e-6, 1.0)
    m = JunctionModel(r, n_phi=601)
    assert m.N >= 1
    m.cos_th = m.cos_th[:1]; m.cs = m.cs[:1]; m.N = 1
    m._levels = None
    Ic = m.Ic(0.001)  # ~T -> 0 ; valley factor 2 included in model
    # model current includes spin (in F) and valley (x2): one orbital mode
    # -> Ic = 2 * e Delta / hbar (two spin-degenerate channels)
    Ican = 2.0 * E_CHARGE * r.Delta / HBAR
    err = abs(Ic - Ican) / Ican
    print(f"ballistic IcRn anchor: rel err = {err:.2e}")
    assert err < 5e-3   # limited by phase-grid resolution near phi=pi


def test_recipe_sanity():
    """Calibration factors should be O(1) for all recipes."""
    for r in RECIPES[:2]:
        m = JunctionModel(r, n_phi=41)
        s = m.calibrate()
        print(f"{r.label}: N={m.N}, scale={s:.3f}")
        assert 0.1 < s < 10.0


def test_free_energy_continuity_at_exit():
    """When a bound level leaves the gap as phi varies, the bound-state
    free energy jumps by about 3 kB T (here Delta/kBT = 3) and the
    continuum term (scattering phase plus threshold boundary term) must
    jump by the opposite amount. Checked for several transparencies and
    lengths on a 4001-point phase grid; the residual must be below 1% of
    the jump."""
    from abs_model import continuum_delta, continuum_phase_grid
    T = Delta / (3.0 * KB)
    x_D = Delta / (2 * KB * T)
    l2c = 2 * KB * T * (x_D + np.log1p(np.exp(-2 * x_D)))
    Efac = continuum_phase_grid(20000)
    E = Delta * Efac
    w = np.tanh(E / (2 * KB * T))

    def parts(phi, tau, c):
        Eb = abs_energies(phi, tau, c, Delta)
        x = np.abs(Eb) / (2 * KB * T)
        Fb = -np.sum(2 * KB * T * (x + np.log1p(np.exp(-2 * x))))
        d = continuum_delta(phi, tau, c, Delta, Efac)
        Fc = l2c * d[0] / np.pi + np.trapezoid(w * d, E) / np.pi
        return Fb, Fc

    worst, n_exits = 0.0, 0
    phis = np.linspace(1e-3, np.pi - 1e-3, 4001)
    for tau, cD in ((0.5, 1.0), (0.9, 1.0), (0.9, 2.0), (0.9, 4.0)):
        c = cD / Delta
        n = [len(abs_energies(p, tau, c, Delta)) for p in phis]
        for i in range(1, len(n)):
            if n[i] != n[i - 1]:
                a, b = parts(phis[i - 1], tau, c), parts(phis[i], tau, c)
                jump = abs(b[0] - a[0])
                resid = abs((b[0] + b[1]) - (a[0] + a[1]))
                worst = max(worst, resid / jump)
                n_exits += 1
                break
    print(f"free-energy continuity at bound-state exit: {n_exits} exits, "
          f"max residual/jump = {worst:.4f}")
    assert n_exits >= 4 and worst < 0.01


def test_bcs_gap_asymptotes():
    """Tabulated BCS gap against its two asymptotes: low-T form
    1 - u = sqrt(2 pi t/A) exp(-A/t) at t = 0.15 (within 1% in 1 - u)
    and the Ginzburg-Landau form u = 1.74 sqrt(1 - t) at t = 0.995
    (within 0.5%)."""
    from abs_model import gap_bcs
    A = 1.764
    t = 0.15
    one_minus_u = 1.0 - gap_bcs(t, 1.0, 1.0)
    asym = np.sqrt(2 * np.pi * t / A) * np.exp(-A / t)
    e_low = abs(one_minus_u - asym) / asym
    t = 0.995
    e_gl = abs(gap_bcs(t, 1.0, 1.0) - 1.74 * np.sqrt(1 - t)) / \
        (1.74 * np.sqrt(1 - t))
    print(f"BCS gap asymptotes: low-T {e_low:.4f}, Ginzburg-Landau {e_gl:.4f}")
    assert e_low < 0.01 and e_gl < 0.005


if __name__ == "__main__":
    test_short_junction_formula()
    test_kulik_levels()
    test_ballistic_IcRn()
    test_recipe_sanity()
    test_free_energy_continuity_at_exit()
    test_bcs_gap_asymptotes()
    print("ALL ABS TESTS PASSED")
