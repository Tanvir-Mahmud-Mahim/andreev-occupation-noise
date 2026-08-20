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


if __name__ == "__main__":
    test_short_junction_formula()
    test_kulik_levels()
    test_ballistic_IcRn()
    test_recipe_sanity()
    print("ALL ABS TESTS PASSED")
