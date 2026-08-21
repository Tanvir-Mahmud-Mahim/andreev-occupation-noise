"""Testbench: generalized (pair-process) noise model and anchors."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from noise_general import channel_generator, sigma_spectrum


def test_singles_limit():
    """Gp = 0 must reproduce the independent-spin Lorentzian:
    var = 2f(1-f), tau_eff = 1/Gs."""
    for f in (0.05, 0.3, 0.5):
        S, S0, var, te = sigma_spectrum(f, Gs=2.0, Gp=0.0,
                                        omegas=np.array([0.0]))
        assert abs(var - 2 * f * (1 - f)) < 1e-12
        assert abs(te - 0.5) < 1e-12
    print("singles limit: var and tau_eff exact")


def test_equilibrium_invariance():
    """Equilibrium variance independent of pair rate."""
    for Gp in (0.0, 0.5, 3.0, 30.0):
        _, _, var, _ = sigma_spectrum(0.2, 1.0, Gp, np.array([0.0]))
        assert abs(var - 2 * 0.2 * 0.8) < 1e-12
    print("equilibrium variance invariant under pair processes")


def test_pairs_only_shorten():
    """tau_eff decreases monotonically with the pair rate, so the
    single-quasiparticle bound is the worst case."""
    taus = []
    for Gp in (0.0, 0.3, 1.0, 3.0, 10.0):
        _, _, _, te = sigma_spectrum(0.15, 1.0, Gp, np.array([0.0]))
        taus.append(te)
    d = np.diff(taus)
    print("tau_eff vs Gp:", [round(t, 4) for t in taus])
    assert np.all(d <= 1e-12)


def test_generator_conservation():
    W = channel_generator(0.23, 1.7, 0.6)
    assert np.allclose(W.sum(axis=0), 0.0, atol=1e-14)
    pi = np.array([(1 - .23) ** 2, .23 * .77, .23 * .77, .23 ** 2])
    assert np.allclose(W @ pi, 0.0, atol=1e-14)
    print("generator conserves probability; Gibbs state stationary")


def test_mc_pair_spectrum():
    """Monte Carlo of the 4-state channel against the analytic S(0)."""
    rng = np.random.default_rng(11)
    f, Gs, Gp, dt, n, M = 0.3, 1.0, 2.0, 0.02, 300_000, 16
    W = channel_generator(f, Gs, Gp)
    P = np.eye(4) + W * dt          # transition kernel (dt small)
    assert np.diag(P).min() > 0.9 and P.min() >= 0
    sig = np.array([1.0, 0.0, 0.0, -1.0])
    state = rng.integers(0, 4, size=M)
    acc = np.zeros(n)
    cum = np.cumsum(P, axis=0)      # cum[i,j] = P(state<=i | from j)
    for k in range(n):
        r = rng.random(M)
        state = (r[None, :] > cum[:, state]).sum(axis=0)
        acc[k] = sig[state].sum()
    # zero-frequency density from window-averages: var(mean_t) = S0/(2t)
    t_int = 40.0
    kk = int(t_int / dt)
    nw = n // kk
    means = acc[:nw * kk].reshape(nw, kk).mean(axis=1)
    S0_mc = 2 * t_int * np.var(means)
    _, S0_an, _, _ = sigma_spectrum(f, Gs, Gp, np.array([0.0]))
    S0_an *= M
    err = abs(S0_mc - S0_an) / S0_an
    print(f"pair-process MC: S0={S0_mc:.3f} vs analytic {S0_an:.3f} "
          f"({err:.1%})")
    assert err < 0.25


def test_bound_random_ensembles():
    """The Cauchy-Schwarz bound dT^2 >= 2 kB T^2 tauA / (C_A t) holds
    for arbitrary couplings over random level sets, with exact equality
    at g proportional to E (spin factor included)."""
    from constants import KB
    rng = np.random.default_rng(3)
    for _ in range(20):
        M = int(rng.integers(2, 60))
        T = float(rng.uniform(0.03, 0.4))
        E = rng.uniform(0.2, 6.0, M) * KB * T
        g = rng.uniform(-2, 2, M)
        tauA, t = 1e-6, 1.0
        f = 1.0 / (np.exp(E / (KB * T)) + 1.0)
        v = f * (1 - f)
        R = np.sum(g * (-2 * E * v / (KB * T ** 2)))
        dT2 = np.sum(g ** 2 * 8 * v * tauA) / (2 * t * R ** 2)
        CA = np.sum(2 * v * E ** 2) / (KB * T ** 2)
        bound = 2 * KB * T ** 2 * tauA / (CA * t)
        assert dT2 >= bound * (1 - 1e-12)
        Ro = np.sum(E * (-2 * E * v / (KB * T ** 2)))
        dT2o = np.sum(E ** 2 * 8 * v * tauA) / (2 * t * Ro ** 2)
        assert abs(dT2o / bound - 1) < 1e-12
    print("random-ensemble bound: inequality and g~E equality exact")


if __name__ == "__main__":
    test_singles_limit()
    test_equilibrium_invariance()
    test_pairs_only_shorten()
    test_generator_conservation()
    test_mc_pair_spectrum()
    test_bound_random_ensembles()
    print("ALL GENERAL-NOISE TESTS PASSED")
