"""Self-consistent single-photon click simulation (revision, 2026).

Two coupled variables describe each trial: the electron energy
U_e = gamma T_e^2 / 2 and the mean occupation n of the N_s = 2 N_ch
spin-resolved Andreev states (all at the energy E = Delta*, the gap
edge, as for inductance readout near zero phase). The readout sees
m = <sigma> = 1 - 2 n.

Mean (deterministic, fully nonlinear) dynamics after a deposit E_gamma
into the electrons:
    dU_e/dt = -Sigma A (T_e^d - T_0^d) - c_E N_s E (f(T_e) - n)/tau(T_e)
    dn/dt   = (f(T_e) - n)/tau(T_e)
with c_E = 1 when the energy carried by the exchanged quasiparticles is
returned to the electrons (energy coupling, the Andreev heat capacity
C_A then enters the heat balance) and c_E = 0 for the earlier one-way
model in which the occupations follow T_e passively. tau(T_e) is either
constant (scenario C) or tau_A exp[-(Delta/kB)(1/T_0 - 1/T_e)], floored
at 1 ns (scenario T).

Fluctuations are treated in the linear-noise (van Kampen) approximation
around the mean trajectory: x = (dU_e, dn) obeys dx = J(t) x dt + noise,
with J the Jacobian of the mean dynamics and the noise covariance rate
    Q = [[S_P/2 + c_E (N_s E)^2 D_n, -c_E N_s E D_n],
         [-c_E N_s E D_n,            D_n]],
    S_P = 2 kB d Sigma A (T_e^(d+1) + T_0^(d+1))  (single-sided),
    D_n = [f (1 - n) + (1 - f) n] / (tau N_s).
With `nonstationary=True`, J and Q are evaluated along the hot
trajectory (instantaneous rates and variances); with False they keep
their bath values (earlier model). The linear SDE is integrated exactly
over each step (Van Loan matrix exponential), which is stable for the
stiff case C_A >> C_e. At equilibrium the stationary covariance gives
var(n) = f(1-f)/N_s and var(U_e) = kB T^2 C_e (tests/test_coupled.py).

Readout: (i) an ideal instantaneous frequency meter with a white floor,
or (ii) a one-port resonator of linewidth kappa driven at its unperturbed
frequency, integrated exactly step by step for the time-dependent
detuning; the output phase carries white noise of single-sided density
S_theta = 2 (1 + 2 n_add) / (n_bar kappa) (n_add = 0: quantum-limited
phase-sensitive readout; n_add = kB T_N / h nu for an amplifier of noise
temperature T_N).
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm, solve_continuous_lyapunov

from constants import KB, HBAR, E_CHARGE


class ClickModel:
    def __init__(self, T0, Delta, Nch, gamma, SigmaA, delta, nu0, Lr,
                 I1_0):
        self.T0, self.E, self.Nch = T0, Delta, Nch
        self.Ns = 2 * Nch
        self.gamma, self.SA, self.d = gamma, SigmaA, delta
        self.nu0, self.Lr = nu0, Lr
        self.f0 = float(self.fermi(T0))
        self.m0 = 1.0 - 2.0 * self.f0
        self.U0 = 0.5 * gamma * T0 ** 2
        self.I1_0 = I1_0
        LJ0 = (HBAR / (2 * E_CHARGE)) / I1_0
        self.Cres = 1.0 / ((2 * np.pi * nu0) ** 2 * (Lr + LJ0))

    # ---- helpers ---------------------------------------------------------
    def fermi(self, T):
        x = np.clip(self.E / (KB * np.maximum(T, 1e-6)), 0, 700)
        return 1.0 / (np.exp(x) + 1.0)

    def Te_of_U(self, U):
        return np.sqrt(2.0 * np.maximum(U, 1e-40) / self.gamma)

    def nu_of_m(self, m):
        I1 = self.I1_0 * np.maximum(m, 1e-6) / self.m0
        LJ = (HBAR / (2 * E_CHARGE)) / I1
        return 1.0 / (2 * np.pi * np.sqrt((self.Lr + LJ) * self.Cres))

    def tau_of(self, Te, scen, tauA):
        if scen == "C":
            return tauA
        expo = -(self.E / KB) * (1.0 / self.T0 - 1.0 / max(Te, 1e-6))
        return max(1e-9, tauA * float(np.exp(min(expo, 0.0))))

    # ---- mean dynamics ---------------------------------------------------
    def drift(self, U, n, scen, tauA, cE):
        Te = float(self.Te_of_U(U))
        f = float(self.fermi(Te))
        tau = self.tau_of(Te, scen, tauA)
        P = self.SA * (Te ** self.d - self.T0 ** self.d)
        aU = -P - cE * self.Ns * self.E * (f - n) / tau
        an = (f - n) / tau
        return aU, an

    def jac(self, U, n, scen, tauA, cE):
        """Jacobian of the mean dynamics (central differences), in the
        scaled variables (U_e / U_0, n) and units of 1/s."""
        hU = 1e-6 * max(U, self.U0)
        hn = 1e-6 * max(n, 1e-3)
        a1 = np.array(self.drift(U + hU, n, scen, tauA, cE))
        a2 = np.array(self.drift(U - hU, n, scen, tauA, cE))
        b1 = np.array(self.drift(U, n + hn, scen, tauA, cE))
        b2 = np.array(self.drift(U, n - hn, scen, tauA, cE))
        J = np.empty((2, 2))
        J[:, 0] = (a1 - a2) / (2 * hU)
        J[:, 1] = (b1 - b2) / (2 * hn)
        Sc = np.diag([1.0 / self.U0, 1.0])
        return Sc @ J @ np.linalg.inv(Sc)

    def Qmat(self, U, n, scen, tauA, cE):
        """Noise covariance rate in the scaled variables."""
        Sc = np.diag([1.0 / self.U0, 1.0])
        return Sc @ self._Qraw(U, n, scen, tauA, cE) @ Sc

    def _Qraw(self, U, n, scen, tauA, cE):
        Te = float(self.Te_of_U(U))
        f = float(self.fermi(Te))
        tau = self.tau_of(Te, scen, tauA)
        T0, d = self.T0, self.d
        SP = 2 * KB * d * self.SA * (Te ** (d + 1) + T0 ** (d + 1))
        Dn = (f * (1 - n) + (1 - f) * n) / (tau * self.Ns)
        k = self.Ns * self.E
        return np.array([[0.5 * SP + cE * k * k * Dn, -cE * k * Dn],
                         [-cE * k * Dn, Dn]])

    def mean_trajectory(self, tgrid, i0, Egamma, scen, tauA, cE):
        """Noiseless (U_e, n) on the record grid; deposit at index i0."""
        U = np.full(len(tgrid), self.U0)
        n = np.full(len(tgrid), self.f0)
        if Egamma <= 0:
            return U, n

        def rhs(t, y):
            return self.drift(max(y[0], 1e-40), y[1], scen, tauA, cE)
        t0 = tgrid[i0]
        sol = solve_ivp(rhs, (t0, tgrid[-1]),
                        [self.U0 + Egamma, self.f0], method="Radau",
                        t_eval=tgrid[i0:], rtol=1e-8,
                        atol=[1e-10 * self.U0, 1e-12])
        U[i0:], n[i0:] = sol.y[0], sol.y[1]
        return U, n

    # ---- stochastic records ----------------------------------------------
    def records(self, ntr, tgrid, i0, Egamma, scen, tauA, rng,
                coupled=True, nonstationary=True):
        """Occupation records m(t) = 1 - 2 n(t) for ntr trials and the
        noiseless template, in the linear-noise approximation."""
        cE = 1.0 if coupled else 0.0
        dt = tgrid[1] - tgrid[0]
        Ubar, nbar = self.mean_trajectory(tgrid, i0, Egamma, scen, tauA,
                                          cE)
        J0 = self.jac(self.U0, self.f0, scen, tauA, cE)
        Q0 = self.Qmat(self.U0, self.f0, scen, tauA, cE)
        S0 = solve_continuous_lyapunov(J0, -Q0)
        L0 = np.linalg.cholesky(S0)
        x = (L0 @ rng.standard_normal((2, ntr)))
        Phi0, C0 = self._step(J0, Q0, dt)
        out = np.empty((ntr, len(tgrid)))
        for k in range(len(tgrid)):
            out[:, k] = 1.0 - 2.0 * (nbar[k] + x[1])
            if k == len(tgrid) - 1:
                break
            if nonstationary and k >= i0 and Egamma > 0:
                Um = 0.5 * (Ubar[k] + Ubar[k + 1])
                nm = 0.5 * (nbar[k] + nbar[k + 1])
                J = self.jac(Um, nm, scen, tauA, cE)
                Q = self.Qmat(Um, nm, scen, tauA, cE)
                Phi, C = self._step(J, Q, dt)
            elif k >= i0 and Egamma > 0:
                # earlier model: rates follow the trajectory in the mean
                # only; the fluctuations keep their bath statistics
                Phi, C = Phi0, C0
            else:
                Phi, C = Phi0, C0
            x = Phi @ x + C @ rng.standard_normal((2, ntr))
        return out, 1.0 - 2.0 * nbar, Ubar

    @staticmethod
    def _step(J, Q, h):
        """Exact discretization of dx = J x dt + noise(Q) over h."""
        M = np.zeros((4, 4))
        M[:2, :2] = -J
        M[:2, 2:] = Q
        M[2:, 2:] = J.T
        E = expm(M * h)
        Phi = E[2:, 2:].T
        Sig = Phi @ E[:2, 2:]
        Sig = 0.5 * (Sig + Sig.T)
        w, V = np.linalg.eigh(Sig)
        C = V @ np.diag(np.sqrt(np.maximum(w, 0.0)))
        return Phi, C


def resonator_phase(dnu, dt, kappa):
    """Output phase of a one-port resonator (linewidth kappa, rad/s),
    driven at its unperturbed frequency, when its frequency is shifted by
    dnu(t) (Hz). Exact step integration of
        d alpha/dt = -(kappa/2 + i Delta) alpha + kappa,
    alpha = 2 in the unperturbed steady state; the output field is
    proportional to alpha - 1 (unit modulus for a lossless one-port).
    Returns the unwrapped phase (rad), zero when unperturbed. For a slow,
    small detuning the phase is -4 Delta / kappa."""
    ntr, nrec = dnu.shape
    Dl = 2 * np.pi * dnu
    a = np.full(ntr, 2.0 + 0j)
    th = np.empty((ntr, nrec))
    for k in range(nrec):
        lam = kappa / 2 + 1j * Dl[:, k]
        ass = kappa / lam
        a = ass + (a - ass) * np.exp(-lam * dt)
        th[:, k] = np.angle(a - 1.0)
    return np.unwrap(th, axis=1)


def matched_filter_stats(Y1, Y0, Ypsd, tpl, thr=0.5):
    """Optimal (whitened) matched filter with an empirical noise PSD.

    Y1, Y0: photon and dark records; Ypsd: independent dark records used
    only to estimate the noise PSD; tpl: noiseless template. The score is
    normalized so that the template scores 1. Returns a dict with the
    Kumar-Carroll index d_KC = (mu1 - mu0)/sqrt((s0^2 + s1^2)/2), which is
    our Monte Carlo signal-to-noise ratio, the efficiency and dark
    fraction at the threshold, and the 1/e correlation lag of the running
    filter output (used for the search-mode dark-count estimate).
    """
    Xp = np.fft.rfft(Ypsd - Ypsd.mean(axis=0), axis=1)
    S = np.mean(np.abs(Xp) ** 2, axis=0)
    S = np.maximum(S, 1e-30 * S.max())
    s_hat = np.fft.rfft(tpl)
    wgt = np.conj(s_hat) / S
    norm = float(np.real(np.sum(wgt * s_hat)))

    def score(Y):
        return np.real(np.fft.rfft(Y, axis=1) @ wgt) / norm

    o1, o0 = score(Y1), score(Y0)
    dkc = float((o1.mean() - o0.mean()) /
                np.sqrt(0.5 * (o0.var() + o1.var())))
    run = np.fft.irfft(np.fft.rfft(Y0, axis=1) * wgt[None, :],
                       n=Y0.shape[1], axis=1) / norm
    run = run - run.mean(axis=1, keepdims=True)
    ac = np.mean(np.fft.irfft(np.abs(np.fft.rfft(run, axis=1)) ** 2,
                              n=Y0.shape[1], axis=1), axis=0)
    ac = ac / ac[0]
    lag_c = int(np.argmax(ac < np.exp(-1.0)))
    return dict(dkc=dkc, mu1=float(o1.mean()), mu0=float(o0.mean()),
                s1=float(o1.std()), s0=float(o0.std()),
                eff=float(np.mean(o1 > thr)),
                dark=float(np.mean(o0 > thr)), corr_lag=lag_c,
                out1=o1, out0=o0, wgt=wgt, norm=norm)
