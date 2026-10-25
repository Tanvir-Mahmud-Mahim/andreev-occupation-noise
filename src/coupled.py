"""Coupled two-node thermal model of the junction (revision, 2026).

The graphene electron gas (heat capacity C_e, cooled to the phonon bath
through G_ep) and the Andreev levels (heat capacity C_A) exchange
quasiparticles, and therefore energy, on the time scale tau_A. In linear
response the Andreev levels behave as a second thermal node attached to
the electrons by the conductance

    G_A = C_A / tau_A .

The readout measures the occupations, i.e. the "temperature" T_A of this
second node. Two independent noise sources act: the electron-phonon
power noise (single-sided PSD 4 kB T^2 G_ep, between electrons and bath)
and the quasiparticle-exchange noise (PSD 4 kB T^2 G_A, between electrons
and Andreev levels). With T_e held fixed, the exchange noise alone gives
S_TA = 4 kB T^2 tau_A / C_A / (1 + w^2 tau_A^2), which is exactly the
Lorentzian occupation noise of Eq. (2) of the main text expressed in
temperature units. The earlier one-way model is the limit C_A -> 0 in the
electron heat balance.

Heat balance (deposit P into the electrons, readout of T_A):
    i w C_e T_e = P - G T_e - G_A (T_e - T_A) + n_ep - n_A
    i w C_A T_A = G_A (T_e - T_A) + n_A
so that T_A = [P + n_ep + n_A (G + i w C_e)/G_A] / D(w), with
    D(w) = [(i w C_e + G + G_A)(G_A + i w C_A) - G_A^2] / G_A .
The noise-equivalent power referred to the deposit is therefore
    NEP^2(w) = 4 kB T^2 [ G + (G^2 + w^2 C_e^2) / G_A ]
              + S_TA,ro |D(w)|^2 ,
where S_TA,ro is a white readout floor in units of T_A. Without the
readout term the optimal-filter energy resolution has the closed form
    sigma_E = 2 sqrt(kB T^2 C_e) [x (1 + x)]^(1/4),  x = G tau_A / C_A,
which reduces to the one-way result 2 sqrt(kB T^2 C_e) sqrt(x) for
x >> 1.
"""

import numpy as np

from constants import KB


def x_ratio(G, tauA, CA):
    """x = G tau_A / C_A = tau_A / tau_A*, with tau_A* = C_A / G."""
    return G * tauA / CA


def D_of_w(w, Ce, CA, G, GA):
    """Transfer denominator: T_A = (P + ...)/D(w)."""
    return ((1j * w * Ce + G + GA) * (GA + 1j * w * CA) - GA ** 2) / GA


def nep2(w, T, Ce, CA, G, tauA, S_TA_ro=0.0, H_ro=None):
    """Input-referred single-sided NEP^2 (W^2/Hz) of the two-node model.

    S_TA_ro: white readout floor in K^2/Hz (referred to T_A).
    H_ro: optional complex readout transfer function (e.g. resonator
    low-pass) applied to the signal before the white floor is added.
    """
    GA = CA / tauA
    out = 4.0 * KB * T ** 2 * (G + (G ** 2 + (w * Ce) ** 2) / GA)
    if S_TA_ro:
        D2 = np.abs(D_of_w(w, Ce, CA, G, GA)) ** 2
        if H_ro is not None:
            D2 = D2 / np.abs(H_ro) ** 2
        out = out + S_TA_ro * D2
    return out


def sigmaE_closed(T, Ce, CA, G, tauA):
    """Closed-form energy resolution (J), Andreev + phonon noise."""
    x = x_ratio(G, tauA, CA)
    return 2.0 * np.sqrt(KB * T ** 2 * Ce) * (x * (1.0 + x)) ** 0.25


def sigmaE_oneway_closed(T, Ce, CA, G, tauA):
    """Earlier one-way closed form (no back-action), Andreev noise only."""
    x = x_ratio(G, tauA, CA)
    return 2.0 * np.sqrt(KB * T ** 2 * Ce) * np.sqrt(x)


def sigmaE_numeric(T, Ce, CA, G, tauA, S_TA_ro=0.0, kappa=None,
                   n_f=20000):
    """Optimal-filter energy resolution (J) from the full NEP.

    sigma_E^-2 = 4 int_0^inf df / NEP^2(f). kappa (rad/s): optional
    resonator linewidth; the readout then follows the signal through
    H = 1/(1 + 2 i w / kappa) (amplitude response of a one-port
    resonator to a slow frequency modulation).
    """
    tau_min = min(Ce / (G + CA / tauA), tauA, Ce / G)
    f = np.logspace(np.log10(1e-3 / (2 * np.pi * max(tauA, CA / G,
                                                    Ce / G))),
                    np.log10(1e3 / (2 * np.pi * tau_min)), n_f)
    w = 2 * np.pi * f
    H = None if kappa is None else 1.0 / (1.0 + 2j * w / kappa)
    integ = 4.0 * np.trapezoid(1.0 / nep2(w, T, Ce, CA, G, tauA,
                                          S_TA_ro, H), f)
    return 1.0 / np.sqrt(integ)


def STA_spectrum(w, T, Ce, CA, G, tauA):
    """Equilibrium fluctuation spectrum of T_A (K^2/Hz, single-sided)
    in the coupled model (electron-phonon + exchange noise)."""
    GA = CA / tauA
    return nep2(w, T, Ce, CA, G, tauA) / np.abs(D_of_w(w, Ce, CA, G,
                                                       GA)) ** 2


def corner_time(T, Ce, CA, G, tauA, n_f=4000):
    """Half-power corner of the occupation-noise spectrum in the coupled
    model, returned as an equivalent time 1/(2 pi f_corner) (s). The
    spectrum is first reduced by its high-frequency phonon tail so that
    only the Lorentzian-like part is considered: we locate the frequency
    where S_TA falls to half of its zero-frequency value."""
    f = np.logspace(np.log10(1e-4 / (2 * np.pi * (tauA + CA / G))),
                    np.log10(1e3 / (2 * np.pi * min(tauA, Ce / G))), n_f)
    S = STA_spectrum(2 * np.pi * f, T, Ce, CA, G, tauA)
    S0 = S[0]
    i = int(np.argmax(S < 0.5 * S0))
    # log interpolation
    f1, f2, s1, s2 = f[i - 1], f[i], S[i - 1], S[i]
    fc = np.exp(np.log(f1) + (np.log(0.5 * S0) - np.log(s1)) *
                (np.log(f2) - np.log(f1)) / (np.log(s2) - np.log(s1)))
    return 1.0 / (2 * np.pi * fc)
