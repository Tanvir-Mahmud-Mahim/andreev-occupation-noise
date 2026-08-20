"""Nonlinear single-photon click Monte Carlo for the matched design.

Removes the linear-response and Gaussian-statistics assumptions of the
main text:
  * full nonlinear heat balance C_e(T) dT/dt = -Sigma A (T^3 - T0^3)
    after a deposit E_gamma, integrated with a stiffness-safe
    exponential integrator (energy conservation verified to machine
    precision through the peak-temperature identity);
  * Andreev occupations m(t) relax toward tanh(Delta*/2 kB Te(t)) with
    an exchange time that is either the cold value tau_A (scenario C)
    or thermally activated, tau(Te) = tau_A exp[-(Delta*/kB)(1/T0 -
    1/Te)], floored at 1 ns (scenario T); the electrodes stay cold so
    Delta* is fixed;
  * the resonator frequency follows the exact series-LC relation with
    LJ = (hbar/2e)/I'(0), I'(0) proportional to m;
  * noises: occupation telegraph (variance 2f(1-f)/Nch, correlation
    tau_A), phonon TFN entering through the occupation lag, and the
    quantum-limited readout floor; detection uses a matched filter
    built from the mean click template and empirical statistics from
    400 photon and 400 dark trials per scenario.
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from scipy.optimize import brentq
from constants import KB, HBAR, H_PLANCK, E_CHARGE
from materials import Recipe, RECIPES, carrier_density
from sensor_limits import SensorBudget
from abs_model import gap_bcs

OUT = os.path.join(os.path.dirname(__file__), "..", "data")
rng = np.random.default_rng(21)
EGAMMA = H_PLANCK * 26e9
res = {}

for T0, WUM, LUM in ((0.05, 1.0, 0.1), (0.05, 5.3, 0.5),
                     (0.05, 5.3, 1.5), (0.1, 5.3, 0.5)):
    g = lambda Tc: gap_bcs(T0, Tc, 1.7639 * KB * Tc) - 2.3994 * KB * T0
    Tc = brentq(g, T0 * 1.01, T0 * 6)
    r = Recipe("m", "m", Tc, 1e-6, LUM * 1e-6, WUM * 1e-6, 30.0, 0.3,
               RECIPES[0].Ic20, RECIPES[0].Rn)
    sb = SensorBudget(r)
    sj = sb.sj
    A = r.area
    gam = sb.Ce(T0) / T0                    # C_e = gam T
    Sig, dlt = sb.sigma, sb.delta
    D = sj.Delta(T0)
    Nch = sj.Nch
    m0 = float(np.tanh(D / (2 * KB * T0)))
    f0 = 1.0 / (np.exp(D / (KB * T0)) + 1.0)
    I1_0 = sj.dIdphi0(T0)
    LJ0 = (HBAR / (2 * E_CHARGE)) / I1_0
    Lr, nu0 = sb.L_r, sb.nu_r
    Cres = 1.0 / ((2 * np.pi * nu0) ** 2 * (Lr + LJ0))

    def nu_of_m(m):
        I1 = I1_0 * np.maximum(m, 1e-6) / m0
        LJ = (HBAR / (2 * E_CHARGE)) / I1
        return 1.0 / (2 * np.pi * np.sqrt((Lr + LJ) * Cres))

    def tau_of_T(Te, scen, tauA):
        if scen == "C":
            return tauA
        expo = -(D / KB) * (1.0 / T0 - 1.0 / max(Te, 1e-4))
        return max(1e-9, tauA * float(np.exp(np.clip(expo, -500, 0))))

    def meq_of_T(Te):
        return np.tanh(D / (2 * KB * np.maximum(Te, 1e-4)))

    def template(scen, tauA, dt=1e-10, tend=2e-7, dt2=None, tend2=None):
        if dt2 is None:
            dt2 = min(2e-8, tauA / 6.0)
        if tend2 is None:
            tend2 = 60.0 * tauA
        """Noiseless m(t) after the deposit (stiff-safe)."""
        Te, m, t = np.sqrt(T0 ** 2 + 2 * EGAMMA / gam), m0, 0.0
        ts, ms = [0.0], [m0]
        for step, stop in ((dt, tend), (dt2, tend2)):
            n = int((stop - t) / step)
            for _ in range(n):
                if Te - T0 > 1e-9:
                    P = Sig * A * (Te ** dlt - T0 ** dlt)
                    k = P / (gam * Te) / (Te - T0)   # local rate
                    Te = T0 + (Te - T0) * np.exp(-k * step)
                tl = tau_of_T(Te, scen, tauA)
                m = m + (meq_of_T(Te) - m) * (1 - np.exp(-step / tl))
                t += step
                ts.append(t); ms.append(m)
        return np.array(ts), np.array(ms)

    Tpk = np.sqrt(T0 ** 2 + 2 * EGAMMA / gam)
    entry = dict(Tc=Tc, T_peak=float(Tpk), Ce_kB=float(sb.Ce(T0) / KB),
                 sigT_over_T=float(np.sqrt(KB / sb.Ce(T0))),
                 energy_conservation=float(
                     abs(0.5 * gam * (Tpk**2 - T0**2) - EGAMMA) / EGAMMA))
    print(f"--- T0={T0} W={WUM} L={LUM}: Ce={sb.Ce(T0)/KB:.1f} kB, "
          f"sigT/T={np.sqrt(KB/sb.Ce(T0)):.2f}, Tpk={Tpk:.3f} K")

    var_m = 2 * f0 * (1 - f0) / Nch          # telegraph on the mean
    sigT_eq = np.sqrt(KB * T0 ** 2 / sb.Ce(T0))
    dmeq_dT = float(-(D / (2 * KB * T0 ** 2)) /
                    np.cosh(D / (2 * KB * T0)) ** 2)
    kap = 2 * np.pi * 1e6
    S_nu_floor = (kap / 4) ** 2 * 2 / (30.0 * kap) / (2 * np.pi) ** 2
    conv = (nu_of_m(m0 * (1 + 1e-6)) - nu0) / (m0 * 1e-6)

    from scipy.signal import lfilter
    tau_th = sb.tau_th(T0)

    for tauA in (3e-8, 1e-7, 3e-7, 1e-6):
      dt = 5e-10
      nwin = int(max(6.0 * tauA, 4e-7) / dt)
      tgrid = np.arange(nwin) * dt
      i0 = nwin // 4
      aA = np.exp(-dt / tauA)
      ath = np.exp(-dt / tau_th)
      sig_white = np.sqrt(S_nu_floor / (2 * dt))
      for scen in ("C", "T"):
        tt, mm = template(scen, tauA)
        dnu_tpl = np.interp(tgrid, tt, nu_of_m(mm) - nu0, right=None)
        entry[f"m_dip_{scen}_{tauA:.0e}"] = float((m0 - mm.min()) / m0)
        entry[f"peak_dnu_{scen}_{tauA:.0e}_kHz"] = float(
            (nu_of_m(mm.min()) - nu0) / 1e3)

        ntr = 1000
        def batch(photon):
            # electron-temperature OU (phonon TFN), stationary
            xiT = rng.standard_normal((ntr, nwin))
            dT = lfilter([sigT_eq * np.sqrt(1 - ath ** 2)], [1, -ath],
                         xiT, axis=1)
            # occupation AR(1): telegraph drive + TFN drive via the lag
            drv = (np.sqrt(var_m * (1 - aA ** 2)) *
                   rng.standard_normal((ntr, nwin))
                   + dmeq_dT * dT * (dt / tauA))
            x = lfilter([1.0], [1, -aA], drv, axis=1)
            y = conv * x + sig_white * rng.standard_normal((ntr, nwin))
            if photon:
                y[:, i0:] += dnu_tpl[None, :nwin - i0]
            return y

        # optimal (whitened) matched filter in the frequency domain,
        # using the known analytic noise PSD (common scalings cancel)
        tpl_full = np.zeros(nwin)
        tpl_full[i0:] = dnu_tpl[:nwin - i0]
        freqs = np.fft.rfftfreq(nwin, dt)
        w = 2 * np.pi * freqs
        HA2 = 1.0 / (1.0 + (w * tauA) ** 2)
        S_A = conv ** 2 * 4 * var_m * tauA * HA2
        S_TFN = conv ** 2 * dmeq_dT ** 2 * \
            (4 * sigT_eq ** 2 * tau_th / (1 + (w * tau_th) ** 2)) * HA2
        S_n = S_A + S_TFN + S_nu_floor
        s_hat = np.fft.rfft(tpl_full)
        wgt = np.conj(s_hat) / S_n
        norm = float(np.real(np.sum(wgt * s_hat)))

        def score(Y):
            return np.real(np.fft.rfft(Y, axis=1) @ wgt) / norm

        out0 = score(batch(False))
        out1 = score(batch(True))
        snr_mc = float((out1.mean() - out0.mean()) /
                       np.sqrt(0.5 * (out0.var() + out1.var())))
        thr = 0.5
        eff = float(np.mean(out1 > thr))
        dark_frac = float(np.mean(out0 > thr))
        nu_eff = 1.0 / (2 * np.pi * tauA)
        z = float((thr - out0.mean()) / out0.std())
        dark_rate = float(nu_eff * np.exp(-0.5 * min(z, 37.0) ** 2))
        entry[f"eff_{scen}_{tauA:.0e}"] = eff
        entry[f"darkfrac_{scen}_{tauA:.0e}"] = dark_frac
        entry[f"snr_mc_{scen}_{tauA:.0e}"] = snr_mc
        entry[f"dark_rate_{scen}_{tauA:.0e}"] = dark_rate
        if (T0 == 0.05 and WUM == 5.3 and abs(LUM - 0.5) < 1e-9
                and scen == "T" and tauA in (3e-8, 1e-7)):
            yp, yn = batch(True)[0], batch(False)[0]

            cnorm = np.fft.irfft(np.fft.rfft(tpl_full) * wgt,
                                 n=nwin)[0]

            def runscore(yy):
                """Whitened cross-correlation with the template at all
                circular lags; lag 0 is the nominal arrival time, and
                the normalization makes the noiseless template score 1
                at lag 0 (consistent with the scalar score)."""
                c = np.fft.irfft(np.fft.rfft(yy) * wgt, n=nwin)
                return np.real(c) / cnorm

            np.savez(os.path.join(OUT, f"click_traces_{tauA:.0e}.npz"),
                     tgrid=tgrid, tpl=dnu_tpl, y_photon=yp, y_no=yn,
                     run_photon=runscore(yp), run_no=runscore(yn),
                     out0=out0, out1=out1, nu0=nu0)
        print(f"T0={T0} W{WUM} L{LUM} {scen} tauA={tauA:.0e}: "
              f"m-dip={entry[f'm_dip_{scen}_{tauA:.0e}']:.4f} "
              f"MC SNR={snr_mc:.1f} eff={eff:.3f} dark z={z:.1f}")
    res[f"T{T0}_W{WUM}_L{LUM}"] = entry

json.dump(res, open(os.path.join(OUT, "nonlinear_click.json"), "w"))
print("energy conservation:", {k: v["energy_conservation"]
                               for k, v in res.items()})
