"""Self-consistent single-photon click Monte Carlo (revision, 2026).

Replaces the earlier exp_nonlinear_click.py results in the paper. For the
matched 5.3 x 0.5 um^2 design at 50 mK and further designs it computes
the Kumar-Carroll index d_KC (our Monte Carlo SNR) of a 26 GHz photon
with the whitened matched filter, for:
  (1) a step-by-step breakdown from the earlier model to the full model
      (nonlinear frequency map applied to the noise, nonstationary
      occupation noise, energy exchange with the Andreev levels);
  (2) the full model with a resonator readout of finite linewidth kappa
      and an amplifier with n_add added noise quanta;
  (3) d_KC versus tau_A for several geometries and both scenarios;
  (4) d_KC for several induced gaps at fixed tau_A;
  (5) d_KC for electron-phonon exponents delta = 3 to 6 at fixed
      conductance at 0.19 K (option --eph-only reruns only this part).
Model: src/click_sim.py (linear-noise approximation around the nonlinear
mean trajectory). Records are long compared with the slow relaxation
time (C_e + C_A)/G_ep so that the frequency-domain filter resolves the
slow noise.
"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm as gauss
from constants import KB, H_PLANCK
from materials import Recipe, RECIPES
from sensor_limits import SensorBudget
from abs_model import gap_bcs
from click_sim import ClickModel, resonator_phase, matched_filter_stats

OUT = os.path.join(os.path.dirname(__file__), "..", "data")
EGAMMA = H_PLANCK * 26e9
NBAR, KAPPA0 = 30.0, 2 * np.pi * 1e6        # floor of the ideal readout
QUICK = "--quick" in sys.argv
NTR = 300 if QUICK else 1000
NPSD = 300 if QUICK else 1000


def design(T0, W, L, y=2.3994, nu=None, delta=3, T_anchor=0.19):
    g = lambda Tc: gap_bcs(T0, Tc, 1.7639 * KB * Tc) - y * KB * T0
    Tc = brentq(g, T0 * 1.01, T0 * 40)
    r = Recipe("m", "m", Tc, 1e-6, L * 1e-6, W * 1e-6, 30.0, 0.3,
               RECIPES[0].Ic20, RECIPES[0].Rn)
    D = gap_bcs(T0, Tc, 1.7639 * KB * Tc)
    if nu is None:     # readout tone well below the pair threshold 2E/h
        nu = 3e9 if 2 * D / H_PLANCK < 12e9 else 6e9
    sb = SensorBudget(r, nu_r=nu)
    if delta != 3:
        # Sigma_delta chosen so that G_ep(T_anchor) equals the delta = 3
        # value (Sigma = 2 W m^-2 K^-3, from the fit of Lee et al.)
        G3 = 3 * sb.sigma * r.area * T_anchor ** 2
        sb = SensorBudget(r, nu_r=nu, delta=delta,
                          sigma=G3 / (delta * r.area *
                                      T_anchor ** (delta - 1)))
    sj = sb.sj
    cm = ClickModel(T0, sj.Delta(T0), sj.Nch, sb.Ce(T0) / T0,
                    sb.sigma * r.area, sb.delta, nu, sb.L_r,
                    sj.dIdphi0(T0))
    s = sj.andreev_sums(1e-4, T0, "L")
    info = dict(T0=T0, W_um=W, L_um=L, Tc=Tc, DkT=D / (KB * T0),
                delta=delta,
                nu_r_GHz=nu / 1e9, Ce_kB=sb.Ce(T0) / KB,
                CA_kB=s["C_A"] / KB, G=sb.Gep(T0),
                tau_slow_ns=(sb.Ce(T0) + s["C_A"]) / sb.Gep(T0) * 1e9,
                pair_GHz=2 * D / H_PLANCK / 1e9)
    return cm, info


def run_case(cm, info, scen, tauA, rng, coupled=True, nonstationary=True,
             readout="ideal", kappa=None, nadd=0.0, keep=False,
             dt=0.5e-9, Tmax=1.2e-5):
    """One Monte Carlo case. The record spans ten slow relaxation times
    (C_e + C_A)/G_ep, but at most Tmax (memory); a record shorter than
    that only removes slow-noise information and so lowers d_KC."""
    T = min(max(10 * info["tau_slow_ns"] * 1e-9, 8 * tauA, 2e-6), Tmax)
    nwin = int(T / dt)
    tg = np.arange(nwin) * dt
    i0 = nwin // 8
    M1, mt, _ = cm.records(NTR, tg, i0, EGAMMA, scen, tauA, rng,
                           coupled, nonstationary)
    M0, _, _ = cm.records(NTR, tg, i0, 0.0, scen, tauA, rng, coupled,
                          nonstationary)
    Mp, _, _ = cm.records(NPSD, tg, i0, 0.0, scen, tauA, rng, coupled,
                          nonstationary)
    f = lambda M: cm.nu_of_m(M) - cm.nu0
    if readout == "ideal":
        Snu = (KAPPA0 / 4) ** 2 * 2 / (NBAR * KAPPA0) / (2 * np.pi) ** 2
        wn = lambda n: np.sqrt(Snu / (2 * dt)) * rng.standard_normal(
            (n, nwin))
        Y1, Y0, Yp = f(M1) + wn(NTR), f(M0) + wn(NTR), f(Mp) + wn(NPSD)
        tpl = f(mt[None, :])[0]
    else:
        Sth = 2 * (1 + 2 * nadd) / (NBAR * kappa)
        wn = lambda n: np.sqrt(Sth / (2 * dt)) * rng.standard_normal(
            (n, nwin))
        ph = lambda M: resonator_phase(f(M), dt, kappa)
        Y1, Y0, Yp = ph(M1) + wn(NTR), ph(M0) + wn(NTR), \
            ph(Mp) + wn(NPSD)
        tpl = ph(mt[None, :])[0]
    st = matched_filter_stats(Y1, Y0, Yp, tpl)
    # search mode: false-alarm rate of the running filter output at the
    # half-photon threshold, from the Gaussian tail of the dark scores
    # and the 1/e correlation time of the running output
    z = (0.5 - st["mu0"]) / st["s0"]
    tc = max(st["corr_lag"], 1) * dt
    rate = float(gauss.sf(z) / tc)
    out = dict(dkc=st["dkc"], eff=st["eff"], dark=st["dark"],
               mu1=st["mu1"], s0=st["s0"], s1=st["s1"],
               corr_time_ns=tc * 1e9, search_rate_hz=rate,
               peak_dnu_MHz=float(-f(mt[None, :]).min() / 1e6),
               record_us=T * 1e6)
    if keep:
        out["_arrays"] = dict(tgrid=tg, tpl=tpl, y1=Y1[0], y0=Y0[0],
                              out0=st["out0"], out1=st["out1"],
                              wgt=st["wgt"], norm=st["norm"], i0=i0)
    return out


def run_eph(rng, t0):
    """Electron-phonon exponent: Sigma_delta anchored at 0.19 K. Slower
    cooling at 50 mK lengthens the slow relaxation time, so the records
    are capped at 12 us (conservative for delta >= 5)."""
    out = []
    for delta in (3, 4, 5, 6):
        cm, info = design(0.05, 5.3, 0.5, delta=delta)
        for scen in ("C", "T"):
            r = run_case(cm, info, scen, 3e-8, rng)
            out.append(dict(design=info, scen=scen, tauA_ns=30.0, **r))
            print(f"[{time.time()-t0:6.0f}s] delta={delta} {scen}: "
                  f"dKC={r['dkc']:.2f} eff={r['eff']:.3f}", flush=True)
    return out


def main():
    rng = np.random.default_rng(2026)
    res = {}
    t0 = time.time()

    # (1) breakdown for the matched 5.3 x 0.5 um^2 design at 50 mK --------
    cm6, info6 = design(0.05, 5.3, 0.5, nu=6e9)
    cm3, info3 = design(0.05, 5.3, 0.5, nu=3e9)
    res["design_main"] = info3
    steps = [("earlier model, 6 GHz", cm6, info6, False, False),
             ("earlier model, 3 GHz", cm3, info3, False, False),
             ("+ nonstationary noise", cm3, info3, False, True),
             ("+ energy exchange (full)", cm3, info3, True, True)]
    bd = []
    ex = {}
    for tauA in (3e-8, 1e-7):
        for scen in ("C", "T"):
            for lab, cm, info, cpl, nst in steps:
                keep = (tauA == 3e-8 and lab.startswith("+ energy"))
                r = run_case(cm, info, scen, tauA, rng, cpl, nst,
                             keep=keep)
                if keep:
                    a = r.pop("_arrays")
                    np.savez(os.path.join(
                        OUT, f"click_coupled_{scen}_3e-08.npz"), **a)
                    ex[scen] = r
                bd.append(dict(step=lab, scen=scen, tauA_ns=tauA * 1e9,
                               **r))
                print(f"[{time.time()-t0:6.0f}s] {lab:28s} {scen} "
                      f"{tauA*1e9:4.0f} ns: dKC={r['dkc']:.2f} "
                      f"eff={r['eff']:.3f} dark={r['dark']:.3f} "
                      f"peak={r['peak_dnu_MHz']:.1f} MHz", flush=True)
    res["breakdown"] = bd

    # (2) finite resonator linewidth and added noise ----------------------
    ro = []
    for scen in ("T", "C"):
        for kHz in (1e6, 1e7, 5e7, 1e8, 3e8):
            for nadd in (0.0, 14.0):
                r = run_case(cm3, info3, scen, 3e-8, rng, readout="res",
                             kappa=2 * np.pi * kHz, nadd=nadd)
                ro.append(dict(scen=scen, kappa_MHz=kHz / 1e6,
                               nadd=nadd, **r))
                print(f"[{time.time()-t0:6.0f}s] resonator {scen} "
                      f"kappa/2pi={kHz/1e6:5.0f} MHz nadd={nadd:4.1f}: "
                      f"dKC={r['dkc']:.2f} eff={r['eff']:.3f} "
                      f"dark={r['dark']:.3f}", flush=True)
    res["readout"] = ro

    # (3) dKC versus tauA, full model, ideal readout -----------------------
    geoms = [(0.05, 5.3, 0.5), (0.05, 2.0, 0.5), (0.05, 5.3, 1.5),
             (0.1, 5.3, 0.5)]
    vs = []
    for (T0, W, L) in geoms:
        cm, info = design(T0, W, L)
        for scen in ("C", "T"):
            for tauA in (1e-8, 3e-8, 1e-7, 3e-7, 1e-6):
                r = run_case(cm, info, scen, tauA, rng)
                vs.append(dict(design=info, scen=scen,
                               tauA_ns=tauA * 1e9, **r))
                print(f"[{time.time()-t0:6.0f}s] {T0} {W}x{L} {scen} "
                      f"{tauA*1e9:5.0f} ns: dKC={r['dkc']:.2f}",
                      flush=True)
    res["vs_tauA"] = vs

    # (4) induced gap at fixed tauA ---------------------------------------
    gp = []
    for (W, L) in ((5.3, 0.5), (2.0, 0.5)):
        for y in (2.4, 4.0, 6.0):
            cm, info = design(0.05, W, L, y=y)
            for scen in ("C", "T"):
                r = run_case(cm, info, scen, 3e-8, rng)
                gp.append(dict(design=info, scen=scen, tauA_ns=30.0,
                               **r))
                print(f"[{time.time()-t0:6.0f}s] gap y={y} {W}x{L} "
                      f"{scen}: dKC={r['dkc']:.2f} eff={r['eff']:.3f}",
                      flush=True)
    res["gap_scan"] = gp

    # (5) electron-phonon exponent (Sigma anchored at 0.19 K) ------------
    res["eph"] = run_eph(rng, t0)

    res["example"] = ex      # = the full-model 30 ns breakdown rows

    json.dump(res, open(os.path.join(OUT, "click_coupled.json"), "w"),
              indent=1)
    print("done in", time.time() - t0, "s")


if __name__ == "__main__":
    if "--eph-only" in sys.argv:
        # append the exponent study to an existing results file
        fn = os.path.join(OUT, "click_coupled.json")
        res = json.load(open(fn))
        res["eph"] = run_eph(np.random.default_rng(7), time.time())
        json.dump(res, open(fn, "w"), indent=1)
    else:
        main()
