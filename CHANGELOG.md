# Changelog

All notable changes to this code are listed here, newest first. The
repository has no releases or version tags; the entries below follow the
git history.

## Documentation of the cooling law (25 October 2026, Zenodo v1.13)

No model, data or result changed; the four test files give output
identical to their reference files.

- `src/materials.py`, `scripts/exp_revision.py`,
  `scripts/exp_click_coupled.py`, `DETAILS.md`: comments now say that
  Lee *et al.* fitted their bolometer with the exponent set to 3 and
  Sigma between 2.04 and 3.30 W m^-2 K^-3 (Sigma = 2 is close to the low
  end), rather than calling these values measured; the 0.6 ns thermal
  time of Lee *et al.* is an estimate.
- Documentation gives the v1.13 DOI (doi:10.5281/zenodo.23648106).

## Figure layout (25 October 2026, Zenodo v1.12)

No model, data or result changed.

- `scripts/fig_device.py`: the schematic is drawn at its printed size
  (0.72 of the text width) with all text at 6.3 pt or more; the device is
  no longer cut at the left edge and the key box no longer at the right.
- `scripts/fig1.py`: more room for the x-axis label.
- `scripts/fig2.py`: panel labels moved clear of the tick labels; the
  "phonon limited" label no longer crosses the boundary curve.
- `scripts/fig3.py`: legend of (a) moved clear of the dotted optimum line
  and the annotation below the curves; the 26 GHz photon line is named in
  a legend instead of a label over the bars of (c).
- `scripts/fig4.py`: the statistics in (b) no longer overlap the
  histogram or the threshold line.
- `scripts/figS2.py`: legends moved clear of the data in all panels.
- Documentation gives the v1.12 DOI (doi:10.5281/zenodo.23647063).

## Revision after peer review (25 October 2026)

Changes made while answering the reviewers of Materials Today Physics.

**New code**
- `src/coupled.py`: electrons and Andreev levels as two coupled thermal
  bodies linked by C_A/tau_A; noise-equivalent power, closed-form energy
  resolution sigma_E = 2 sqrt(kB T^2 C_e) [x(1+x)]^(1/4) with
  x = G_ep tau_A / C_A, numerical filter with readout floor and resonator
  linewidth, corner time of the noise spectrum.
- `src/click_sim.py`: self-consistent click simulation. The mean response
  is fully nonlinear and returns the energy carried by exchanged
  quasiparticles to the electrons; fluctuations use the instantaneous
  rates and variances (linear-noise approximation, exact matrix
  exponential per step); optional resonator readout with added noise;
  whitened matched filter and Kumar-Carroll index.
- `scripts/exp_click_coupled.py` replaces `scripts/exp_nonlinear_click.py`
  (removed, together with `data/nonlinear_click.json` and both
  `click_traces_*.npz`). New outputs `click_coupled.json` and
  `click_coupled_{T,C}_3e-08.npz`.
- `scripts/exp_revision.py`: electron-phonon exponents 3 to 6;
  level-dependent (activated) exchange and the gap optimum when exchange is
  mediated by thermal quasiparticles; two-body energy resolution of
  recipes and matched designs; comparison with reported 1/f noise; Andreev
  transition frequencies; energy balance of the photon response. Output
  `revision.json`.
- `tests/test_coupled.py` with reference output.

**Changed**
- `src/sensor_limits.py`: the numerical energy-resolution filter now
  integrates to 1000/(2 pi min(tau_A, tau_th)) with 20000 points. The
  earlier limit of 10/(2 pi min(...)) truncated the integral and
  overstated sigma_E by about 4% at tau_A = 1 us and about 10% at 30 ns.
- `scripts/exp_matched_points.py`, `exp_calorimetry.py`: matched designs
  are now 5.3 x 0.2 um at 100 mK and 5.3 x 0.5 um at 50 mK (C_e about 8 to
  10 kB). The 1 x 0.1 um design used before has C_e < 1 kB and is outside
  linear response. The 50 mK designs are read at 3 GHz, below the pair
  transition 2 Delta*/h = 5 GHz (6 GHz before).
- `scripts/exp_design.py`, `make_numbers.py`, figure scripts: closed-form
  two-body energy resolution; new and renamed numbers (90 values, 77 LaTeX
  macros); figures redrawn at printed size with larger fonts; new Fig. 4
  and Fig. S2.
- `run_all.sh`: adds `test_coupled.py`, `exp_revision.py` and
  `exp_click_coupled.py`.
- `README.md` and `DETAILS.md` updated for the new code; the Zenodo
  archive of this version is v1.11 (doi:10.5281/zenodo.23646142).

## Revised manuscript (2 October 2026)

- Documentation follows the revised manuscript title, "Andreev-level
  occupation noise in graphene Josephson thermal detectors: a sensitivity
  floor and how to measure it". No model, data or result changed.
- `scripts/make_numbers.py`: `LxiMoRe` is now computed from the tabulated
  L and coherence length; new derived numbers for the revised text
  (`tauVisNs`, `tauAstarTiAlAuNs`, `occSharePct`, `nlDark100`,
  `nlShiftMHz`); `nlCeSmall` is given to two decimals (0.38). 56 macros.
- `tests/test_abs.py`: two new checks, free-energy continuity when a level
  leaves the gap (residual 0.4% of the jump) and the BCS gap against its
  low-temperature and Ginzburg-Landau limits; reference output updated.
- `src/abs_model.py`, `src/materials.py`, `scripts/exp_nonlinear_click.py`:
  docstrings corrected (where checks are made; Sigma range 2.0 to 3.3;
  thermal time 0.8 ns versus the 0.6 ns of Lee et al.; the click Monte
  Carlo is an idealized readout and runs 1000 trials per case).
- `data/calibration.json`: `s_full` for Ti/Al(thick) (0.659 -> 0.658) and
  MoRe (0.501 -> 0.493) set to the recomputed values.
- README updated accordingly; `CITATION.cff` gives the new title.

## Fixes (30 September 2026)

- `scripts/make_numbers.py` now writes the LaTeX macros to
  `data/numbers.tex`, next to `data/numbers.json`. Before, it wrote
  `paper/numbers.tex`, in a `paper/` folder that is not in the repository
  and that neither the script nor `run_all.sh` created, so the last step of
  `run_all.sh` stopped with `FileNotFoundError` unless `mkdir -p paper` had
  been run first. The content of the file is unchanged (checked: identical
  byte for byte to the file written by the previous version from the same
  data).
- `.gitignore`: removed the six `paper/...` patterns; added
  `data/numbers.tex`, so the generated macro file is not committed.
- `requirements.txt`: `numpy>=1.24` changed to `numpy>=2.0`, because the
  code calls `numpy.trapezoid`, which exists only from NumPy 2.0. With that,
  `scipy>=1.10` became `scipy>=1.13` and `matplotlib>=3.7` became
  `matplotlib>=3.8.4`, the first releases that work with NumPy 2.
- README: removed the `mkdir -p paper` steps and the warning about the
  missing `paper/` folder, gave the new file location and minimum versions,
  and recorded a check with exactly these minimum versions.

## Documentation (30 September 2026)

Documentation only; no code, data or results changed.

- README rewritten as a step-by-step guide: installation, three ways to use
  the code (with run times measured on a shared 2-core machine), a table of
  the scripts, a table of which script makes which figure, the parameter
  sources, and the built-in checks.
- README now notes two known problems, without changing the code:
  `make_numbers.py` needs a `paper/` folder that is not in the repository,
  and the code needs NumPy 2.0 or newer although `requirements.txt` allows
  `numpy>=1.24`.
- Added this CHANGELOG and `CITATION.cff`.

## Initial code (21 August 2026)

Code for the manuscript "Andreev occupation noise sets the sensitivity
limit of proximity Josephson thermal detectors". All of the following was
added on 21 August 2026.

**Model and checks**
- `src/`: constants, the six measured junction recipes (Jung *et al.*),
  the exact finite-length Andreev-level solver, the short-junction model,
  finite-length deficits, sensor noise budgets and the telegraph Monte
  Carlo.
- `tests/test_abs.py` and `tests/test_noise.py` with reference outputs.
- `scripts/`: BCS gap table, the four calculation scripts
  (`exp_universal.py`, `exp_design.py`, `exp_matched_points.py`,
  `exp_calorimetry.py`), figure scripts `fig1.py` to `fig4.py` and
  `figS1.py`, and `make_numbers.py`.
- `data/`: `bcs_gap_table.npz`, `calibration.json`, `numbers.json`.
- `README.md`, `LICENSE` (Apache-2.0), `requirements.txt`, `run_all.sh`.

**Generalized noise model and limits study**
- `src/noise_general.py`: pair processes, activated exchange time and
  non-thermal filling; `tests/test_general.py` with reference output.
- `scripts/exp_limits.py` (checks of the approximations) and
  `scripts/exp_nonlinear_click.py` (nonlinear single-photon click Monte
  Carlo); revised figures 1 to 4; new `figS2.py`.
- Stored outputs `limits.json`, `nonlinear_click.json` and the two click
  trace files; `numbers.json` refreshed; README and `run_all.sh` updated.

**Device schematic**
- `scripts/fig_device.py`: 3D device schematic, added to `run_all.sh`,
  refined in several steps, and combined with `fig1.py` as panels (a) to (c)
  of Fig. 1.

**Further limits checks**
- Diffusive (Dorokhov) transparency distribution and the time needed to
  measure tau_A from the spectral knee (`exp_limits.py`).
- Random-ensemble test of the Cauchy-Schwarz bound (`test_general.py`).
- `abs_model.py`: split of the current into below-gap and above-gap parts
  (`continuum_share`, `free_energy_components`).
- Multi-level knee worst case and degraded readout floor checks
  (`exp_limits.py`); `limits.json` and `numbers.json` refreshed.
