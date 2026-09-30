# Changelog

All notable changes to this code are listed here, newest first. The
repository has no releases or version tags; the entries below follow the
git history.

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
