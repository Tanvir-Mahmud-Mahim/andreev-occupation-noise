# Andreev-level occupation noise in graphene Josephson thermal detectors

Simulation code for the manuscript

> T. M. Mahim, A. S. M. Mohsin, and M. M. Rahman,
> "Andreev-level occupation noise in graphene Josephson thermal
> detectors: a sensitivity floor and how to measure it".

(Earlier versions of this archive accompanied a draft titled "Andreev
occupation noise sets the sensitivity limit of proximity Josephson
thermal detectors".)

The pipeline is a ballistic short-junction Andreev-bound-state model
(one current-scale calibration per recipe, which cancels from the main
results) of the six graphene Josephson junction contact recipes measured by
Jung *et al.*, [Phys. Rev. Applied **26**, 014078 (2026)](https://doi.org/10.1103/9lsg-mdb8),
from which it derives the occupation-noise floor of Andreev
thermometry, the regime map against the electron-phonon thermal floor,
the predicted resonator frequency-noise spectra and a simulated
measurement of the exchange time tau_A from their corner frequency, the
heat-capacity-peak design condition E = 2.40 kBT, and
single-microwave-photon calorimetry budgets.

## Layout

| path | content |
|---|---|
| `src/constants.py` | physical constants |
| `src/materials.py` | measured recipes (Jung *et al.*), graphene thermodynamics |
| `src/abs_model.py` | exact finite-length ABS solver, continuum free energy, BCS gap |
| `src/short_junction.py` | closed-form short-junction ensemble, occupation sums, bound |
| `src/finiteL.py` | bound-saturation deficits of the finite-length model |
| `src/sensor_limits.py` | noise budgets, frequency-noise spectra, matched filter |
| `src/montecarlo.py` | telegraph-noise Monte Carlo |
| `src/noise_general.py` | pair-process master equation, activated exchange, nonequilibrium floor |
| `tests/` | analytic testbench (run `python tests/test_abs.py`, `test_noise.py`, `test_general.py`) |
| `scripts/make_bcs_table.py` | tabulates the universal BCS gap function (run first) |
| `scripts/exp_*.py` | experiments that generate all data in `data/` |
| `scripts/exp_limits.py` | quantitative resolution of the model approximations (SM Sec. VII) |
| `scripts/exp_nonlinear_click.py` | nonlinear single-photon click Monte Carlo (SM Sec. VIII) |
| `scripts/fig*.py` | figure generation |
| `scripts/make_numbers.py` | regenerates every number quoted in the manuscript |

## Reproducing everything

```bash
pip install -r requirements.txt
bash run_all.sh
```

`run_all.sh` runs, in order: the BCS gap tabulation, the three
testbench suites, the experiment scripts (including the
approximation-resolution study and the nonlinear click Monte Carlo),
the figures, and `make_numbers.py`. Total runtime is a few minutes
to a few tens of minutes depending on the machine.

## Testbench pass criteria

Short-junction and Kulik level anchors at machine precision; ballistic
e Ic Rn = pi Delta anchor to 5e-3; free-energy continuity across
bound-state exit into the continuum; analytic occupation responsivity
to 1e-5; the Cauchy-Schwarz bound inequality on randomized level
ensembles, with exact equality for couplings proportional to the level
energies; exact bound saturation for uniform-transparency ensembles to
1e-9; Monte Carlo Lorentzian plateau to 10%, knee
frequency to 12%, and the variance convention var = S/2t to 25%
(statistics limited); pair-process generalization: exact singles limit,
equilibrium-variance invariance, monotone shortening of the effective
correlation time, generator probability conservation, and the
zero-frequency spectral density against a four-state Monte Carlo;
nonlinear click Monte Carlo: energy conservation of the exponential
integrator to machine precision.

## License

Apache-2.0. Data and archived outputs: CC-BY-4.0 (see Zenodo record).
