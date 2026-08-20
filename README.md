# Andreev occupation noise: the sensitivity limit of proximity Josephson thermal detectors

Simulation code for the manuscript

> T. M. Mahim, A. S. M. Mohsin, and M. M. Rahman,
> "Andreev occupation noise sets the sensitivity limit of proximity
> Josephson thermal detectors".

The pipeline is a parameter-free ballistic Andreev-bound-state model of
the six graphene Josephson junction contact recipes measured by
Jung *et al.*, [Phys. Rev. Applied **26**, 014078 (2026)](https://doi.org/10.1103/9lsg-mdb8),
from which it derives the intrinsic occupation-noise limit of Andreev
thermometry, the matched-level design condition E = 2.40 kBT, predicted
resonator frequency-noise spectra, and single-microwave-photon
calorimetry budgets.

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
| `tests/` | analytic testbench (run `python tests/test_abs.py`, `python tests/test_noise.py`) |
| `scripts/make_bcs_table.py` | tabulates the universal BCS gap function (run first) |
| `scripts/exp_*.py` | experiments that generate all data in `data/` |
| `scripts/fig*.py` | figure generation |
| `scripts/make_numbers.py` | regenerates every number quoted in the manuscript |

## Reproducing everything

```bash
pip install -r requirements.txt
bash run_all.sh
```

`run_all.sh` runs, in order: the BCS gap tabulation, both testbench
suites, the three experiment scripts, the matched-design operating
points, the four figures, and `make_numbers.py`. Total runtime is a few
minutes on a laptop.

## Testbench pass criteria

Short-junction and Kulik level anchors at machine precision; ballistic
e Ic Rn = pi Delta anchor to 5e-3; free-energy continuity across
bound-state exit into the continuum; analytic occupation responsivity
to 1e-5; exact Cauchy-Schwarz bound saturation for uniform-transparency
ensembles to 1e-9; Monte Carlo Lorentzian plateau to 10%, knee
frequency to 12%, and the variance convention var = S/2t to 25%
(statistics limited).

## License

Apache-2.0. Data and archived outputs: CC-BY-4.0 (see Zenodo record).
