# Andreev Occupation Noise: Simulation Code

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://www.apache.org/licenses/LICENSE-2.0)

Code for the manuscript **"Andreev occupation noise sets the sensitivity
limit of proximity Josephson thermal detectors"**
by T. M. Mahim, A. S. M. Mohsin, and M. M. Rahman.

- Repository: https://github.com/Tanvir-Mahmud-Mahim/andreev-occupation-noise
- Measured devices that the model is built on: W. Jung *et al.*,
  "Engineering Andreev bound states for thermal sensing in proximity
  Josephson junctions", Phys. Rev. Applied **26**, 014078 (2026),
  https://doi.org/10.1103/9lsg-mdb8 (preprint: arXiv:2503.06850)

The scripts produce the data files in `data/`, the figures of the
manuscript, and (according to `make_numbers.py`) every number quoted in
the manuscript (`data/numbers.json`, see
[Section 5](#5-the-scripts-step-by-step)).

---

## Contents

1. [The idea in one minute](#1-the-idea-in-one-minute)
2. [What is in this repository](#2-what-is-in-this-repository)
3. [Installation](#3-installation)
4. [Quick start: three ways to use the code](#4-quick-start-three-ways-to-use-the-code)
5. [The scripts, step by step](#5-the-scripts-step-by-step)
6. [Which script makes which figure](#6-which-script-makes-which-figure)
7. [The Python modules](#7-the-python-modules)
8. [Where the numbers come from](#8-where-the-numbers-come-from)
9. [Built-in checks](#9-built-in-checks)
10. [Notes on the calculations](#10-notes-on-the-calculations)
11. [Version history](#11-version-history)
12. [How to cite](#12-how-to-cite)
13. [License and contact](#13-license-and-contact)

---

## 1. The idea in one minute

A **proximity Josephson junction** is a short strip of a normal conductor
(here graphene) between two superconducting contacts. The superconductors
make a supercurrent flow through the graphene. Inside the junction, this
current is carried by discrete energy levels called
**Andreev bound states**. Each level can be empty or can hold stray
excitations (quasiparticles); a filled level carries less current. Since
the chance that a level is filled depends on temperature, the current and
the junction's inductance change with temperature, and the junction works
as a thermometer and, in principle, as a photon detector.

This project points out that the same filling of the levels also **flips
randomly** in time, even at a perfectly steady temperature. This random
"occupation noise" sets a floor on how finely the junction can measure
temperature. The time over which a level keeps its filling is called the
**exchange time** `tau_A` in the code; it is not known for these devices,
so the code treats it as a free parameter and scans it.

Starting from the measured properties of six graphene junctions (Jung *et
al.*, above), the code computes:

- the energy levels and supercurrent of each junction, with an exact
  finite-length solver and a simpler short-junction formula,
- the temperature response (responsivity) and the occupation noise of the
  levels, and a lower bound on the temperature resolution set by that noise
  (a Cauchy-Schwarz inequality over the levels),
- a comparison with the other intrinsic noise, the thermal fluctuations of
  the graphene electrons exchanging heat with the lattice (called phonon
  thermal-fluctuation noise, "phonon TFN", in the code),
- the design condition for the best energy resolution, a level energy of
  about 2.40 kB T (kB T is the thermal energy at temperature T),
- the frequency noise that the occupation noise would produce in a
  microwave resonator used to read out the junction (a measurable
  prediction),
- energy-resolution budgets for detecting a single 26 GHz microwave photon,
  including a full nonlinear Monte Carlo (random-sampling) simulation of
  single-photon "clicks",
- checks of the approximations (pair processes, energy-dependent exchange
  times, non-thermal filling, spread of contact transparency, the part of
  the current carried above the gap).

---

## 2. What is in this repository

```
andreev-occupation-noise/
|-- README.md             this guide
|-- CHANGELOG.md          what changed, from the git history
|-- CITATION.cff          citation details (drives the "Cite this repository" button)
|-- LICENSE               Apache-2.0 license
|-- requirements.txt      Python packages to install
|-- run_all.sh            runs every step in order (bash)
|-- src/
|   |-- constants.py      physical constants (CODATA 2018) and graphene values
|   |-- materials.py      the six measured junction recipes; graphene heat capacity and cooling
|   |-- abs_model.py      exact finite-length Andreev-level solver, free energy, BCS gap
|   |-- short_junction.py closed-form short-junction model: responsivity, noise sums, bound
|   |-- finiteL.py        how far the finite-length levels fall short of the bound
|   |-- sensor_limits.py  noise budgets, resonator frequency-noise spectra, energy resolution
|   |-- montecarlo.py     random telegraph (two-state) simulation and spectrum estimate
|   `-- noise_general.py  pair processes, energy-dependent exchange time, non-thermal filling
|-- tests/
|   |-- test_abs.py            checks of the level solver against exact limits
|   |-- test_noise.py          checks of the noise model (bound, responsivity, Monte Carlo)
|   |-- test_general.py        checks of the generalized noise model
|   `-- reference_output_*.txt expected screen output of each test file
|-- scripts/
|   |-- make_bcs_table.py        tabulate the temperature dependence of the superconducting gap
|   |-- exp_universal.py         responsivity versus temperature for many transparencies
|   |-- exp_design.py            recipe budgets, noise spectra, regime map, design scans
|   |-- exp_matched_points.py    best-design operating points at 50 and 100 mK
|   |-- exp_calorimetry.py       energy resolution versus tau_A; simulated detection trace
|   |-- exp_limits.py            quantitative checks of the model approximations
|   |-- exp_nonlinear_click.py   nonlinear single-photon click Monte Carlo
|   |-- fig_device.py            3D device schematic (Fig. 1a)
|   |-- fig1.py ... fig4.py      main-text figures
|   |-- figS1.py, figS2.py       supplementary figures
|   |-- figstyle.py              shared figure style
|   `-- make_numbers.py          collect every number quoted in the manuscript
`-- data/                  results stored in the repository (see below)
```

**Stored in `data/`:** `bcs_gap_table.npz` (gap table), `limits.json`,
`nonlinear_click.json`, `click_traces_3e-08.npz`, `click_traces_1e-07.npz`
(the outputs of the two slowest steps), `numbers.json` (the quoted numbers)
and `calibration.json` (see [Section 10](#10-notes-on-the-calculations)).
The other result files (`universal.json`, `design.json`,
`matched_points.json`, `calorimetry.json`) are not stored; they are made in
a few seconds (Way B below). `make_numbers.py` also writes the LaTeX macros
`numbers.tex` into `data/`; it is not stored either (it is listed in
`.gitignore`). Figures go to `figures/`, which is also not stored.

---

## 3. Installation

The repository does not state a minimum Python version. It was tested here
with **Python 3.11**.

```
pip install -r requirements.txt
```

This installs `numpy` (2.0 or newer), `scipy` (1.13 or newer) and
`matplotlib` (3.8.4 or newer), with their own dependencies. The test here
installed numpy 2.4.6, scipy 1.17.1 and matplotlib 3.11.2.

**Why these minimum versions.** The code calls `numpy.trapezoid`, which
exists only in NumPy 2.0 and newer (NumPy 1.26.4 has only the older
`numpy.trapz`). SciPy 1.13 and Matplotlib 3.8.4 are the first releases
that work with NumPy 2. Older SciPy releases and Matplotlib 3.7.3 to 3.8.3
declare `numpy<2` or a similar limit; Matplotlib 3.7.0 to 3.7.2 declare
none but fail to import with NumPy 2.

**Checked with the minimum versions.** On 30 September 2026 the code was
also run with exactly numpy 2.0.0, scipy 1.13.0 and matplotlib 3.8.4
(Python 3.11): the three test files printed exactly their reference output,
`bash run_all.sh` ran through to `ALL DONE`, and `bcs_gap_table.npz`,
`nonlinear_click.json`, both click-trace files and `numbers.json` came out
identical to the stored files. `limits.json` differed in 18 of its numbers
(fitted quantities of the knee checks) by at most 7.2e-8 relative, which
does not change any rounded value in `numbers.json`.

---

## 4. Quick start: three ways to use the code

Run all commands from the repository folder. Times below were measured on
a shared 2-core Linux machine that was also running other jobs.

### Way A: check that everything works (about 20 seconds)

```
python tests/test_abs.py
python tests/test_noise.py
python tests/test_general.py
```

Each file ends with `ALL ... TESTS PASSED`. The printed lines should match
`tests/reference_output_abs.txt`, `reference_output_noise.txt` and
`reference_output_general.txt`; in the test here they matched exactly.

### Way B: make the fast results and all figures (about 25 seconds)

The slow results (`limits.json`, `nonlinear_click.json`, the click traces)
and the gap table are already stored in `data/`, so only the fast steps are
needed:

```
mkdir -p data figures
python scripts/exp_universal.py
python scripts/exp_design.py
python scripts/exp_matched_points.py
python scripts/exp_calorimetry.py
python scripts/fig_device.py
python scripts/fig1.py
python scripts/fig2.py
python scripts/fig3.py
python scripts/fig4.py
python scripts/figS1.py
python scripts/figS2.py
python scripts/make_numbers.py
```

The figures appear in `figures/` as PDF files. The figure scripts do not
create `figures/` themselves, so the `mkdir` line is needed. In the test here, this sequence
took 23 seconds from a fresh copy, and `data/numbers.json` came out
identical to the stored file.

### Way C: recompute everything from scratch (about 10 minutes)

```
bash run_all.sh
```

`run_all.sh` runs the gap table, the three test files, the seven
calculation scripts, the seven figure scripts and `make_numbers.py`, in
that order. The earlier README gives the total time as "some tens of
minutes on a laptop". Here the whole run took about 9.5 minutes, most of it
in `exp_limits.py` (8 minutes).

Running Way C overwrites the stored files in `data/`. In the test here,
every stored file (`bcs_gap_table.npz`, `limits.json`,
`nonlinear_click.json`, both click-trace files and `numbers.json`) came out
identical to the stored copy.

---

## 5. The scripts, step by step

| Step | Command | What it does | Time* | Results |
|---|---|---|---|---|
| 0 | `python scripts/make_bcs_table.py` | Solves the weak-coupling BCS gap equation at 270 reduced temperatures T/Tc and stores the gap ratio; prints a check value `u(0.5)=0.956887` | 2 s | `data/bcs_gap_table.npz` |
| 1 | `python tests/test_abs.py`, `test_noise.py`, `test_general.py` | The built-in checks ([Section 9](#9-built-in-checks)) | 11 s, 5 s, 4 s | printed on screen |
| 2 | `python scripts/exp_universal.py` | Responsivity of the short-junction model versus T/Tc* for transparencies 0.1 to 1.0; split into occupation part and gap part; activation energies | 1.5 s | `data/universal.json` |
| 3 | `python scripts/exp_design.py` | For each recipe: noise budget at T = Tc*/6 and 0.1 K; resonator frequency-noise spectra (first four recipes); map of occupation noise versus phonon noise over temperature and tau_A (Ta/Ti/Au); energy resolution versus Tc*/T; phase-bias route; small-area design points | 1 s | `data/design.json` |
| 4 | `python scripts/exp_matched_points.py` | Best-design operating points at 50 and 100 mK with the readout noise floor included; signal-to-noise for a 26 GHz photon and a dark-count estimate | 1 s | `data/matched_points.json` |
| 5 | `python scripts/exp_calorimetry.py` | Energy resolution versus tau_A (numerical matched filter, with and without readout floor, and the analytic occupation-noise-only formula); a simulated detection trace | 1.5 s | `data/calorimetry.json` |
| 6 | `python scripts/exp_limits.py` | Checks of the approximations (the earlier README links this to Sec. VII of the manuscript's supplementary material): pair processes, part of the current above the gap for each recipe, energy-dependent exchange time, non-thermal filling, spread of transparency (random and diffusive "Dorokhov" distributions), time needed to measure tau_A from the noise spectrum, independence from the current calibration | 8 min | `data/limits.json` |
| 7 | `python scripts/exp_nonlinear_click.py` | Nonlinear click Monte Carlo (the earlier README links this to Sec. VIII of the supplementary material): single 26 GHz photon deposited in four designs; nonlinear heating and cooling, level filling and resonator shift; 1000 photon and 1000 dark trials per case; matched-filter detection | 44 s | `data/nonlinear_click.json`, `data/click_traces_3e-08.npz`, `data/click_traces_1e-07.npz` |
| 8 | `python scripts/fig_device.py`, `fig1.py`, `fig2.py`, `fig3.py`, `fig4.py`, `figS1.py`, `figS2.py` | Draw the figures ([Section 6](#6-which-script-makes-which-figure)) | 2 to 6 s each | `figures/*.pdf` |
| 9 | `python scripts/make_numbers.py` | Collects every number quoted in the manuscript (61 values) and writes 51 LaTeX macros | 2 s | `data/numbers.json`, `data/numbers.tex` |

\*Measured on a shared 2-core Linux machine with Python 3.11 while other
jobs were running; a free modern computer is usually faster.

Every random number generator uses a fixed seed, so the results repeat
exactly: the stored data files were reproduced bit for bit in the test
here.

**A few of the results** (from `data/numbers.json`, key names in brackets;
all depend on the assumed tau_A). Most values in `numbers.json` are computed
from the result files, but three are typed directly into `make_numbers.py`:
`ystar` (2.40), `LxiMoRe` (0.43) and `dark50` (the text
`3\times10^{-11}`; `exp_matched_points.py` prints 2.78e-11 per second for
the 50 mK, tau_A = 1 microsecond point).

| Result | Value |
|---|---|
| Ta/Ti/Au junction at 100 mK, tau_A = 1 microsecond, 1 s of averaging: temperature resolution limited by occupation noise [`dTA100`] versus by phonon noise [`dTph100`] | 50 microkelvin versus 2.1 microkelvin |
| Design condition for the level energy, used in the scripts as 2.3994 kB T; the stored value [`ystar`] is typed into `make_numbers.py`, not computed. Numerical optimum found by the design scan in `exp_design.py` [`optDkT`] | 2.40 kB T; 2.39 kB T |
| Energy resolution (as a frequency, sigma_E/h) of the best design with tau_A = 1 microsecond and readout floor, at 100 mK [`sigEmatched100`] and 50 mK [`sigEmatched50`] | 7.9 GHz and 1.53 GHz |
| Predicted resonator frequency noise at 1 Hz for Ta/Ti/Au at Tc*/6 [`SnuTa`], versus the assumed quantum-limited readout floor [`SnuFloor`] | 198 versus 26 Hz per root hertz |
| Nonlinear click Monte Carlo, 50 mK, 5.3 x 0.5 micrometre channel, tau_A = 30 ns, activated exchange: signal-to-noise [`nlSNR30`] and detection efficiency [`nlEff30`] | 9.1 and 1.000 |

---

## 6. Which script makes which figure

Figure numbers are those written in each script's description.

| Figure | Content | Data from | Drawn by |
|---|---|---|---|
| Fig. 1 (a) | 3D schematic of the Ta/Ti/Au graphene junction, with an inset of the Andreev levels | none (drawing only) | `fig_device.py` -> `figures/fig_device.pdf` |
| Fig. 1 (b, c) | Responsivity versus T/Tc* for several transparencies; occupation versus gap part | step 2 | `fig1.py` -> `figures/fig1.pdf` |
| Fig. 2 | Predicted frequency-noise spectra; temperature-resolution budgets (Ta/Ti/Au); regime map | step 3 (panel b is computed when drawn) | `fig2.py` -> `figures/fig2.pdf` |
| Fig. 3 | Energy resolution versus Tc*/T (design law); phase-bias route; recipes and designs | steps 3, 4 | `fig3.py` -> `figures/fig3.pdf` |
| Fig. 4 | Energy resolution versus tau_A; running matched-filter output; score histograms | steps 5, 7 | `fig4.py` -> `figures/fig4.pdf` |
| Supplementary Fig. S1 | Monte Carlo check of the telegraph-noise spectrum | computed when drawn | `figS1.py` -> `figures/figS1.pdf` |
| Supplementary Fig. S2 | Monte Carlo signal-to-noise versus tau_A; raw single-shot records | step 7 | `figS2.py` -> `figures/figS2.pdf` |

`fig1.py` and `fig_device.py` make separate PDF files; according to
`fig1.py`, they are combined into one figure in the LaTeX source of the
manuscript, which is not in this repository.

---

## 7. The Python modules

| File | What it contains |
|---|---|
| `src/constants.py` | Physical constants (SI units, CODATA 2018); graphene Fermi velocity 1.0e6 m/s; BCS gap ratio 1.764 |
| `src/materials.py` | The six junction recipes (`RECIPES`); gate-voltage-to-density conversion; graphene Fermi energy, density of states, electron heat capacity; electron-phonon cooling power and thermal conductance; number of transverse modes |
| `src/abs_model.py` | Exact Andreev-level energies of a ballistic channel of finite length with one barrier (`abs_energies`); contribution of the states above the gap; `JunctionModel` (free energy, current-phase relation, critical current, calibration); BCS gap versus temperature (`gap_bcs`); split of the current into below-gap and above-gap parts (`continuum_share`) |
| `src/short_junction.py` | `ShortJunction`: closed-form short-junction levels E = Delta sqrt(1 - tau sin^2(phi/2)); responsivities; occupation-noise sums, Andreev heat capacity C_A, and achieved versus bound temperature resolution |
| `src/finiteL.py` | `FiniteLJunction`: the same sums from the exact finite-length levels; ratio of achieved to bound resolution |
| `src/sensor_limits.py` | `SensorBudget`: electron heat capacity, thermal conductance and thermal time; occupation and phonon temperature resolution; Josephson inductance and participation in the readout resonator; frequency-noise spectrum; numerical matched-filter and analytic energy resolution |
| `src/montecarlo.py` | Random two-state (telegraph) occupation traces and a one-sided spectrum estimate |
| `src/noise_general.py` | Four-state model of one level (both spins) with single and pair processes (exact spectrum); activated exchange time tau(E) = tau0 exp[(Delta - E)/kB T]; penalty for a non-thermal filling q |
| `scripts/figstyle.py` | Shared APS-style figure settings (STIX fonts, Okabe-Ito colour-blind-safe colours, single- and double-column widths) |

---

## 8. Where the numbers come from

**Measured junction recipes.** `src/materials.py` takes these from Table I
of Jung *et al.* (cited in the code as arXiv:2503.06850, published as Phys.
Rev. Applied 26, 014078 (2026)). Contact stacks with film thicknesses in nm:

| Recipe (label) | Tc* (K) | Coherence length (um) | L (um) | W (um) | Back gate (V) | Transparency | Ic at 20 mK (uA) | Rn (ohm) |
|---|---|---|---|---|---|---|---|---|
| Ta(10)/Ti(60)/Au(5) (Ta/Ti/Au) | 0.57 | 7.2 | 0.2 | 5.3 | 30 | 0.30 | 1.08 | 43.3 |
| Ti(6)/Al(60)/Au(5) (Ti/Al/Au) | 0.75 | 5.1 | 0.2 | 1.8 | 30 | 0.78 | 2.11 | 42.0 |
| Ti(6)/Al(70) (Ti/Al(thin)) | 0.99 | 4.4 | 0.3 | 2.9 | 20 | 0.53 | 2.11 | 64.9 |
| Ti(6)/Al(200) (Ti/Al(thick)) | 1.17 | 3.6 | 0.3 | 2.9 | 45 | 0.42 | 3.17 | 48.2 |
| Ti(6)/Nb(5)/NbN(50) (Ti/Nb/NbN) | 2.5 | 0.38 | 0.1 | 1.6 | 30 | 0.58 | 0.985 | 72.9 |
| MoRe(50) (MoRe) | 7.4 | 0.46 | 0.2 | 1.8 | 30 | 0.27 | 3.13 | 133.0 |

**Derived values.** The induced gap is Delta* = 1.764 kB Tc*. Carrier
density comes from a parallel-plate model of the 280 nm SiO2 gate
(relative permittivity 3.9) with a charge-neutrality offset of -2 V,
"chosen so that Vbg = 20 V gives n ~ 1.7e16 m^-2 as quoted by Jung et al."
(comment in `materials.py`). The number of orbital transverse modes is
kF W / pi (`n_modes`). `ShortJunction` counts two channels per orbital mode
(valley x orbital, spin-degenerate): Nch = 2 kF W / pi. `JunctionModel`
counts orbital modes and puts the valley factor 2 into the current instead.
This is why `calibration.json` lists 938 channels for Ta/Ti/Au while
`test_abs.py` prints N=469.

**Heat flow.** The electron-phonon cooling power is Sigma A (Te^3 - Tp^3)
with Sigma = 2.0 W m^-2 K^-3. The code comment cites Lee *et al.*, Nature
586, 42 (2020) (https://doi.org/10.1038/s41586-020-2752-4), with measured
values 2.1 to 3.3.

**One calibration per recipe.** The computed current is multiplied by one
overall factor so that the critical current at 20 mK equals the measured
value (`calibrate()` in `abs_model.py` and `short_junction.py`).
`exp_limits.py` checks that the bound-saturation ratio and C_A do not
depend on this factor.

**Assumed values** (not measured for these devices; set in the scripts):

| Value | Used | Where |
|---|---|---|
| Exchange time tau_A | 1 microsecond by default; scanned from 1 ns to 1 ms (regime map) and from 10^-8.5 to 10^-4.5 s, about 3 ns to 32 us (Fig. 4a) | `exp_design.py`, `exp_calorimetry.py` |
| Readout resonator | 6 GHz, external inductance 2 nH | `SensorBudget` defaults |
| Readout noise floor | quantum-limited phase readout with 30 photons and linewidth 1 MHz | `exp_design.py`, `exp_matched_points.py`, `exp_calorimetry.py`, `exp_nonlinear_click.py` |
| Photon energy | h x 26 GHz | calorimetry and click scripts |
| Best-design channel | W x L = 1.0 x 0.1 um, transparency 0.3, Tc* chosen so that Delta*(T) = 2.3994 kB T | `exp_matched_points.py`, `exp_calorimetry.py` |
| Click Monte Carlo designs | W x L at 50 mK: 1.0 x 0.1, 5.3 x 0.5, 5.3 x 1.5 um; at 100 mK: 5.3 x 0.5 um | `exp_nonlinear_click.py` |

---

## 9. Built-in checks

- **`tests/test_abs.py`** checks the level solver against exact results:
  the short-junction formula (to 1e-12 of the gap), the Kulik levels of a
  fully transparent long channel (to 1e-10 rad), the ballistic result
  e Ic Rn = pi Delta (to 5e-3), and that the calibration factor of the
  first two recipes lies between 0.1 and 10.
- **`tests/test_noise.py`** checks that a short junction with one
  transparency reaches the bound exactly (to 1e-9) for current and
  inductance readout (first three recipes); that the analytic occupation
  responsivity matches a numerical derivative (to 1e-5); and, with a
  Monte Carlo of telegraph noise, the zero-frequency noise (within 10%),
  the half-height point of the spectrum (ratio 0.5 within 0.12), and the
  rule "variance of a t-second average = S/(2t)" (within 25%).
- **`tests/test_general.py`** checks the four-state model: the
  single-process limit, that pair processes do not change the equilibrium
  variance and only shorten the effective correlation time, conservation of
  probability, the zero-frequency noise against a four-state Monte Carlo
  (within 25%), and the bound (with exact equality when the coupling is
  proportional to the level energy) over 20 random level sets.
- **`scripts/make_bcs_table.py`** prints `u(0.5)` next to the expected
  value 0.956887.
- **`scripts/exp_calorimetry.py`** prints the ratio of the numerical to
  the analytic energy resolution at tau_A = 1 microsecond (1.12 here,
  stored as `anchorRatio`).
- **`scripts/exp_nonlinear_click.py`** prints the energy-conservation error
  of the heat-balance integrator (0 to 2e-16 here).
- **`scripts/exp_limits.py`** records, among others, that the
  single-process case is the worst case (`pair_worst_is_singles`) and that
  the results do not depend on the calibration (`scale_independent_ok`);
  both are `true` in the stored `limits.json`.

---

## 10. Notes on the calculations

- **Units.** All calculations use SI units. Results are stored in
  convenient units, marked in the key names: energy resolution as sigma_E/h
  in GHz (`sigE_GHz`), temperature resolution in microkelvin (`_uK`), times
  in ns (`_ns`), heat capacities in units of kB (`_kB`).
- **Noise convention.** Spectra are one-sided. The variance of a
  t-second average of white noise with spectrum S is S/(2t); this
  convention is checked by Monte Carlo in `test_noise.py`.
- **Two junction models.** `ShortJunction` uses the closed-form levels of
  a short junction and is used for the noise budgets. `JunctionModel` and
  `FiniteLJunction` use the exact levels at the real channel length. The
  largest length-to-coherence-length ratio in the recipe set is 0.43
  (MoRe). The finite-length sums in `finiteL.py` leave out the states above
  the gap; `exp_limits.py` estimates their share for each recipe. For MoRe
  this share is large, and `make_numbers.py` reports the largest share over
  the other five recipes only (`contShareMaxPct`).
- **Gap versus temperature.** `gap_bcs` uses the tabulated BCS solution
  above T/Tc = 0.08 and a low-temperature asymptotic formula below it (the
  code calls it the "exact low-T asymptote"), so
  `data/bcs_gap_table.npz` must exist (it is stored).
- **Two values of the BCS ratio.** `constants.py` uses 1.764 (for the
  recipes); `abs_model.py` and the design scripts use 1.7639.
- **Grazing modes.** Transverse modes with cos(theta) < 0.05 are left out
  of `JunctionModel` (`COS_TH_MIN`).
- **Unused constants.** `constants.py` also defines a graphene mass
  density, sound velocity and deformation potential; no script uses them.
- **`data/calibration.json`** lists, for each recipe, the number of
  channels, the carrier density and two calibration factors (`s_short`,
  `s_full`). No script in the repository writes or reads it. Recomputed
  here, the channel numbers, densities and `s_short` match, and `s_full`
  matches `JunctionModel(recipe).calibrate()` to within 0.001 except for
  MoRe (0.501 stored, 0.493 recomputed; the value depends on the phase
  grid).
- **Docstring versus code in the click Monte Carlo.** The description at
  the top of `exp_nonlinear_click.py` says 400 photon and 400 dark trials per
  scenario; the code uses 1000 of each (`ntr = 1000`), and the Fig. 4 legend
  also says 1000 trials.
- **Checks that are claimed but not in the tests.** The earlier README
  listed "free-energy continuity across bound-state exit into the
  continuum" among the testbench pass criteria, and the docstrings in
  `abs_model.py` say this sign convention is "verified in the testbench".
  No test in `tests/` checks it. Likewise, `gap_bcs` says it is "validated
  against both asymptotes in the testbench", but no test calls `gap_bcs`
  directly (only the check value printed by `make_bcs_table.py` exists).
- **LaTeX macros.** `make_numbers.py` writes the macros to
  `data/numbers.tex`, next to `data/numbers.json`. That file is not stored
  in the repository (it is listed in `.gitignore`); the manuscript source
  itself is not included.

---

## 11. Version history

The repository has no releases or version tags. All code was added on
21 August 2026 (git history). On 30 September 2026 `make_numbers.py` was
changed to write `data/numbers.tex` instead of `paper/numbers.tex`, and
`requirements.txt` was corrected to NumPy 2.0 or newer. Details are in
[CHANGELOG.md](CHANGELOG.md).

---

## 12. How to cite

Please cite the manuscript. GitHub also shows a **"Cite this repository"**
button in the right-hand column, which reads `CITATION.cff`.

> T. M. Mahim, A. S. M. Mohsin, and M. M. Rahman, "Andreev occupation
> noise sets the sensitivity limit of proximity Josephson thermal
> detectors" (manuscript).

If you use the junction parameters, please also cite the measurements:

> W. Jung *et al.*, "Engineering Andreev bound states for thermal sensing
> in proximity Josephson junctions", Phys. Rev. Applied 26, 014078 (2026),
> https://doi.org/10.1103/9lsg-mdb8

---

## 13. License and contact

Code: Apache License 2.0 (see `LICENSE`). The earlier README states
"Data and archived outputs: CC-BY-4.0 (see Zenodo record)"; the Zenodo
record is not identified in this repository.

Questions and bug reports: please open an issue on this repository, or
contact Tanvir M. Mahim, BRAC University (tanvir.mahim@bracu.ac.bd).
