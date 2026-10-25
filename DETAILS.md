# Andreev Occupation Noise: details

This file has the full notes for this repository. For the overview, installation and quick start,
see [README.md](README.md).

## Contents

- [Extra notes for Sections 1 to 4](#extra-notes-for-sections-1-to-4)
  - [More on the introduction: an earlier title](#more-on-the-introduction-an-earlier-title)
  - [More on Section 1: everything the code computes](#more-on-section-1-everything-the-code-computes)
  - [More on Section 2: the full file tree and the stored data](#more-on-section-2-the-full-file-tree-and-the-stored-data)
  - [More on Section 3: versions and checks](#more-on-section-3-versions-and-checks)
  - [More on Section 4: run times and test results](#more-on-section-4-run-times-and-test-results)
- [5. The scripts, step by step](#5-the-scripts-step-by-step)
- [6. Which script makes which figure](#6-which-script-makes-which-figure)
- [7. The Python modules](#7-the-python-modules)
- [8. Where the numbers come from](#8-where-the-numbers-come-from)
- [9. Built-in checks](#9-built-in-checks)
- [10. Notes on the calculations](#10-notes-on-the-calculations)
- [11. Version history](#11-version-history)

---

## Extra notes for Sections 1 to 4

### More on the introduction: an earlier title

(An earlier draft of the manuscript was titled "Andreev
occupation noise sets the sensitivity limit of proximity Josephson thermal
detectors".)

### More on Section 1: everything the code computes

Starting from the measured properties of six graphene junctions (Jung *et
al.*, cited at the top of [README.md](README.md)), the code computes:

- the energy levels and supercurrent of each junction, with an exact
  finite-length solver and a simpler short-junction formula;
- the temperature response and the occupation noise of the levels, and the
  floor on the temperature resolution set by that noise (it follows from
  equilibrium statistics and the Cauchy-Schwarz inequality);
- when this floor exceeds the usual floor from heat exchange between the
  graphene electrons and the lattice (called phonon thermal-fluctuation
  noise, "phonon TFN", in the code);
- the frequency noise that the occupation noise would produce in a readout
  resonator, from which tau_A can be measured, and a comparison with the
  1/f noise reported for graphene junctions;
- the best induced gap (levels at about 2.40 kB T) and how it shifts when
  the exchange is carried by thermal quasiparticles;
- a two-body thermal model in which the electrons and the Andreev levels
  exchange quasiparticles and energy, with a closed-form energy resolution;
- a self-consistent simulation of a single 26 GHz photon: the full
  nonlinear mean response, noise that follows the hot electrons, and an
  ideal or a resonator readout;
- checks of the approximations (pair processes, level-dependent exchange
  times, non-thermal filling, spread of contact transparency, the part of
  the current carried above the gap, the electron-phonon exponent).

### More on Section 2: the full file tree and the stored data

```
andreev-occupation-noise/
|-- README.md             overview, installation, quick start
|-- DETAILS.md            these notes
|-- CHANGELOG.md          what changed
|-- CITATION.cff          citation details (drives the "Cite this repository" button)
|-- LICENSE               Apache-2.0 license
|-- requirements.txt      Python packages to install
|-- run_all.sh            runs every step in order (bash)
|-- src/
|   |-- constants.py      physical constants (CODATA 2018) and graphene values
|   |-- materials.py      the six measured recipes; graphene heat capacity and cooling
|   |-- abs_model.py      exact finite-length Andreev-level solver, free energy, BCS gap
|   |-- short_junction.py closed-form short-junction model: response, noise sums, floor
|   |-- finiteL.py        how far the finite-length levels fall short of the floor
|   |-- sensor_limits.py  noise budgets, resonator frequency noise, energy resolution
|   |-- coupled.py        two-body thermal model (electrons + Andreev levels)
|   |-- click_sim.py      self-consistent photon click simulation, resonator readout
|   |-- montecarlo.py     random telegraph simulation and spectrum estimate
|   `-- noise_general.py  pair processes, level-dependent exchange, non-thermal filling
|-- tests/
|   |-- test_abs.py            level solver against exact limits
|   |-- test_noise.py          noise model (floor, response, Monte Carlo)
|   |-- test_general.py        generalized noise model
|   |-- test_coupled.py        two-body model and click simulation
|   `-- reference_output_*.txt expected screen output of each test file
|-- scripts/
|   |-- make_bcs_table.py        temperature dependence of the superconducting gap
|   |-- exp_universal.py         response versus temperature for many transparencies
|   |-- exp_design.py            recipe budgets, noise spectra, regime map, design scans
|   |-- exp_matched_points.py    matched-gap designs at 50 and 100 mK
|   |-- exp_calorimetry.py       energy resolution versus tau_A (Fig. 4a)
|   |-- exp_limits.py            checks of the model approximations
|   |-- exp_revision.py          analyses added in the revision (Section 5)
|   |-- exp_click_coupled.py     self-consistent single-photon click simulation
|   |-- fig_device.py            device schematic (Fig. 1a)
|   |-- fig1.py ... fig4.py      main-text figures
|   |-- figS1.py, figS2.py       supplementary figures
|   |-- figstyle.py              shared figure style
|   `-- make_numbers.py          every number quoted in the manuscript
`-- data/                  results stored in the repository (see below)
```

**Stored in `data/`:** `bcs_gap_table.npz` (gap table),
`calibration.json` (see [Section 10](#10-notes-on-the-calculations)), the
outputs of the two slow steps (`limits.json`; `click_coupled.json`,
`click_coupled_T_3e-08.npz` and `click_coupled_C_3e-08.npz`) and
`numbers.json` (the quoted numbers). The other result files
(`universal.json`, `design.json`, `matched_points.json`,
`calorimetry.json`, `revision.json`) are not stored. You can make them in
a few seconds (Way B in [README.md](README.md)). `make_numbers.py` also
writes the LaTeX macros `numbers.tex` into `data/`. That file is not stored
either (it is listed in `.gitignore`). Figures go to `figures/`, which is
also not stored.

### More on Section 3: versions and checks

These notes go with the install command in
[Section 3 of README.md](README.md#3-installation).

**Why these minimum versions.** The code calls `numpy.trapezoid`, which
exists only in NumPy 2.0 and newer (NumPy 1.26.4 has only the older
`numpy.trapz`). SciPy 1.13 and Matplotlib 3.8.4 are the first releases
that work with NumPy 2. Older SciPy releases and Matplotlib 3.7.3 to 3.8.3
declare `numpy<2` or a similar limit. Matplotlib 3.7.0 to 3.7.2 declare no
limit, but they fail to import with NumPy 2.

**Checked with the minimum versions.** On 30 September 2026 I ran the
code of that time with exactly numpy 2.0.0, scipy 1.13.0 and
matplotlib 3.8.4 (Python 3.11), and every stored file was reproduced. I
have not repeated this check for the October 2026 revision; the revision
was tested with Python 3.11 and current releases of the three packages.

### More on Section 4: run times and test results

These notes go with the commands in
[Section 4 of README.md](README.md#4-quick-start-three-ways-to-use-the-code).
I measured the times on a shared 2-core Linux machine that was also running
other jobs.

**Way A** takes about 20 seconds. In my test the printed lines of all four
test files matched their reference output exactly.

**Way B** takes well under a minute. The figure scripts do not create
`figures/` themselves, so you need the `mkdir` line.

**Way C** took about 14 minutes on 25 October 2026, most of it in
`exp_click_coupled.py` (8.6 minutes) and `exp_limits.py` (4.4 minutes).
Running Way C overwrites the stored files in `data/`. In my run from a
fresh copy every result file came out identical to the stored copy, except
`limits.json`, whose fitted values can differ at the 1e-8 relative level
(see [Section 10](#10-notes-on-the-calculations)).

---

## 5. The scripts, step by step

| Step | Command | What it does | Time* | Results |
|---|---|---|---|---|
| 0 | `python scripts/make_bcs_table.py` | Solves the weak-coupling BCS gap equation at 270 reduced temperatures and stores the gap ratio; prints `u(0.5)=0.956887` | 1 s | `data/bcs_gap_table.npz` |
| 1 | `python tests/test_abs.py`, `test_noise.py`, `test_general.py`, `test_coupled.py` | The built-in checks ([Section 9](#9-built-in-checks)) | 10 s, 4 s, 3 s, 3 s | printed on screen |
| 2 | `python scripts/exp_universal.py` | Response of the short-junction model versus T/Tc* for transparencies 0.1 to 1; occupation part and gap part | 1 s | `data/universal.json` |
| 3 | `python scripts/exp_design.py` | Budgets of the six recipes at Tc*/6 and 0.1 K; resonator frequency-noise spectra; regime map; energy resolution versus Tc*/T; phase-bias route | 1 s | `data/design.json` |
| 4 | `python scripts/exp_matched_points.py` | Matched-gap designs: 5.3 x 0.2 um at 100 mK (read at 6 GHz) and 5.3 x 0.5 um at 50 mK (read at 3 GHz), with an ideal readout and with a 1 MHz resonator | 1 s | `data/matched_points.json` |
| 5 | `python scripts/exp_calorimetry.py` | Energy resolution versus tau_A for the 50 mK design (Fig. 4a): closed form, one-way form, readout floor, 1 and 10 MHz resonators | 1 s | `data/calorimetry.json` |
| 6 | `python scripts/exp_limits.py` | Checks of the approximations: pair processes, states above the gap, level-dependent exchange, non-thermal filling, spread of transparency, measurement of tau_A from the spectrum, independence from the calibration | 4 to 8 min | `data/limits.json` |
| 7 | `python scripts/exp_revision.py` | Electron-phonon exponents 3 to 6; activated exchange and the gap optimum for exchange by thermal quasiparticles; two-body energy resolution of recipes and matched designs; comparison with 1/f noise; transition frequencies; energy balance of the photon response | 2 s | `data/revision.json` |
| 8 | `python scripts/exp_click_coupled.py` | Self-consistent click simulation of a 26 GHz photon: breakdown from the earlier model to the full model, resonator readout (1 to 300 MHz, 0 and 14 added quanta), index versus tau_A for four designs, induced-gap scan, electron-phonon exponents | 6 to 10 min | `data/click_coupled.json`, `data/click_coupled_{T,C}_3e-08.npz` |
| 9 | `python scripts/fig_device.py`, `fig1.py` to `fig4.py`, `figS1.py`, `figS2.py` | Draw the figures ([Section 6](#6-which-script-makes-which-figure)) | 1 to 5 s each | `figures/*.pdf` |
| 10 | `python scripts/make_numbers.py` | Collects every number quoted in the manuscript (90 values) and writes 77 LaTeX macros | 1 s | `data/numbers.json`, `data/numbers.tex` |

\*Measured on a shared 2-core Linux machine with Python 3.11 while other
jobs were running.

Every random number generator uses a fixed seed, so the results repeat
exactly.

**A few of the results** are listed below. They come from
`data/numbers.json` (key names in brackets), and all of them depend on the
assumed tau_A. `ystar` (2.40) is typed into `make_numbers.py`; the other
values are computed from the result files.

| Result | Value |
|---|---|
| Ta/Ti/Au at 100 mK, tau_A = 1 us, 1 s of averaging: occupation-noise floor versus phonon floor [`dTA100`, `dTph100`] | 50 uK versus 2.1 uK |
| Exchange time above which the occupation floor dominates (Ta/Ti/Au, 100 mK) [`tauAstar`] | 1.8 ns |
| Matched level energy (used in the scripts as 2.3994 kB T) [`ystar`] | 2.40 kB T |
| Energy resolution sigma_E/h of the matched designs at tau_A = 1 us, ideal readout / 1 MHz resonator [`sigEmatched100`, `sigEmatched100k`, `sigEmatched50`, `sigEmatched50k`] | 100 mK: 41.1 / 131 GHz; 50 mK: 18.9 / 52 GHz |
| Frequency noise at 1 Hz, Ta/Ti/Au at Tc*/6, versus the readout floor [`SnuTa`, `SnuFloor`] | 198 versus 26 Hz per root hertz |
| Click simulation, 5.3 x 0.5 um at 50 mK, tau_A = 30 ns, ideal readout: Kumar-Carroll index, scenario T / C [`kcFullT`, `kcFullC`] | 4.1 / 3.3 |
| Same with a 50 MHz resonator and 30 photons [`kcResFifty`] | 1.6 |
| Best index found with the T^3 cooling law [`kcBest`] | 5.0 |

---

## 6. Which script makes which figure

| Figure | Content | Data from | Drawn by |
|---|---|---|---|
| Fig. 1 (a) | Device schematic of the Ta/Ti/Au junction, with an inset of the Andreev levels | none (drawing only) | `fig_device.py` -> `figures/fig_device.pdf` |
| Fig. 1 (b, c) | Response versus T/Tc*; occupation and gap parts | step 2 | `fig1.py` -> `figures/fig1.pdf` |
| Fig. 2 | Frequency-noise spectra; temperature-resolution budgets; regime map | step 3 | `fig2.py` -> `figures/fig2.pdf` |
| Fig. 3 | Energy resolution versus Tc*/T; phase-bias route; recipes and matched designs | steps 3, 4, 7 | `fig3.py` -> `figures/fig3.pdf` |
| Fig. 4 | Energy resolution versus tau_A; click score histograms; breakdown of the index | steps 5, 8 | `fig4.py` -> `figures/fig4.pdf` |
| Supplementary Fig. S1 | Monte Carlo check of the telegraph-noise spectrum | computed when drawn | `figS1.py` -> `figures/figS1.pdf` |
| Supplementary Fig. S2 | Index versus tau_A; index versus resonator linewidth; raw records | step 8 | `figS2.py` -> `figures/figS2.pdf` |

`fig1.py` and `fig_device.py` make separate PDF files, which the
manuscript places one above the other as Fig. 1.

---

## 7. The Python modules

| File | What it contains |
|---|---|
| `src/constants.py` | Physical constants (SI units, CODATA 2018); graphene Fermi velocity 1.0e6 m/s; BCS gap ratio 1.764 |
| `src/materials.py` | The six junction recipes (`RECIPES`); gate-voltage-to-density conversion; graphene Fermi energy, density of states, electron heat capacity; electron-phonon cooling power and thermal conductance; number of transverse modes |
| `src/abs_model.py` | Exact Andreev levels of a finite-length ballistic channel with one barrier (`abs_energies`); states above the gap; `JunctionModel` (free energy, current-phase relation, critical current, calibration); BCS gap versus temperature (`gap_bcs`); split of the current into bound and continuum parts (`continuum_share`) |
| `src/short_junction.py` | `ShortJunction`: closed-form levels E = Delta sqrt(1 - tau sin^2(phi/2)); responses; occupation-noise sums, C_A, and achieved versus floor |
| `src/finiteL.py` | `FiniteLJunction`: the same sums from the exact finite-length levels |
| `src/sensor_limits.py` | `SensorBudget`: electron heat capacity, thermal conductance and time; temperature resolutions; Josephson inductance and participation in the readout resonator; frequency-noise spectrum; energy resolution |
| `src/coupled.py` | Two-body thermal model: noise-equivalent power, closed-form energy resolution sigma_E = 2 sqrt(kB T^2 C_e) [x(1+x)]^(1/4) with x = G_ep tau_A / C_A, numerical filter with readout floor and resonator linewidth, corner time of the noise spectrum |
| `src/click_sim.py` | `ClickModel`: nonlinear mean response of electrons and Andreev occupations after a photon; fluctuations in the linear-noise approximation (exact matrix exponential per step); records with the nonlinear frequency map. `resonator_phase`: exact cavity response. `matched_filter_stats`: whitened matched filter and Kumar-Carroll index |
| `src/montecarlo.py` | Random two-state (telegraph) occupation traces and a one-sided spectrum estimate |
| `src/noise_general.py` | Four-state model of one level (both spins) with single and pair processes; activated exchange time tau(E) = tau0 exp[(Delta - E)/kB T]; penalty for a non-thermal filling q |
| `scripts/figstyle.py` | Shared figure settings (STIX fonts, Okabe-Ito color-blind-safe colors, figure widths) |

---

## 8. Where the numbers come from

**Measured junction recipes.** `src/materials.py` takes these from Table I
of Jung *et al.* (Phys. Rev. Applied 26, 014078 (2026)). The contact stacks
are given with film thicknesses in nm:

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
(relative permittivity 3.9) with a charge-neutrality offset of -2 V. A
comment in `materials.py` says this offset was "chosen so that Vbg = 20 V
gives n ~ 1.7e16 m^-2 as quoted by Jung et al." `ShortJunction` counts two
channels per orbital mode (valley x orbital, spin-degenerate):
Nch = 2 kF W / pi. `JunctionModel` counts orbital modes and puts the valley
factor 2 into the current instead. This is why `calibration.json` lists
938 channels for Ta/Ti/Au while `test_abs.py` prints N=469.

**Heat flow.** The electron-phonon cooling power is Sigma A (Te^3 - T0^3)
with Sigma = 2.0 W m^-2 K^-3. This is the law with which Lee *et al.*,
Nature 586, 42 (2020) (https://doi.org/10.1038/s41586-020-2752-4) fitted
a graphene Josephson bolometer operated at 0.19 K. They set the exponent
to 3 and fitted Sigma between 2.04 and 3.30 W m^-2 K^-3; 2.0 is close to
the low end of that range. With Sigma = 2.0, the code gives a thermal time
of 0.8 ns at 0.19 K and n = 2e16 m^-2, the same order as the 0.6 ns
estimated by Lee *et al.* `exp_revision.py` and
`exp_click_coupled.py` also use exponents 4 to 6, with Sigma chosen so that
the thermal conductance at 0.19 K is unchanged.

**One calibration per recipe.** The computed current is multiplied by one
overall factor, so that the critical current at 20 mK equals the measured
value. `exp_limits.py` checks that the bound-saturation ratio and C_A do
not depend on this factor.

**Assumed values.** These are not measured for these devices. They are set
in the scripts:

| Value | Used | Where |
|---|---|---|
| Exchange time tau_A | 1 us by default; scanned from about 3 ns to 30 us (Fig. 4a) and from 10 ns to 1 us (click simulation) | most scripts |
| Readout resonator | 6 GHz for the recipes and the 100 mK design; 3 GHz for the 50 mK designs, below their pair transition 2 Delta*/h = 5 GHz; external inductance 2 nH | `SensorBudget`, `exp_click_coupled.py` |
| Readout noise | quantum-limited phase readout with 30 photons and 1 MHz linewidth; resonator readout with 1 to 300 MHz linewidth and 0 or 14 added quanta | `exp_design.py`, `exp_matched_points.py`, `exp_click_coupled.py` |
| Photon energy | h x 26 GHz | calorimetry and click scripts |
| Matched designs | transparency 0.3, Tc* chosen so that Delta*(T) = 2.3994 kB T; 5.3 x 0.2 um at 100 mK, 5.3 x 0.5 um at 50 mK | `exp_matched_points.py`, `exp_revision.py`, `exp_click_coupled.py` |
| Exchange scenarios | C: constant tau_A; T: tau = tau_A exp[-(Delta*/kB)(1/T0 - 1/Te)], at least 1 ns | `click_sim.py` |

---

## 9. Built-in checks

- **`tests/test_abs.py`** checks the level solver against exact results:
  the short-junction formula (to 1e-12 of the gap), the Kulik levels of a
  fully transparent long channel (to 1e-10 rad) and e Ic Rn = pi Delta (to
  5e-3). It checks that the calibration factors are of order one, that the
  free energy stays continuous when a level leaves the gap (to 0.4%), and
  that the tabulated BCS gap matches its two limits (1% and 0.3%).
- **`tests/test_noise.py`** checks that a single-transparency short
  junction reaches the floor exactly (to 1e-9), that the analytic response
  matches a numerical derivative (to 1e-5), and, with a telegraph Monte
  Carlo, the zero-frequency noise (10%), the half-height point (0.12) and
  the rule "variance of a t-second average = S/(2t)" (25%).
- **`tests/test_general.py`** checks the four-state model: the
  single-process limit, that pair processes leave the variance unchanged,
  conservation of probability, the zero-frequency noise against a
  four-state Monte Carlo (25%), and the floor over 20 random level sets.
- **`tests/test_coupled.py`** checks the closed-form energy resolution
  against the numerical filter (better than 4e-4; tolerance 2e-3), that the
  one-way and two-body models give the same noise-equivalent power, that
  the linear-noise model reproduces the equilibrium variances (to 1e-10)
  and the closed form in its linear limit (5e-4; tolerance 3e-3), and that
  the corner time equals tau_A for Ta/Ti/Au (0.2%; tolerance 0.5%).
- **`scripts/make_bcs_table.py`** prints `u(0.5)` next to the expected
  value 0.956887.
- **`scripts/exp_calorimetry.py`** prints the agreement of the numerical
  filter with the closed form (3e-4).
- **`scripts/exp_revision.py`** prints the energy balance of the mean
  photon response (residual 3e-5 or less).
- **`scripts/exp_limits.py`** records that the single-process case is the
  worst case (`pair_worst_is_singles`) and that the results do not depend
  on the calibration (`scale_independent_ok`); both are `true` in the
  stored `limits.json`.

---

## 10. Notes on the calculations

- **Units.** All calculations use SI units. Stored results mark their units
  in the key names: sigma_E/h in GHz (`_GHz`), microkelvin (`_uK`), ns
  (`_ns`), heat capacities in units of kB (`_kB`).
- **Noise convention.** Spectra are one-sided. The variance of a t-second
  average of white noise with spectrum S is S/(2t).
- **Two junction models.** `ShortJunction` uses the closed-form levels and
  is used for the noise budgets. `JunctionModel` and `FiniteLJunction` use
  the exact levels at the real channel length. The finite-length sums leave
  out the states above the gap; `exp_limits.py` estimates their share. For
  MoRe this share is large, so `make_numbers.py` reports the largest share
  over the other five recipes only (`contShareMaxPct`).
- **Linear-noise approximation.** The click simulation keeps the full
  nonlinear mean response and treats fluctuations to linear order around
  it. This needs an electron heat capacity of order 10 kB or more. Designs
  with C_e < 1 kB (for example 1 x 0.1 um) are outside its validity and are
  not simulated.
- **Record length.** Click records span ten slow relaxation times
  (C_e + C_A)/G_ep, at most 12 us. For electron-phonon exponents 5 and 6
  this cap is shorter than ten relaxation times, so those values are less
  certain.
- **Energy-resolution filter.** Since October 2026 the numerical filter in
  `sensor_limits.py` integrates to 1000/(2 pi min(tau_A, tau_th)). The
  earlier limit of 10/(2 pi min(...)) truncated the integral and overstated
  sigma_E by up to about 10%.
- **Gap versus temperature.** `gap_bcs` uses the tabulated BCS solution
  above T/Tc = 0.08 and a low-temperature asymptotic formula below it, so
  `data/bcs_gap_table.npz` must exist (it is stored).
- **Two values of the BCS ratio.** `constants.py` uses 1.764 (for the
  recipes); `abs_model.py` and the design scripts use 1.7639.
- **Grazing modes.** Transverse modes with cos(theta) < 0.05 are left out
  of `JunctionModel` (`COS_TH_MIN`).
- **Unused constants.** `constants.py` also defines a graphene mass
  density, sound velocity and deformation potential; no script uses them.
- **`data/calibration.json`** lists, for each recipe, the number of
  channels, the carrier density and two calibration factors (`s_short`,
  `s_full`). No script writes or reads it.
- **Reproducibility.** `limits.json` holds fitted quantities that can
  differ at the 1e-8 relative level between runs; no rounded number in
  `numbers.json` changes.
- **LaTeX macros.** `make_numbers.py` writes the macros to
  `data/numbers.tex`, next to `data/numbers.json`. That file is not stored
  in the repository; the manuscript source itself is not included.

---

## 11. Version history

All code was first added on 21 August 2026. On 30 September and
2 October 2026 the documentation, `make_numbers.py`, `requirements.txt`,
the tests and two calibration factors were corrected. On 25 October 2026
the code was extended for the revised manuscript: a two-body thermal model,
a self-consistent photon click simulation that replaces
`exp_nonlinear_click.py`, the analyses requested in review, a new test file,
and a fix of the energy-resolution filter. The archived version for the
revised manuscript is v1.13 on Zenodo (doi:10.5281/zenodo.23648106). It
differs from v1.11 only in the layout of the figures (v1.12) and in
comments and documentation about the cooling law of Lee *et al.* (v1.13). The
details are in [CHANGELOG.md](CHANGELOG.md).
