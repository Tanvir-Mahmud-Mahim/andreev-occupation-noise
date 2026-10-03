# Andreev Occupation Noise: Simulation Code

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://www.apache.org/licenses/LICENSE-2.0)

This repository holds the code for our manuscript **"Andreev-level
occupation noise in graphene Josephson thermal detectors: a sensitivity
floor and how to measure it"**, which I wrote with A. S. M. Mohsin and
M. M. Rahman. An earlier draft had another title (see
[DETAILS.md](DETAILS.md#more-on-the-introduction-an-earlier-title)).

- Repository: https://github.com/Tanvir-Mahmud-Mahim/andreev-occupation-noise
- Measured devices that the model is built on: W. Jung *et al.*,
  "Engineering Andreev bound states for thermal sensing in proximity
  Josephson junctions", Phys. Rev. Applied **26**, 014078 (2026),
  https://doi.org/10.1103/9lsg-mdb8 (preprint: arXiv:2503.06850)

The scripts produce the data files in `data/` and the figures of the
manuscript. According to `make_numbers.py`, they also produce every number
quoted in the manuscript (`data/numbers.json`, see
[Section 5](DETAILS.md#5-the-scripts-step-by-step)).

---

## Contents

1. [The idea in one minute](#1-the-idea-in-one-minute)
2. [What is in this repository](#2-what-is-in-this-repository)
3. [Installation](#3-installation)
4. [Quick start: three ways to use the code](#4-quick-start-three-ways-to-use-the-code)
5. [The scripts, step by step](DETAILS.md#5-the-scripts-step-by-step)
6. [Which script makes which figure](DETAILS.md#6-which-script-makes-which-figure)
7. [The Python modules](DETAILS.md#7-the-python-modules)
8. [Where the numbers come from](DETAILS.md#8-where-the-numbers-come-from)
9. [Built-in checks](DETAILS.md#9-built-in-checks)
10. [Notes on the calculations](DETAILS.md#10-notes-on-the-calculations)
11. [Version history](DETAILS.md#11-version-history)
12. [How to cite](#12-how-to-cite)
13. [License and contact](#13-license-and-contact)

Sections 5 to 11 are in [DETAILS.md](DETAILS.md), together with
[extra notes for Sections 1 to 4](DETAILS.md#extra-notes-for-sections-1-to-4).
A short guide to them is in [More details](#more-details-sections-5-to-11).

---

## 1. The idea in one minute

A **proximity Josephson junction** is a short strip of a normal conductor
(here graphene) between two superconducting contacts. The superconductors
make a supercurrent flow through the graphene. Inside the junction, this
current is carried by discrete energy levels called
**Andreev bound states**. Each level can be empty or can hold stray
excitations (quasiparticles). A filled level carries less current. The
chance that a level is filled depends on temperature. So the current and
the junction's inductance change with temperature, and the junction works
as a thermometer and, in principle, as a photon detector.

In our manuscript we point out that the same filling of the levels also
**flips randomly** in time, even at a perfectly steady temperature. This
random "occupation noise" sets a floor on how finely the junction can
measure temperature. The time over which a level keeps its filling is
called the **exchange time** `tau_A` in the code. It is not known for these
devices, so the code treats it as a free parameter and scans it.

Starting from the measured properties of six graphene junctions (Jung *et
al.*, above), the code computes, among other things:

- the energy levels and supercurrent of each junction;
- the occupation noise of the levels, and a lower bound on the temperature
  resolution set by that noise;
- the frequency noise that the occupation noise would produce in a
  microwave resonator used to read out the junction (a measurable
  prediction);
- energy-resolution budgets for detecting a single 26 GHz microwave photon.

The full list is in
[DETAILS.md](DETAILS.md#more-on-section-1-everything-the-code-computes).

**A few of the results.** They come from `data/numbers.json` (key names in
brackets). All of them depend on the assumed tau_A.

| Result | Value |
|---|---|
| Ta/Ti/Au junction at 100 mK, tau_A = 1 microsecond, 1 s of averaging: temperature resolution limited by occupation noise [`dTA100`] versus by phonon noise [`dTph100`] | 50 microkelvin versus 2.1 microkelvin |
| Predicted resonator frequency noise at 1 Hz for Ta/Ti/Au at Tc*/6 [`SnuTa`], versus the assumed quantum-limited readout floor [`SnuFloor`] | 198 versus 26 Hz per root hertz |

More results are in [Section 5](DETAILS.md#5-the-scripts-step-by-step).

---

## 2. What is in this repository

- **src/**: the Python modules (see
  [Section 7](DETAILS.md#7-the-python-modules)).
- **scripts/**: the calculation and figure scripts (see
  [Section 5](DETAILS.md#5-the-scripts-step-by-step)).
- **tests/**: three test files, with the expected screen output of each.
- **data/**: results stored in the repository.
- `run_all.sh`: runs every step in order (bash).
- At the top level you also find `requirements.txt`, `LICENSE`,
  `CITATION.cff` and [CHANGELOG.md](CHANGELOG.md).

Some result files are not stored. You can make them in a few seconds
(Way B below). Figures go to `figures/`, which is also not stored.

The full annotated file tree, and the list of what is stored in `data/`, are
in [DETAILS.md](DETAILS.md#more-on-section-2-the-full-file-tree-and-the-stored-data).

---

## 3. Installation

The repository does not state a minimum Python version. I tested it with
**Python 3.11**.

```
pip install -r requirements.txt
```

This installs `numpy` (2.0 or newer), `scipy` (1.13 or newer) and
`matplotlib` (3.8.4 or newer), with their own dependencies.

Why these are the minimum versions, and my check with exactly these
versions, are in [DETAILS.md](DETAILS.md#more-on-section-3-versions-and-checks).

---

## 4. Quick start: three ways to use the code

Run all commands from the repository folder. I measured the times below on
a shared 2-core Linux machine that was also running other jobs.

### Way A: check that everything works (about 20 seconds)

```
python tests/test_abs.py
python tests/test_noise.py
python tests/test_general.py
```

Each file ends with `ALL ... TESTS PASSED`. The printed lines should match
`tests/reference_output_abs.txt`, `reference_output_noise.txt` and
`reference_output_general.txt`.

### Way B: make the fast results and all figures (about 25 seconds)

The slow results (`limits.json`, `nonlinear_click.json`, the click traces)
and the gap table are already stored in `data/`. So you only need the fast
steps:

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
create `figures/` themselves, so you need the `mkdir` line.

### Way C: recompute everything from scratch (about 10 minutes)

```
bash run_all.sh
```

`run_all.sh` runs the gap table, the three test files, the seven
calculation scripts, the seven figure scripts and `make_numbers.py`, in
that order. Be aware that running Way C overwrites the stored files in
`data/`.

My test results and run times for each way are in
[DETAILS.md](DETAILS.md#more-on-section-4-run-times-and-test-results).

---

## More details (Sections 5 to 11)

The full notes are in [DETAILS.md](DETAILS.md). Here is what each section
holds:

- [5. The scripts, step by step](DETAILS.md#5-the-scripts-step-by-step): every
  step with its command, run time and result files, and a table of a few
  results.
- [6. Which script makes which figure](DETAILS.md#6-which-script-makes-which-figure):
  each figure, its content, its data, and the script that draws it.
- [7. The Python modules](DETAILS.md#7-the-python-modules): what each file
  in the code contains.
- [8. Where the numbers come from](DETAILS.md#8-where-the-numbers-come-from):
  the measured junction recipes, the derived values, and the assumed values.
- [9. Built-in checks](DETAILS.md#9-built-in-checks): what the three test
  files and the checks inside the scripts cover.
- [10. Notes on the calculations](DETAILS.md#10-notes-on-the-calculations):
  units, conventions, the two junction models, and other notes.
- [11. Version history](DETAILS.md#11-version-history): what changed, by
  date.

---

## 12. How to cite

Please cite our manuscript. On GitHub you will also see a
**"Cite this repository"** button in the right-hand column, which reads
`CITATION.cff`.

> T. M. Mahim, A. S. M. Mohsin, and M. M. Rahman, "Andreev-level
> occupation noise in graphene Josephson thermal detectors: a sensitivity
> floor and how to measure it" (manuscript).

If you use the junction parameters, please also cite the measurements:

> W. Jung *et al.*, "Engineering Andreev bound states for thermal sensing
> in proximity Josephson junctions", Phys. Rev. Applied 26, 014078 (2026),
> https://doi.org/10.1103/9lsg-mdb8

---

## 13. License and contact

Code: Apache License 2.0 (see `LICENSE`). Archived copies of the code
with all generated data are on Zenodo under CC-BY-4.0. You can find every
version listed at https://zenodo.org/records/22040477 (the DOI of the
version used for the manuscript is given in its Data availability
statement).

If you have a question or find a bug, please open an issue on this
repository, or contact me, Tanvir M. Mahim, BRAC University
(tanvir.mahim@bracu.ac.bd).
