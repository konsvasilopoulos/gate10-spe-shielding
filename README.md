# gate10-spe-shielding

[![tests](https://github.com/konsvasilopoulos/gate10-spe-shielding/actions/workflows/tests.yml/badge.svg)](https://github.com/konsvasilopoulos/gate10-spe-shielding/actions/workflows/tests.yml)

Simulation and analysis code for the article

> K. Vasilopoulos ([ORCID 0009-0001-7245-1404](https://orcid.org/0009-0001-7245-1404)), K. Chatzipapas, *The spectral slope of solar particle events sets the effectiveness of shielding materials*, manuscript in preparation for Life Sciences in Space Research (2026).

The code uses [GATE 10](https://github.com/OpenGATE/opengate) (opengate 10.0.0, Geant4 11.2.p01) to transport solar particle event (SPE) protons through a shield slab into a 30 g/cm² water phantom, and computes absorbed dose and dose equivalent (ICRP 60). It covers:

- validation with 150 MeV protons against NIST PSTAR and measured reaction cross sections;
- the NASA code-comparison benchmark (Webber 1956 SPE, 20 g/cm² Al; Heinbockel et al., NASA/TP-2009-215560);
- a comparison of nine shield materials (polyethylene, lithium hydride, water, PMMA, Al, Ti, Fe, Cu, Pb) at 5, 10 and 20 g/cm² for five SPE spectra: the classical Webber 1956 and King 1972 fits and the Koldobskiy et al. (2021) reconstructions of February 1956, August 1972 and September–October 1989;
- an analytic CSDA transmission-threshold model and an EM-only/full-physics decomposition;
- a response matrix (25 energy windows, 30 MeV–5 GeV) that folds any SPE spectrum without new simulations, applied to the 58 ground-level enhancement (GLE) spectra reconstructed by Koldobskiy et al. (2021), with a threshold rule that gives the dose ratio of a material to Al from the slope of the spectrum near the transmission threshold;
- a comparison with the deterministic code HZETRN (NASA OLTARIS) for the same slabs and four spectra, including the Tylka et al. reconstruction of October 1989 and the LaRC fits of February 1956 and August 1972 (Wilson & Denn 1976);
- sensitivity tests: isotropic incidence, SPE α particles, particle-specific RBE (NCRP 132), the neutron share of the BFO dose, a cortical-bone layer in front of the water and borated polyethylene.

## A note for GATE 10 users: `spectrum_histogram` weights

In opengate 10.0.0 the `spectrum_histogram` energy source treats `spectrum_weights` as a **probability density per unit energy**: the weights are multiplied by the bin width before the cumulative distribution is built (`GateGenericSource.cpp`). If you supply the fluence *inside* each bin with non-uniform (e.g. logarithmic) bins, the sampled spectrum is tilted by a factor proportional to the bin width. `run_spectrum.py` therefore divides the per-bin fluence by the bin width (`weights_as_density: true` in `setup.json`). `analysis/check_spectrum_sampling.py` demonstrates the effect in about one minute.

## Installation

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

The versions are the ones used for the paper (Python 3.9.6). The scripts also run with the latest release, opengate 10.1.1 (Geant4 11.4.2, Python 3.10-3.14), which changes the filter API; `run_release_check.sh` compares the two releases (absorbed dose within 1.2 %, ratios between materials within 1 %, dose equivalent 0-5.5 % higher with 10.1.1). The first `import opengate` downloads the Geant4 data sets. All commands below are run from the repository root.

## Tests

```bash
pip install pytest
pytest
```

The tests (a few seconds, no GATE needed) check the spectrum normalisations, the response matrix and its folding against the stored validation results and the NASA-limit table (Table 4 of the paper) and the comparison with OLTARIS (Table 2). They run on every push (GitHub Actions).

## Repository layout

```
run_material.py            150 MeV monoenergetic runs (validation geometry)
run_spectrum.py            one SPE energy group on a shield + water phantom (main simulation)
spectra.py                 SPE spectra: Webber 1956, King 1972, Koldobskiy et al. 2021, auxiliary log-flat
run_*.sh                   batch scripts (each resumes if interrupted; caffeinate on macOS)
analyze_spectrum.py        combines energy groups into the event dose; plane fluence of all species
analyze_dose_eq.py         dose equivalent from the LET actors (ICRP 60 Q(L))
compare_materials_spe.py   per-spectrum material table (D, H at skin, 5 g/cm2, 10 g/cm2, mean 0-30)
compare_nasa.py            comparison with the NASA benchmark
validation/                150 MeV analysis (PSTAR comparison, seed statistics, entrance dose) and results
analysis/                  paper analyses, figures and their CSV outputs
analysis/response_matrix.npz   the response matrix (see below)
benchmark/                 digitised NASA/TP-2009-215560 curves (WebPlotDigitizer)
oltaris/                   HZETRN results from OLTARIS (TARIS 5.01) for Al, PE and Pb slabs (see analysis/oltaris_compare.py)
pstar/                     NIST PSTAR tables for the Geant4 materials (lithium hydride: scaled from PE, analysis/pstar_scaled.py)
figures/                   paper figures (PDF)
tests/                     pytest checks of the spectra, the response matrix and the paper tables
```

## Reproducing the results

Run times are for 4 cores of an Apple M5; set `JOBS` to the number of fast cores.

| Step | Command | Time |
|---|---|---|
| 150 MeV validation | `bash run_validation_150MeV.sh`, then `python validation/compare_pstar.py`, `python validation/seed_stats.py runs` | ~1 h |
| Hadronic origin, nuclear attenuation | `bash run_hadronic_origin.sh`, `python analysis/hadronic_origin.py`, `python analysis/nuclear_attenuation.py` | ~10 min |
| Koldobskiy 1956/1972 materials + NASA benchmark | `bash run_set1_kold_nasa.sh` | 7–10 h |
| King 1972 and Webber 1956 materials | `bash run_set2_king_webber.sh` | 8–10 h |
| Koldobskiy 1989 materials + QBBC test | `bash run_set3_1989_qbbc.sh` | 3.5–4 h |
| Response matrix, 7 materials (525 runs) | `bash run_set4_response.sh`, then `python analysis/response_matrix.py build` | 2.5–3 h |
| EM-only decomposition | `bash run_em_only.sh`, then `python analysis/decompose_em_nuclear.py` | ~2 h |
| 30 and 40 g/cm2 shields (response matrix + direct check) | `bash run_set5_thick.sh`, then `python analysis/response_matrix.py build_thick`, `python analysis/nasa_limits_thick.py` | ~1-2.5 h |
| SPE alpha particles (Al, PE, Pb, 5-20 g/cm2) | `bash run_set6_alpha.sh`, then `python analysis/alpha_contribution.py` | ~2.2 h |
| Neutron share of the BFO dose (Al, PE, Pb at 20 g/cm2, with and without neutron hadronic processes) | `bash run_neutron_share.sh`, then `python analysis/neutron_share.py` | ~45 min |
| Lithium hydride and water (direct runs, 5 spectra) | `bash run_set7_lih_water.sh`, then `compare_materials_spe.py` for each spectrum | ~3 h |
| Lithium hydride and water: response matrix 5–40 g/cm2, α particles; cortical-bone layer | `bash run_set8_lih_water_rm.sh`, then `python analysis/response_matrix.py add`, `python analysis/alpha_contribution.py`, `python analysis/bone_bpe.py` | ~3.3 h |
| Isotropic-incidence test (Al, PE, Pb at 20 g/cm2) | `bash run_iso_test.sh`, then `python analysis/iso_test.py` | ~1.5 h |
| Lateral-size and layered-shield tests | `bash run_quick_tests.sh`, then `python analysis/quick_tests.py` | ~35 min |
| Comparisons with external references repeated with GATE 10.1.1 | `bash run_release_benchmark.sh`, then `python analysis/release_benchmark.py` | ~2.5 h |
| Release check with GATE 10.1.1 / Geant4 11.4.2 (separate environment) | `bash run_release_check.sh`, then `python validation/seed_stats.py runs_gate101` (results: `analysis/release_check.csv`) | ~1.5 h |

After each material set: `python compare_materials_spe.py spe_runs/materials_<spectrum>` (Webber: `spe_runs/materials`).

The borated-polyethylene test (5 % boron by mass, 1.04 g/cm³, I = 58.35 eV from Bragg additivity; Al, PE and borated PE at 20 g/cm² for the Koldobskiy 1956 and 1972 spectra) needs a user-defined material with its own ionisation potential, which opengate 10.0.0 does not support. It was run with opengate 10.1.1 / Geant4 11.4.2 and a custom-material definition that is not part of this repository; only ratios within that release are used (`analysis/bone_bpe.csv`).

Analyses and figures (from the CSV summaries):

```bash
python analysis/pstar_scaled.py            # PSTAR-format table for lithium hydride (scaled from PE)
python analysis/thresholds_pstar.py        # analytic threshold model vs GATE
python analysis/uncertainty_budget.py
python analysis/nasa_limits.py             # NASA-STD-3001 30-day limits (direct runs, 5-20 g/cm2)
python analysis/nasa_limits_thick.py       # same with the response matrix, 5-40 g/cm2 (Table 4)
python analysis/fold_spectra.py crossval   # group folding; predict1989 / check1989 / scan
python analysis/response_matrix.py validate
python analysis/response_matrix.py scan    # all 58 GLEs -> analysis/rm_gle_scan.csv
python analysis/gle_rule.py
python analysis/spectral_uncertainty.py    # 68 % envelope of the reconstructed spectra
python analysis/iso_test.py               # isotropic vs normal incidence (after run_iso_test.sh)
python analysis/alpha_contribution.py     # SPE alpha particles relative to protons (after run_set6_alpha.sh)
python analysis/oltaris_compare.py        # GATE vs HZETRN (OLTARIS) and Tylka vs Koldobskiy for October 1989
python analysis/oltaris_larc.py           # OLTARIS LaRC spectra (Wilson & Denn 1976) for 1956 and 1972
python analysis/rbe_check.py              # particle-specific RBE from the HZETRN particle doses
python analysis/bone_bpe.py               # bone layer and borated PE (ratios to Al / PE)
python analysis/gle_rule_td.py            # threshold rule on the Tylka & Dietrich Band fits of 59 GLEs (Raukunen et al. 2018)
python analysis/energy_contrib.py         # primary energies that give the BFO dose and the Al-PE difference
python analysis/make_figures.py; python analysis/make_fig_geometry.py; python analysis/fig_response.py
```

The CSV files in `analysis/` and `validation/results/` are the outputs used in the paper.

## Using the response matrix

`analysis/response_matrix.npz` holds, for each material `M` (POLYETHYLENE, LITHIUM_HYDRIDE, WATER, PLEXIGLASS, Al, Ti, Fe, Cu, Pb) and areal density `T` (5, 10, 20, 30, 40 g/cm²):

- `D_M_T`, `H_M_T`: arrays of shape (25, 300) with the absorbed dose (Gy cm²) and dose equivalent (Sv cm²) per unit incident fluence, for 25 logarithmic energy windows (`edges`, MeV) and 300 depth voxels of 1 mm in water (`z`, g/cm²);
- normal incidence on an infinite slab, scored in the central 30 × 30 cm² of a wide water phantom.

```bash
python analysis/response_matrix.py fold 5          # GLE 5 (Koldobskiy et al. 2021 parameters)
python analysis/response_matrix.py fold kold1989   # any spectrum defined in spectra.py
```

To fold your own spectrum, define its integral fluence Φ(>E) in `spectra.py` and call `fold`. For the five directly simulated spectra the folded mean dose agrees with the GATE runs within ±1.5 %; the mean dose equivalent agrees within ±1.5 % except for King 1972 behind 10–20 g/cm², where the direct runs are statistically limited (−3 % to +3 %). `python analysis/response_matrix.py validate` writes the full comparison.

## Data

The raw dose and LET maps (MetaImage) of all runs are archived on Zenodo: DOI to be added at publication. The phase-space ROOT files are not archived because of their size; the plane fluences derived from them are included in the Zenodo archive.

## How to cite

If you use this code or the response matrix, please cite the article:

> K. Vasilopoulos, K. Chatzipapas, The spectral slope of solar particle events sets the effectiveness of shielding materials, *Life Sciences in Space Research* (manuscript in preparation, 2026).

and the archived software version:

> K. Vasilopoulos, gate10-spe-shielding (v1.0.0), Zenodo, doi: *to be added at release*.

`CITATION.cff` provides the citation in machine-readable form (GitHub: "Cite this repository").

## Funding

This work was supported by the Andreas Mentzelopoulos Foundation through a doctoral scholarship.

## Licence

The code is released under the MIT licence.

The spectral parameters come from Koldobskiy et al. (2021), A&A 647, A132; the stopping-power tables from NIST PSTAR; the benchmark curves from Heinbockel et al. (2009), NASA/TP-2009-215560.
