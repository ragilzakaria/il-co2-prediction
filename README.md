# CO₂ solubility in ionic liquids from DFT ion descriptors

Machine-learning prediction of CO₂ solubility (mole fraction) in ionic liquids (ILs). Inputs are temperature, pressure and quantum-chemical descriptors of the cation and anion.

The project started as an undergraduate thesis (UGM Chemistry, 2026) and continues as a BRIN internship project.

## Status (November 2026)

- **Thesis model.** A Random Forest trained on 11,256 points for 113 ILs.
  - It scores R² = 0.986 on a random 80/20 split.
  - Grouped validation shows this is optimistic for new ILs: holding out each IL in turn gives R² ≈ 0.94, and holding out each anion gives R² ≈ 0.83 (`paper/data/results.json`).
- **Rebuilt dataset.** 12,297 points for 109 ILs from NIST ILThermo, each point with its citation. Not yet used for training.
- **Fixes in progress.**
  - Six ions were calculated as the wrong species. Their inputs have been rebuilt; the calculations still need to be run.
  - The descriptor parser has been corrected.
  - Next: re-extract the descriptors, then retrain with grouped validation.
- **Next stage** (BRIN HPC): a better QM protocol (conformer search, ion pairs, CO₂ complexes, COSMO-RS), MD of the bulk liquids, and ML with uncertainty estimates.
- **Validation rule.** Every new computed output is checked against published values in `data/literature_reference/`.

## Layout

```
data/
  raw/                    source data as received (rawdata.xlsx = thesis Tables S5/S6; ILThermo matches)
  ions/                   ion master list: names, SMILES, atom counts (84 ions)
  processed/              current pipeline outputs: referenced dataset, descriptors, training table
  thesis/                 frozen thesis datasets (final_model_dataset_v3.csv = thesis training set)
  references/             bibliography of the experimental data
  literature_reference/   published computed values for validating our calculations (see its README)
calculations/
  structures/             3D starting structures (.xyz) for the 84 ions
  ions/                   ORCA r2SCAN-3c inputs and outputs, one per ion
    old_wrong_species/    outputs of the six ions that were built from wrong SMILES (kept for the record)
  candidates/             new anions screened after the thesis (oleate, laurate)
scripts/                  pipeline code: build and merge datasets, rebuild structures, parse ORCA, validate
notebooks/
  thesis/                 the original thesis notebooks 01–20, kept as a record (paths updated)
models/
  thesis/                 thesis Random Forest, scaler and feature names (used by demo/)
paper/                    manuscript (LaTeX); paper/code/make_figures.py makes its figures
slides/                   BRIN sharing-session deck (source + standalone index.html)
demo/                     Flask web demo of the thesis model
images/                   thesis figures, 2D structures and result plots
```

## Pipeline

```
data/raw/rawdata.xlsx ──► notebooks/thesis/01–03 ──► data/ions/ (names, SMILES)
data/ions/ions_with_smiles.csv ──► calculations/structures/ ──► calculations/ions/*.inp ──► ORCA
calculations/ions/*.out ──► scripts/extract_descriptors.py ──► data/processed/dft_descriptors_r2scan3c.csv
data/raw/matched_ilthermo_targets.csv ──► scripts/build_referenced_dataset.py ──► data/processed/referenced_co2_solubility_dataset.csv
                                     ──► scripts/merge_referenced_model_dataset.py ──► data/processed/final_referenced_model_dataset.csv
data/processed/final_referenced_model_dataset.csv ──► scripts/check_validation.py ──► paper/data/results.json
```

## Setup

```bash
pip install -r requirements.txt
```

Run all scripts from the repository root, for example `python scripts/check_validation.py`.

ORCA 6 is needed for the calculations. It is free for academic use and must be installed separately.

## Where new work goes

- **New QM protocols:** a new folder per protocol under `calculations/`, e.g. `calculations/r2scan3c_crest/` or `calculations/ion_pairs/`. Keep `calculations/ions/` as the thesis record.
- **New datasets and descriptor tables:** `data/processed/`. Never overwrite anything in `data/raw/` or `data/thesis/`.
- **Reusable code:** `scripts/`. Notebooks for exploration go in `notebooks/`, numbered from 21.
- **Trained models:** `models/<name>/`. Files over about 50 MB should not be committed. The thesis model is 79 MB and already in git history.
- **New computed outputs:** compare against `data/literature_reference/reference_values.csv` before using them.

## Known issues

- `calculations/ions/BMIM.inp` is a leftover B3LYP/def2-SVP single-point input. It does not match `BMIM.out`, which is the r2SCAN-3c optimisation actually used.
- The thesis notebooks refer to some folders that no longer exist (`orca_calculations/`, `calculations_individual/`). These were mapped to `calculations/candidates/` and `calculations/ions/`. The notebooks have not been re-run since the move.
- `scripts/build_full_refs.py` and `scripts/generate_perfect_pdf.py` use the `workspace/` folder next to this repository, which is not part of the repo.
