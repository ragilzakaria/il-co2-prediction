# CO₂ solubility in ionic liquids from DFT ion descriptors

Machine-learning prediction of CO₂ solubility (mole fraction) in ionic liquids (ILs). Inputs are temperature, pressure and quantum-chemical descriptors of the cation and anion.

The project started as an undergraduate thesis (UGM Chemistry, 2026) and continues as a BRIN internship project.

## Status (November 2026)

- **Thesis model.** A Random Forest trained on 11,256 points for 113 ILs.
  - It scores R² = 0.986 on a random 80/20 split.
  - Grouped validation shows this is optimistic for new ILs: holding out each IL in turn gives R² ≈ 0.94, and holding out each anion gives R² ≈ 0.83 (`paper/data/results.json`).
- **Rebuilt dataset.** 12,297 points for 109 ILs from NIST ILThermo, each point with its citation. Not yet used for training.
- **Data fixed (25 November 2026).**
  - All 9 broken ions were rerun with ORCA 6.1.1: six built from the wrong species, [PF6] (dissociated), and [HSO4] and [C3F7CO2] (imaginary modes). All are now minima with no imaginary modes.
  - The descriptors were re-extracted with the corrected parser (`data/processed/dft_descriptors_r2scan3c.csv`).
  - A corrected training table has been built: `data/processed/final_referenced_model_dataset_r2scan3c.csv`.
  - Next: retrain with grouped validation (roadmap step 1.4).
- **Next stage** (BRIN HPC): a better QM protocol (conformer search, ion pairs, CO₂ complexes, COSMO-RS), MD of the bulk liquids, and ML with uncertainty estimates.
- **Validation rule.** Every new computed output is checked against published values in `data/literature_reference/`.

**Plan:** see [docs/ROADMAP.md](docs/ROADMAP.md) for the phases, steps and decision points.

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
                                                         ──► scripts/check_descriptors_vs_literature.py (vs data/literature_reference/)
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

- **[PF6]⁻ broke apart in its thesis r2SCAN-3c optimisation. Fixed 25 November 2026.**
  - The thesis run stopped at the 50-cycle limit with one P–F bond at 2.87 Å (PF₅ + F⁻) and a HOMO of +2.4 eV. That affected 1,320 of the 11,256 thesis training points and 1,557 of the 12,297 referenced points.
  - The rerun (ORCA 6.1.1, started from an octahedron) is an intact octahedron: all six P–F bonds at 1.642 Å, no imaginary modes, HOMO −3.79 eV.
  - The old run is in `calculations/ions/old_failed/`.
  - `data/processed/` has not yet been re-extracted with the new [PF6] values.
- **[C3F7CO2]⁻ and [HSO4]⁻ had imaginary frequencies. Fixed.** [HSO4] was re-optimised with its O–H rotated; [C3F7CO2] with TightOpt and DefGrid3. The six wrong-species ions are also done.
  - [HMMIM], [P4,4,4,4] and [C12PhSO3] were optimised on the laptop and their frequencies run on Mahameru (`calculations/hpc/batch01/`), because the macOS ORCA crashes on them. [DBPO4] and [C3F7CO2] ran entirely on Mahameru.
- **No frequency check yet for [N1,8,8,8], [P6,6,6,14] and [oleyl-MIM]** (the thesis frequency step crashed). They are not in the training data. Check them on Mahameru before using them.
- **The thesis HOMO values for anions do not follow the published ordering at all** (Spearman ρ = −0.25). The corrected parser gives ρ = 0.83 once [PF6] is excluded. See `scripts/check_descriptors_vs_literature.py`.
- `calculations/ions/BMIM.inp` is a leftover B3LYP/def2-SVP single-point input. It does not match `BMIM.out`, which is the r2SCAN-3c optimisation actually used.
- The thesis notebooks refer to some folders that no longer exist (`orca_calculations/`, `calculations_individual/`). These were mapped to `calculations/candidates/` and `calculations/ions/`. The notebooks have not been re-run since the move.
- `scripts/build_full_refs.py` and `scripts/generate_perfect_pdf.py` use the `workspace/` folder next to this repository, which is not part of the repo.
