# Data sources

Where every data file in this folder comes from and what to cite for it. Checked 26 November 2026:
every DOI below resolves to the work named (Crossref, or DataCite for the Zenodo and TU/e records).

## Experimental CO₂ solubility

| File | Source | Cite |
|---|---|---|
| `raw/rawdata.xlsx` | Supporting Information (Tables S1–S6) of Laakso et al. 2025, CC BY. S6 is the compilation of Song et al. 2020 (10,116 points, 124 ILs) with UNIFAC values; S1 is the same data without duplicates (9,684 points); S2 the removed duplicates; S3 ion names; S4/S5 the 80/20 train/test split; σ-profile columns sc1–sa4 are theirs. | Laakso, Ebrahimpoor Gorji, Alopaeus, Uusi-Kyyny, *Chem. Eng. Sci.* **307**, 121226 (2025), [10.1016/j.ces.2025.121226](https://doi.org/10.1016/j.ces.2025.121226); Song, Shi, Zhang, Zhou, *Chem. Eng. Sci.* **223**, 115752 (2020), [10.1016/j.ces.2020.115752](https://doi.org/10.1016/j.ces.2020.115752) |
| `thesis/cleaned_experimental_data_combined.csv`, `thesis/final_model_dataset*.csv` | Tables S5 + S6 of `rawdata.xlsx`, concatenated (`notebooks/thesis/02`), then merged with descriptors. **Contains duplicates:** S5 is a subset of S6, so 2,013 of the 11,256 rows of `final_model_dataset_v3.csv` repeat an earlier row. | as above |
| `thesis/thesis_point_references.csv` | Per-point trace of the thesis set: 8,778 of 11,256 points (78 %) matched to a publication of the referenced dataset (`scripts/trace_thesis_points.py`). | the matched publication |
| `raw/matched_ilthermo_targets.csv` | The 462 NIST ILThermo data sets of CO₂ solubility for ILs made of our ions. | Dong et al., *J. Chem. Eng. Data* **52**, 1151 (2007), [10.1021/je700171f](https://doi.org/10.1021/je700171f) |
| `raw/ilthermo/<set>.json` | Those 462 data sets as downloaded from the ILThermo API (26 Nov 2026), unchanged. | ILThermo and the publication of each set |
| **`processed/curated_co2_solubility_dataset.csv`** | **The dataset used for the main results.** Built from `raw/ilthermo/` by `scripts/build_curated_dataset.py` (rules in its docstring): IL-rich liquid phase only, units converted with molar masses of the optimised ions, duplicates and republished points removed, inconsistent publications removed. 12,026 points, 110 ILs, 167 publications, 411 data sets; each point has its ILThermo set id, citation, **verified** DOI, method and original unit. | the publication in `ref_citation` / `ref_doi`, and ILThermo |
| `processed/curated_exclusions.csv` | Every excluded data set or point (711 points) with the reason. | — |
| `curation/excluded_sets.csv` | Data sets removed by hand after the consistency check, with the evidence: Kassim 2016, Sistla and Khanna 2014, Lim et al. 2009 (150 points). | — |
| `processed/checks/consistency_by_set.csv`, `consistency_by_point.csv` | Output of `scripts/check_data_consistency.py`: monotonicity and agreement with other publications for every set and point. | — |
| `processed/model_dataset_curated.csv`, `model_dataset_curated_thesisdesc.csv` | The curated dataset merged with the corrected (or the thesis) descriptors and the corrected atom counts (`scripts/merge_referenced_model_dataset.py`). | as above |
| `processed/referenced_co2_solubility_dataset.csv`, `processed/final_referenced_model_dataset*.csv` | First ILThermo build (`scripts/build_referenced_dataset.py`), superseded by the curated dataset: its parser took any column mentioning CO₂ (so CO₂-rich-phase points with x up to 0.9999 entered) and it kept duplicates. Kept to reproduce earlier results; its DOIs are verified. | as above |

## References of the experimental data

| File | What it is |
|---|---|
| `references/doi_verification.csv` | One row per publication of the referenced dataset: the DOI first written by the build script, the verified DOI, the check status and the Crossref record it was checked against. |
| `references/doi_from_search.csv` | The DOIs as first written (top hit of an unchecked Crossref search), kept for the record. |
| `references/doi_overrides.csv` | Four publications set by hand, with the reason (ILThermo errors, a PhD thesis). |
| `references/dataset_bibliography.csv`, `references/dataset_references.bib` | Bibliography of the 171 publications, generated from the verified DOIs. |
| `references/thesis_references_extracted.txt` | Reference list of the thesis (text). |

DOI check result (`scripts/verify_dataset_references.py`): 119 publications correct; 18 pointed to a
Supporting Information file (`.s001`); 33 (2,218 points) pointed to a different work; 77 points had no
DOI. All corrected by `scripts/apply_verified_references.py`. To repeat the check after the dataset
is rebuilt, run both scripts, then `scripts/merge_referenced_model_dataset.py` (with and without
`--thesis`).

## Ions

| File | Source |
|---|---|
| `ions/unique_ions_with_names.csv` | Ion abbreviations and names from Table S3 of `rawdata.xlsx` (Laakso et al. 2025). |
| `ions/ions_with_smiles.csv` | SMILES looked up in PubChem by ion name (`notebooks/thesis/04`; Kim et al., *Nucleic Acids Res.* **51**, D1373 (2023), [10.1093/nar/gkac956](https://doi.org/10.1093/nar/gkac956)). Six lookups returned the wrong species and were corrected by hand (see the paper, Table V). |
| `ions/ions_with_atomic_counts.csv` | Atom counts from the thesis SMILES (RDKit); wrong for the six wrong-species ions. Kept for the thesis tables. |
| `ions/ion_atom_counts_r2scan3c.csv` | Atom counts of the optimised (corrected) structures (`scripts/make_ion_atom_counts.py`); used for the curated model tables. |

## Computed in this project

| File | Source |
|---|---|
| `processed/dft_descriptors_r2scan3c.csv` | This work: ORCA 6.0 / 6.1.1, r²SCAN-3c, outputs in `calculations/ions/` and `calculations/hpc/` (`scripts/extract_descriptors.py`). Cite ORCA (Neese, *WIREs Comput. Mol. Sci.* **12**, e1606 (2022), [10.1002/wcms.1606](https://doi.org/10.1002/wcms.1606)) and r²SCAN-3c (Grimme et al., *J. Chem. Phys.* **154**, 064103 (2021), [10.1063/5.0040021](https://doi.org/10.1063/5.0040021)). |
| `processed/orca_descriptors_individual.csv` | Thesis descriptors (parsed with the faulty thesis parser; kept for comparison only). |
| `processed/new_candidates_descriptors.csv` | Thesis descriptors of the bio-based candidate anions. |
| `processed/checks/` | Comparison of our descriptors with literature values. |
| `processed/md_forcefield_coverage.csv` | Which training ILs have CL&P force-field parameters: A = both ions have a ready CL&P structure file, B = all atom types exist but a new structure file is needed, C = parameters missing. Canongia Lopes and Pádua, *Theor. Chem. Acc.* **131**, 1129 (2012), [10.1007/s00214-012-1129-7](https://doi.org/10.1007/s00214-012-1129-7). |

## Literature reference values for validating calculations

`literature_reference/`: 518 published computed values from seven sources, each row with its source,
DOI and the table it was read from. See `literature_reference/README.md`.
