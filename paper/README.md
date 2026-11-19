# Manuscript: ML prediction of CO₂ solubility in ionic liquids

**Title**: *Machine-Learning Prediction of CO₂ Solubility in Ionic Liquids from DFT Ion Descriptors*

**Format**: one-column manuscript (IEEEtran, 11pt, A4). It is not yet a submission.

## Layout

```
paper/
├── main.tex                  # IEEEtran, one column, 11pt, A4
├── references.bib
├── IEEEtran.cls / IEEEtran.bst
├── sections/
│   ├── 01_introduction.tex
│   ├── 02_data.tex
│   ├── 03_computational.tex  # includes the TikZ workflow figure
│   ├── 04_models.tex
│   ├── 05_results.tex
│   ├── 06_screening.tex
│   ├── 07_issues.tex         # limitations of the current workflow
│   └── 08_summary.tex
├── figures/
│   ├── fig_dataset.pdf       # make_figures.py, from data/final_model_dataset_v3.csv
│   ├── fig_importance.pdf    # make_figures.py, from data/full_feature_importance_comparison.csv
│   ├── fig_parity_thesis.png # copied from images/results/ (notebook 17)
│   └── fig_bmim_oleate.png   # copied from images/data_plots/
├── code/
│   ├── make_figures.py
│   ├── extract_descriptors.py  # corrected parser for calculations/ions/*.out
│   └── check_validation.py     # grouped (leave-IL/cation/anion-out) CV check
└── data/
    ├── dft_descriptors_r2scan3c.csv  # output of extract_descriptors.py
    ├── results.json                  # output of check_validation.py
    └── per_il_error_rf_leave_il_out.csv
```

The numbers in Sections 5 and 6 come from the existing notebooks (`notebooks/15`, `17`, `18`, `20`). The checks in Section 7 come from `code/`.

## Build

```bash
cd paper
pdflatex main && bibtex main && pdflatex main && pdflatex main
```

## Regenerate figures and checks

The scripts need `pandas`, `numpy`, `scikit-learn` and `matplotlib`. Run them from the repository root:

```bash
python paper/code/make_figures.py
python paper/code/extract_descriptors.py
python paper/code/check_validation.py   # about 5 minutes on a laptop
```

## Placeholders to fill

- The citation for the source of `data/rawdata.xlsx` (Tables S5 and S6) in `sections/02_data.tex`.

Placeholders show in red in the PDF.
