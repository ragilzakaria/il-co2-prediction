#!/usr/bin/env python3
"""
Figures for the current-state manuscript.

  fig_dataset     coverage of the training set (data/final_model_dataset_v3.csv)
  fig_importance  permutation importances from the thesis model comparison
                  (data/full_feature_importance_comparison.csv, notebook 17)

The parity plot (figures/fig_parity_thesis.png) is copied unchanged from
images/results/ (notebook 17).

Run from the repository root:
    python paper/code/make_figures.py
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "../.."))
FIG = os.path.abspath(os.path.join(HERE, "../figures"))

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "savefig.dpi": 300, "savefig.bbox": "tight",
})
BLUE, ORANGE, GREY = "#2a6f97", "#d1495b", "#8d99ae"


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"))
    plt.close(fig)
    print(f"[fig] {name}")


def fig_dataset():
    df = pd.read_csv(os.path.join(ROOT, "data", "final_model_dataset_v3.csv"))
    df["IL"] = df["cation"] + df["anion"]
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.6), gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axes[0]
    sc = ax.scatter(df["Pressure_kPa"] / 1000, df["CO2_solubility"], c=df["Temperature_K"],
                    s=3, alpha=0.6, cmap="viridis", linewidths=0, rasterized=True)
    ax.set_xscale("log")
    ax.set_xlabel("Pressure (MPa)")
    ax.set_ylabel(r"CO$_2$ mole fraction $x_{\mathrm{CO_2}}$")
    ax.set_title("(a) Measurements", loc="left")
    cb = fig.colorbar(sc, ax=ax, pad=0.02)
    cb.set_label("T (K)")
    ax = axes[1]
    counts = df["IL"].value_counts().to_numpy()
    ax.bar(np.arange(len(counts)), counts, width=1.0, color=GREY)
    ax.set_yscale("log")
    ax.set_xlabel("Ionic liquid (ranked)")
    ax.set_ylabel("Data points")
    ax.set_title(f"(b) Points per IL ({len(counts)} ILs)", loc="left")
    fig.tight_layout()
    save(fig, "fig_dataset")


def fig_importance():
    imp = pd.read_csv(os.path.join(ROOT, "data", "full_feature_importance_comparison.csv"),
                      index_col=0)
    imp = imp[~imp.index.str.contains("_num_")]
    models = ["Random Forest", "Gradient Boosting"]
    order = imp["Random Forest"].sort_values().index
    pretty = {"Pressure_kPa": "Pressure", "Temperature_K": "Temperature"}
    names = [pretty.get(f, f.replace("Energy_Hartree", "total energy").replace("_Eh", "")
                        .replace("_Debye", "").replace("HOMO_LUMO_Gap", "HOMO-LUMO gap")
                        .replace("Dipole_Moment", "dipole").replace("_", " ")) for f in order]
    fig, ax = plt.subplots(figsize=(5.6, 3.4))
    y = np.arange(len(order))
    h = 0.4
    for k, (m, c) in enumerate(zip(models, (BLUE, ORANGE))):
        vals = imp.loc[order, m].clip(lower=1e-5)
        ax.barh(y + (k - 0.5) * h, vals, h, color=c, label=m)
    ax.set_yticks(y, names)
    ax.set_xscale("log")
    ax.set_xlabel(r"Permutation importance (drop in test-set $R^2$)")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    save(fig, "fig_importance")


def main():
    os.makedirs(FIG, exist_ok=True)
    fig_dataset()
    fig_importance()


if __name__ == "__main__":
    main()
