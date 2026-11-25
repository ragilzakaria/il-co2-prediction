"""Roadmap step 1.4: do the corrected DFT descriptors help, especially for unseen anions?

Trains a Random Forest and a histogram gradient-boosting model on the referenced dataset
(12,297 points, 109 ILs) with six feature sets:

  TP                   temperature and pressure only
  TP+counts            + elemental atom counts of both ions (ion identity, no QM)
  TP+DFT_thesis        + the 10 thesis DFT descriptors (wrong HOMO lines, broken PF6, x-dipoles)
  TP+DFT_r2scan3c      + the 10 corrected r2SCAN-3c descriptors
  TP+DFT_thesis+counts       the thesis feature set
  TP+DFT_r2scan3c+counts     the same with corrected descriptors

and evaluates each under four 5-fold schemes (point-wise, leave-IL-out, leave-cation-out,
leave-anion-out; the folds are identical for every feature set). For leave-anion-out it also
reports the out-of-fold error per anion.

Writes paper/data/descriptor_comparison.csv, paper/data/descriptor_comparison_per_anion.csv
and images/checks/descriptor_comparison.png. Takes about an hour on a laptop.

Run from the repository root:
    python scripts/compare_descriptor_sets.py
"""
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor, RandomForestRegressor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_validation import DFT, SEED, STATE, TARGET, cross_validate, load, splitters  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
THESIS = os.path.join(ROOT, "data", "processed", "final_referenced_model_dataset.csv")
CORRECTED = os.path.join(ROOT, "data", "processed", "final_referenced_model_dataset_r2scan3c.csv")
OUT = os.path.join(ROOT, "paper", "data")
FIG = os.path.join(ROOT, "images", "checks", "descriptor_comparison.png")


def model_fns():
    return {
        "Random Forest": lambda: RandomForestRegressor(n_estimators=300, random_state=SEED, n_jobs=-1),
        "Gradient Boosting": lambda: HistGradientBoostingRegressor(max_iter=500, learning_rate=0.05,
                                                                   random_state=SEED),
    }


def main():
    th, _, counts = load(THESIS)
    co, _, _ = load(CORRECTED)
    keys = ["cation", "anion", "Temperature_K", "Pressure_kPa", TARGET]
    assert (th[keys].values == co[keys].values).all(), "the two tables must hold the same points in the same order"

    df = th[keys + ["IL"] + counts].copy()
    for c in DFT:
        df[f"{c}__thesis"] = th[c].to_numpy()
        df[f"{c}__r2scan3c"] = co[c].to_numpy()
    dft_th = [f"{c}__thesis" for c in DFT]
    dft_co = [f"{c}__r2scan3c" for c in DFT]
    sets = {
        "TP": STATE,
        "TP+counts": STATE + counts,
        "TP+DFT_thesis": STATE + dft_th,
        "TP+DFT_r2scan3c": STATE + dft_co,
        "TP+DFT_thesis+counts": STATE + dft_th + counts,
        "TP+DFT_r2scan3c+counts": STATE + dft_co + counts,
    }
    y = df[TARGET].to_numpy()
    folds = splitters(df)
    print(f"{len(df)} points, {df['IL'].nunique()} ILs, {len(counts)} varying atom-count features")

    rows, per_anion = [], []
    t0 = time.time()
    for mname, fn in model_fns().items():
        for sname, cols in sets.items():
            X = df[cols].to_numpy(dtype=float)
            for scheme, f in folds.items():
                oof, fold_summary, pooled = cross_validate(fn, X, y, f)
                rows.append(dict(model=mname, features=sname, n_features=len(cols), scheme=scheme,
                                 R2_mean=fold_summary["R2"]["mean"], R2_sd=fold_summary["R2"]["sd"],
                                 MAE_mean=fold_summary["MAE"]["mean"], MAE_sd=fold_summary["MAE"]["sd"],
                                 R2_pooled=pooled["R2"], MAE_pooled=pooled["MAE"], RMSE_pooled=pooled["RMSE"]))
                if scheme == "leave_anion_out":
                    err = pd.DataFrame({"anion": df["anion"], "abs_err": np.abs(oof - y)})
                    for an, g in err.groupby("anion"):
                        per_anion.append(dict(model=mname, features=sname, anion=an, n_points=len(g),
                                              MAE=g["abs_err"].mean()))
                print(f"[{time.time() - t0:6.0f}s] {mname:17s} {sname:24s} {scheme:17s} "
                      f"R2 {fold_summary['R2']['mean']:+.3f} ± {fold_summary['R2']['sd']:.3f}  "
                      f"MAE {fold_summary['MAE']['mean']:.4f}", flush=True)

    res = pd.DataFrame(rows)
    os.makedirs(OUT, exist_ok=True)
    res.to_csv(os.path.join(OUT, "descriptor_comparison.csv"), index=False, float_format="%.4f")
    pd.DataFrame(per_anion).to_csv(os.path.join(OUT, "descriptor_comparison_per_anion.csv"), index=False,
                                   float_format="%.4f")

    # Figure: fold-mean R2 per scheme for each feature set (Random Forest)
    rf = res[res["model"] == "Random Forest"]
    schemes = list(folds)
    fig, ax = plt.subplots(figsize=(10, 4.5))
    w = 0.13
    for i, s in enumerate(sets):
        sub = rf[rf["features"] == s].set_index("scheme").loc[schemes]
        ax.bar(np.arange(len(schemes)) + (i - 2.5) * w, sub["R2_mean"], w, yerr=sub["R2_sd"], label=s, capsize=2)
    ax.set_xticks(np.arange(len(schemes)))
    ax.set_xticklabels([s.replace("_", " ") for s in schemes])
    ax.set_ylabel("R² (mean of 5 folds)")
    ax.set_ylim(0, 1)
    ax.set_title("Random Forest, referenced dataset: feature sets under four validation schemes")
    ax.legend(fontsize=7, ncol=3, frameon=False, loc="lower left")
    fig.tight_layout()
    os.makedirs(os.path.dirname(FIG), exist_ok=True)
    fig.savefig(FIG, dpi=160)

    print("\nFold-mean R2 (Random Forest):")
    print(rf.pivot(index="features", columns="scheme", values="R2_mean")[schemes].loc[list(sets)].round(3).to_string())
    print("\nFold-mean R2 (Gradient Boosting):")
    gb = res[res["model"] == "Gradient Boosting"]
    print(gb.pivot(index="features", columns="scheme", values="R2_mean")[schemes].loc[list(sets)].round(3).to_string())
    print(f"\nTotal {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
