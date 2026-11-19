#!/usr/bin/env python3
"""
Experiments for the IL-CO2 solubility paper.

Trains Ridge, Random Forest, Gradient Boosting and an MLP on the referenced
dataset (data/final_referenced_model_dataset.csv) and evaluates them under two
cross-validation schemes:

  * point-wise 5-fold CV (random rows, as in the original thesis split), and
  * leave-IL-out 5-fold CV (GroupKFold, all points of an ionic liquid held out
    together), which measures transfer to unseen ionic liquids, and
  * leave-cation-out and leave-anion-out 5-fold CV, which hold out every IL
    that contains a given ion.

It also runs a feature-set ablation with the Random Forest and computes
permutation importances. Results go to paper/data/.

Run from the repository root:
    python paper/code/check_validation.py
"""

import json
import os
import time

import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, r2_score, root_mean_squared_error
from sklearn.model_selection import GroupKFold, KFold, train_test_split
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "../.."))
DATA_IN = os.path.join(ROOT, "data", "final_referenced_model_dataset.csv")
LEGACY_IN = os.path.join(ROOT, "data", "final_model_dataset_v3.csv")
OUT_DIR = os.path.abspath(os.path.join(HERE, "../data"))

SEED = 42
N_FOLDS = 5
TARGET = "CO2_solubility"
STATE = ["Temperature_K", "Pressure_kPa"]
DFT = [f"{ion}_{d}" for ion in ("cation", "anion") for d in
       ("Energy_Hartree", "HOMO_Eh", "LUMO_Eh", "HOMO_LUMO_Gap_Eh", "Dipole_Moment_Debye")]


def models():
    return {
        "Ridge": make_pipeline(StandardScaler(), Ridge(alpha=1.0)),
        "Random Forest": RandomForestRegressor(n_estimators=300, random_state=SEED, n_jobs=-1),
        "Gradient Boosting": GradientBoostingRegressor(n_estimators=300, random_state=SEED),
        "MLP": make_pipeline(StandardScaler(), MLPRegressor(
            hidden_layer_sizes=(128, 64, 32), early_stopping=True, max_iter=500,
            random_state=SEED)),
    }


def metrics(y, p):
    return {"R2": float(r2_score(y, p)), "MAE": float(mean_absolute_error(y, p)),
            "RMSE": float(root_mean_squared_error(y, p))}


def load(path):
    df = pd.read_csv(path)
    df["IL"] = df["cation"] + df["anion"]
    counts = [c for c in df.columns if "_num_" in c]
    counts_used = [c for c in counts if df[c].nunique() > 1]
    return df, counts, counts_used


def splitters(df):
    il, cat, an = (df[c].to_numpy() for c in ("IL", "cation", "anion"))
    return {
        "pointwise": list(KFold(N_FOLDS, shuffle=True, random_state=SEED).split(il)),
        "leave_IL_out": list(GroupKFold(N_FOLDS).split(il, groups=il)),
        "leave_cation_out": list(GroupKFold(N_FOLDS).split(il, groups=cat)),
        "leave_anion_out": list(GroupKFold(N_FOLDS).split(il, groups=an)),
    }


def cross_validate(model_fn, X, y, folds):
    oof = np.zeros_like(y)
    per_fold = []
    for tr, te in folds:
        m = model_fn()
        m.fit(X[tr], y[tr])
        oof[te] = m.predict(X[te])
        per_fold.append(metrics(y[te], oof[te]))
    summary = {k: {"mean": float(np.mean([f[k] for f in per_fold])),
                   "sd": float(np.std([f[k] for f in per_fold]))} for k in per_fold[0]}
    return oof, summary, metrics(y, oof)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    df, counts_all, counts = load(DATA_IN)
    features = STATE + DFT + counts
    X = df[features].to_numpy(float)
    y = df[TARGET].to_numpy(float)
    folds = splitters(df)
    main_folds = {k: folds[k] for k in ("pointwise", "leave_IL_out", "leave_anion_out")}

    out = {"dataset": {
        "n_points": int(len(df)), "n_ILs": int(df["IL"].nunique()),
        "n_cations": int(df["cation"].nunique()), "n_anions": int(df["anion"].nunique()),
        "n_sources": int(df["ref_doi"].nunique()),
        "year_min": int(df["ref_year"].min()), "year_max": int(df["ref_year"].max()),
        "T_range": [float(df.Temperature_K.min()), float(df.Temperature_K.max())],
        "P_range_kPa": [float(df.Pressure_kPa.min()), float(df.Pressure_kPa.max())],
        "x_range": [float(y.min()), float(y.max())],
        "points_per_IL_median": float(df["IL"].value_counts().median()),
        "imidazolium_fraction": float(df["cation"].str.contains("IM").mean()),
        "n_feature_columns": len(STATE + DFT + counts_all),
        "n_features_varying": len(features),
        "atom_count_features": counts,
    }}

    # 1. Model comparison under both CV schemes ---------------------------------
    oof_frame = df[["cation", "anion", "IL", "Temperature_K", "Pressure_kPa", TARGET]].copy()
    out["models"] = {}
    for name in models():
        out["models"][name] = {}
        for scheme, f in folds.items():
            t0 = time.time()
            oof, summ, pooled = cross_validate(lambda: models()[name], X, y, f)
            out["models"][name][scheme] = {"folds": summ, "pooled": pooled,
                                           "seconds": time.time() - t0}
            oof_frame[f"{name}|{scheme}"] = oof
            print(f"[models] {name:18s} {scheme:13s} R2={summ['R2']['mean']:.3f}"
                  f"±{summ['R2']['sd']:.3f}  MAE={summ['MAE']['mean']:.4f}")
    # Out-of-fold predictions are kept in memory only.

    # 2. Per-IL error under leave-IL-out (Random Forest) -------------------------
    col = "Random Forest|leave_IL_out"
    per_il = (oof_frame.assign(abs_err=(oof_frame[col] - oof_frame[TARGET]).abs())
              .groupby("IL").agg(n=("abs_err", "size"), MAE=("abs_err", "mean"),
                                 x_mean=(TARGET, "mean")).sort_values("MAE", ascending=False))
    per_il.to_csv(os.path.join(OUT_DIR, "per_il_error_rf_leave_il_out.csv"))
    out["per_IL_worst"] = per_il.head(8).reset_index().to_dict(orient="records")
    out["per_IL_MAE_median"] = float(per_il["MAE"].median())

    # 3. Feature-set ablation with Random Forest --------------------------------
    il_onehot = pd.get_dummies(df["IL"], dtype=float)
    sets = {
        "T, P only": df[STATE],
        "T, P + IL identity (one-hot)": pd.concat([df[STATE], il_onehot], axis=1),
        "T, P + atom counts": df[STATE + counts],
        "T, P + DFT descriptors": df[STATE + DFT],
        "T, P + DFT + atom counts (full)": df[features],
    }
    out["ablation"] = {}
    for label, frame in sets.items():
        Xa = frame.to_numpy(float)
        out["ablation"][label] = {}
        for scheme, f in main_folds.items():
            _, summ, _ = cross_validate(
                lambda: RandomForestRegressor(n_estimators=300, random_state=SEED, n_jobs=-1),
                Xa, y, f)
            out["ablation"][label][scheme] = summ
            print(f"[ablation] {label:34s} {scheme:13s} R2={summ['R2']['mean']:.3f}")

    # 4. Permutation importance (Random Forest, 80/20 split as in the thesis) ----
    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.2, random_state=SEED)
    rf = RandomForestRegressor(n_estimators=300, random_state=SEED, n_jobs=-1).fit(Xtr, ytr)
    out["holdout_80_20_rf"] = metrics(yte, rf.predict(Xte))
    pi = permutation_importance(rf, Xte, yte, n_repeats=10, random_state=SEED, n_jobs=-1)
    out["permutation_importance"] = sorted(
        [{"feature": f, "mean": float(m), "sd": float(s)}
         for f, m, s in zip(features, pi.importances_mean, pi.importances_std)],
        key=lambda r: -r["mean"])

    # 5. Legacy thesis dataset, same 80/20 split, for continuity -----------------
    ldf, _, lcounts = load(LEGACY_IN)
    lX = ldf[STATE + DFT + lcounts].to_numpy(float)
    ly = ldf[TARGET].to_numpy(float)
    lXtr, lXte, lytr, lyte = train_test_split(lX, ly, test_size=0.2, random_state=SEED)
    lrf = RandomForestRegressor(n_estimators=300, random_state=SEED, n_jobs=-1).fit(lXtr, lytr)
    lg = ldf["IL"].to_numpy()
    _, lsumm, _ = cross_validate(
        lambda: RandomForestRegressor(n_estimators=300, random_state=SEED, n_jobs=-1),
        lX, ly, list(GroupKFold(N_FOLDS).split(lX, groups=lg)))
    out["legacy_dataset"] = {"n_points": int(len(ldf)), "n_ILs": int(ldf["IL"].nunique()),
                             "rf_holdout_80_20": metrics(lyte, lrf.predict(lXte)),
                             "rf_leave_IL_out": lsumm}
    print("[legacy]", out["legacy_dataset"])

    with open(os.path.join(OUT_DIR, "results.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print("[done] wrote", os.path.join(OUT_DIR, "results.json"))


if __name__ == "__main__":
    main()
