"""Roadmap step 2.3, model side: do wB97X-D3/def2-TZVPD descriptors predict new anions better?

scripts/per_il_errors.py found that the r2SCAN-3c descriptors make leave-one-anion-out predictions
worse for the small oxyanions ([MeSO4], [TFA], [NO3], [MeSO3], [DEPO4]), whose r2SCAN-3c HOMO is
barely bound. This script repeats that test with the same 10 per-ion descriptors (energy, HOMO,
LUMO, gap, dipole of cation and anion) taken from the wB97X-D3/def2-TZVPD single points
(data/processed/dft_descriptors_wb97xd3.csv) instead of r2SCAN-3c, on the curated dataset.

Feature sets (all with T, P and the atom counts of both ions):
  counts      no DFT descriptors
  r2scan3c    + r2SCAN-3c descriptors
  wb97xd3     + wB97X-D3/def2-TZVPD descriptors
Validation: leave-one-anion-out (23 folds), Random Forest and gradient boosting (settings of
scripts/compare_descriptor_sets.py); with --il-out also leave-one-IL-out (110 folds, slower).

Output: paper/data/orbital_level_comparison.csv  (pooled MAE per model, features, scheme)
        paper/data/orbital_level_per_anion.csv   (leave-one-anion-out MAE per held-out anion)

Run from the repository root, after scripts/extract_wb97xd3_descriptors.py:
    python scripts/compare_orbital_levels.py [--il-out] [--allow-missing]
--allow-missing fills ions without a wB97X-D3 result with their r2SCAN-3c values (for testing the
script before the HPC job has finished; the results are then not meaningful).
"""
import os
import sys
import time

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor, RandomForestRegressor
from sklearn.model_selection import LeaveOneGroupOut

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_validation import DFT, SEED, STATE, TARGET, load  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
DATA = os.path.join(ROOT, "data", "processed", "model_dataset_curated.csv")
WB97X = os.path.join(ROOT, "data", "processed", "dft_descriptors_wb97xd3.csv")
OUT = os.path.join(ROOT, "paper", "data")
PROPS = ["Energy_Hartree", "HOMO_Eh", "LUMO_Eh", "HOMO_LUMO_Gap_Eh", "Dipole_Moment_Debye"]


def models():
    return {
        "RF": lambda: RandomForestRegressor(n_estimators=300, random_state=SEED, n_jobs=-1),
        "GB": lambda: HistGradientBoostingRegressor(max_iter=500, learning_rate=0.05, random_state=SEED),
    }


def main():
    allow_missing, il_out = "--allow-missing" in sys.argv, "--il-out" in sys.argv
    df, _, counts = load(DATA)
    wb = pd.read_csv(WB97X).set_index("ion_abbreviation")
    needed = set(df.cation) | set(df.anion)
    missing = sorted(needed - set(wb.index))
    if missing and not allow_missing:
        sys.exit(f"no wB97X-D3 result for {len(missing)} ions: {missing}. Run with --allow-missing to test only.")
    for side in ("cation", "anion"):
        for prop in PROPS:
            col = f"{side}_{prop}"
            df[f"{col}__wb97xd3"] = df[side].map(wb[prop]).fillna(df[col] if allow_missing else np.nan)
    assert not df[[f"{c}__wb97xd3" for c in DFT]].isna().any().any()
    if missing:
        print(f"WARNING: {len(missing)} ions use r2SCAN-3c values in the wb97xd3 set (test run only)")

    feats = {"counts": STATE + counts,
             "r2scan3c": STATE + DFT + counts,
             "wb97xd3": STATE + [f"{c}__wb97xd3" for c in DFT] + counts}
    y = df[TARGET].to_numpy(float)
    schemes = {"leave_one_anion_out": df.anion.to_numpy()}
    if il_out:
        schemes["leave_one_IL_out"] = df.IL.to_numpy()
    pooled, per_anion = [], []
    t0 = time.time()
    for mname, fn in models().items():
        for sname, groups in schemes.items():
            for fname, cols in feats.items():
                X = df[cols].to_numpy(float)
                oof = np.zeros_like(y)
                for tr, te in LeaveOneGroupOut().split(X, y, groups):
                    oof[te] = fn().fit(X[tr], y[tr]).predict(X[te])
                ae = np.abs(oof - y)
                pooled.append(dict(model=mname, scheme=sname, features=fname, MAE=ae.mean(),
                                   RMSE=np.sqrt(((oof - y) ** 2).mean())))
                if sname == "leave_one_anion_out":
                    for an, idx in df.groupby("anion").groups.items():
                        per_anion.append(dict(model=mname, features=fname, anion=an, n_points=len(idx),
                                              MAE=ae[idx].mean()))
                print(f"[{time.time() - t0:5.0f}s] {mname} {sname:20s} {fname:9s} MAE {ae.mean():.4f}", flush=True)
    pa = pd.DataFrame(per_anion).pivot_table(index=["model", "anion", "n_points"], columns="features",
                                              values="MAE").reset_index()
    pa["wb97xd3_minus_r2scan3c"] = pa.wb97xd3 - pa.r2scan3c
    pa["wb97xd3_minus_counts"] = pa.wb97xd3 - pa.counts
    if not missing:   # a test run with filled-in values writes nothing
        os.makedirs(OUT, exist_ok=True)
        pd.DataFrame(pooled).to_csv(os.path.join(OUT, "orbital_level_comparison.csv"), index=False,
                                    float_format="%.4f")
        pa.to_csv(os.path.join(OUT, "orbital_level_per_anion.csv"), index=False, float_format="%.4f")

    print("\nPooled MAE:")
    print(pd.DataFrame(pooled).pivot_table(index=["model", "scheme"], columns="features", values="MAE")
          [list(feats)].round(4).to_string())
    for m in models():
        s = pa[pa.model == m]
        print(f"\n{m}, leave-one-anion-out, per anion (sorted by wB97X-D3 minus r2SCAN-3c):")
        print(s.sort_values("wb97xd3_minus_r2scan3c")[["anion", "n_points", "counts", "r2scan3c", "wb97xd3",
                                                       "wb97xd3_minus_r2scan3c"]].round(3).to_string(index=False))
        print(f"  wB97X-D3 better than r2SCAN-3c for {(s.wb97xd3 < s.r2scan3c).sum()} of {len(s)} anions; "
              f"better than counts for {(s.wb97xd3 < s.counts).sum()}")
    if missing:
        print(f"\nWARNING: test run, {len(missing)} ions filled with r2SCAN-3c values")


if __name__ == "__main__":
    main()
