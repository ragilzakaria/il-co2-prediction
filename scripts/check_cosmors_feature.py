"""Does the COSMO-RS CO2 Henry constant help the model predict new anions?

Feature: ln K_H at 25 C from openCOSMO-RS (scripts/cosmors_henry.py, calculations/cosmors/henry.csv),
one value per ionic liquid, computed from quantum chemistry only (no experimental input). Its
temperature dependence is not used (wrong sign in openCOSMO-RS); the model has T as an input.

Test, gradient boosting with the settings of scripts/check_descriptor_noise.py, leave-one-anion-out
(23 folds), curated dataset:
  counts            T, P and atom counts (baseline)
  counts + K_H      the same plus ln K_H
  wB97X-D3          T, P, the 10 wB97X-D3 ion descriptors and atom counts (current best)
  wB97X-D3 + K_H    the same plus ln K_H
Null model: ln K_H shuffled among the 110 ionic liquids (each IL gets the value of another IL), so
the column keeps its distribution but carries no chemistry. A real gain is one outside the spread of
the shuffled draws.

Output: paper/data/cosmors_feature_null.csv      pooled and per-anion MAE of every null draw
        paper/data/cosmors_feature_summary.csv   per anion: MAE of each feature set, null 5-95 %

Run from the repository root (34 fits of 23 folds, spread over --workers processes with one thread
each, because the OpenMP threads of gradient boosting oversubscribe the machine otherwise):
    python scripts/check_cosmors_feature.py [--draws 20] [--draws-wb 10] [--workers 7]
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_validation import DFT, SEED, STATE, TARGET, load  # noqa: E402
from check_descriptor_noise import PROPS, per_anion_mae  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
DATA = os.path.join(ROOT, "data", "processed", "model_dataset_curated.csv")
WB97X = os.path.join(ROOT, "data", "processed", "dft_descriptors_wb97xd3.csv")
HENRY = os.path.join(ROOT, "calculations", "cosmors", "henry.csv")
MANIFEST = os.path.join(ROOT, "calculations", "conformers_crest", "manifest.csv")
OUT = os.path.join(ROOT, "paper", "data")


def arg(name, default):
    return int(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default


def prepare():
    """Curated dataset with the wB97X-D3 columns and ln K_H, and the four feature sets."""
    df, _, counts = load(DATA)
    y = df[TARGET].to_numpy(float)

    wb = pd.read_csv(WB97X).set_index("ion_abbreviation")
    for side in ("cation", "anion"):
        for p in PROPS:
            df[f"{side}_{p}__wb"] = df[side].map(wb[p])
    wbcols = [f"{c}__wb" for c in DFT]

    man = pd.read_csv(MANIFEST).set_index("ion_abbreviation").file
    h = pd.read_csv(HENRY)
    h = h[h.T_K.round(2) == 298.15].set_index("system").KH_MPa
    df["il"] = [man[c] + "__" + man[a] for c, a in zip(df.cation, df.anion)]
    df["lnKH"] = np.log(df.il.map(h))
    assert df.lnKH.notna().all(), "missing K_H for " + str(df.il[df.lnKH.isna()].unique())
    sets = {"counts": STATE + counts, "counts_KH": STATE + counts + ["lnKH"],
            "wb97xd3": STATE + wbcols + counts, "wb97xd3_KH": STATE + wbcols + counts + ["lnKH"]}
    return df, y, sets


_W = {}


def _init():
    _W["df"], _W["y"], _W["sets"] = prepare()
    _W["table"] = _W["df"].groupby("il").lnKH.first()


def _fit(task):
    """task: (set name, None) for a real set, or (base set, permutation of the IL table) for a null."""
    name, perm = task
    df, y, sets = _W["df"], _W["y"], _W["sets"]
    if perm is None:
        return per_anion_mae(df, sets[name], y)
    t = _W["table"]
    d = df.assign(lnKH_null=df.il.map(pd.Series(t.to_numpy()[perm], index=t.index)))
    return per_anion_mae(d, sets[name] + ["lnKH_null"], y)


def main():
    n_draws, n_draws_wb, workers = arg("--draws", 20), arg("--draws-wb", 10), arg("--workers", 7)
    df, y, sets = prepare()
    n_il = df.il.nunique()

    # null: ln K_H permuted among the ionic liquids
    rng = np.random.default_rng(SEED)
    tasks = {("real", k, 0): (k, None) for k in sets}
    for base, n in (("counts", n_draws), ("wb97xd3", n_draws_wb)):
        for d in range(n):
            tasks[("null", base, d)] = (base, rng.permutation(n_il))

    t0, res = time.time(), {}
    with ProcessPoolExecutor(workers, initializer=_init) as ex:
        futs = {ex.submit(_fit, t): key for key, t in tasks.items()}
        for f in as_completed(futs):
            key = futs[f]
            res[key] = f.result()
            print(f"[{time.time() - t0:5.0f}s] {len(res):2d}/{len(tasks)} {key[0]} {key[1]} {key[2]}: "
                  f"pooled {res[key][1]:.4f}", flush=True)

    real = {k: res[("real", k, 0)] for k in sets}
    rows = [dict(base=key[1], draw=key[2], pooled=v[1], **v[0].to_dict())
            for key, v in sorted(res.items()) if key[0] == "null"]
    null = pd.DataFrame(rows)
    os.makedirs(OUT, exist_ok=True)
    null.to_csv(os.path.join(OUT, "cosmors_feature_null.csv"), index=False, float_format="%.5f")

    anions = real["counts"][0].index
    s = pd.DataFrame({"n_points": df.groupby("anion").size()}).loc[anions]
    for k in sets:
        s[f"MAE_{k}"] = real[k][0]
    for base, k in (("counts", "counts_KH"), ("wb97xd3", "wb97xd3_KH")):
        nb = null[null.base == base]
        s[f"effect_{k}"] = s[f"MAE_{k}"] - s[f"MAE_{base}"]
        eff = nb[anions].sub(s[f"MAE_{base}"], axis=1)
        s[f"null_{k}_p05"], s[f"null_{k}_p95"] = eff.quantile(0.05), eff.quantile(0.95)
        s[f"outside_{k}"] = np.where(s[f"effect_{k}"] < s[f"null_{k}_p05"], "better",
                                     np.where(s[f"effect_{k}"] > s[f"null_{k}_p95"], "worse", ""))
    s.reset_index().rename(columns={"index": "anion"}).to_csv(
        os.path.join(OUT, "cosmors_feature_summary.csv"), index=False, float_format="%.4f")

    print("\nPooled leave-one-anion-out MAE in x(CO2):")
    for base, k in (("counts", "counts_KH"), ("wb97xd3", "wb97xd3_KH")):
        nb = null[null.base == base].pooled
        print(f"  {base:8s} {real[base][1]:.4f}   + ln K_H {real[k][1]:.4f}   + shuffled ln K_H: mean "
              f"{nb.mean():.4f}, range {nb.min():.4f}-{nb.max():.4f} ({len(nb)} draws); real better than "
              f"{(nb > real[k][1]).mean():.0%} of them")
        print(f"    better than the noise band: {list(s.index[s[f'outside_{k}'] == 'better'])}")
        print(f"    worse than the noise band:  {list(s.index[s[f'outside_{k}'] == 'worse'])}")


if __name__ == "__main__":
    main()
