"""How much does the pooled leave-one-anion-out MAE move with the random seed of gradient boosting?

The feature tests (check_cosmors_feature.py, check_sigma_moments_feature.py) fit every real feature set
once, with SEED. HistGradientBoostingRegressor uses its random_state for the early-stopping validation
split (on by default above 10 000 training rows), so two feature sets that differ by less than the
seed-to-seed spread are not distinguishable. This refits the real sets with several seeds.

Sets: counts, counts + ln K_H, counts + σ, wB97X-D3, wB97X-D3 + ln K_H, wB97X-D3 + σ,
      wB97X-D3 + ln K_H + σ (σ = the 14 σ-moment columns, check_sigma_moments_feature.py).

Output: paper/data/seed_spread.csv   pooled MAE per set and seed

Run from the repository root:
    python scripts/check_seed_spread.py [--seeds 5] [--workers 7]
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import LeaveOneGroupOut

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_cosmors_feature import OUT, arg  # noqa: E402
from check_sigma_moments_feature import prepare, sig_cols  # noqa: E402
from check_validation import SEED  # noqa: E402

_W = {}


def feature_sets(base):
    return {"counts": base["counts"], "counts_KH": base["counts"] + ["lnKH"],
            "counts_sig": base["counts"] + sig_cols(), "wb97xd3": base["wb97xd3"],
            "wb97xd3_KH": base["wb97xd3_KH"], "wb97xd3_sig": base["wb97xd3"] + sig_cols(),
            "wb97xd3_KH_sig": base["wb97xd3_KH"] + sig_cols()}


def _init():
    df, y, base, _ = prepare()
    _W.update(df=df, y=y, sets=feature_sets(base))


def _fit(task):
    name, seed = task
    df, y = _W["df"], _W["y"]
    X = df[_W["sets"][name]].to_numpy(float)
    oof = np.zeros_like(y)
    for tr, te in LeaveOneGroupOut().split(X, y, df.anion):
        m = HistGradientBoostingRegressor(max_iter=500, learning_rate=0.05, random_state=seed)
        oof[te] = m.fit(X[tr], y[tr]).predict(X[te])
    return np.abs(oof - y).mean()


def main():
    n_seeds, workers = arg("--seeds", 5), arg("--workers", 7)
    seeds = [SEED] + list(range(1, n_seeds))
    names = list(feature_sets({"counts": [], "wb97xd3": [], "wb97xd3_KH": []}))
    tasks = [(n, s) for s in seeds for n in names]
    t0, rows = time.time(), []
    with ProcessPoolExecutor(workers, initializer=_init) as ex:
        futs = {ex.submit(_fit, t): t for t in tasks}
        for f in as_completed(futs):
            n, s = futs[f]
            rows.append(dict(set=n, seed=s, pooled=f.result()))
            print(f"[{time.time() - t0:5.0f}s] {len(rows):2d}/{len(tasks)} {n} seed {s}: {rows[-1]['pooled']:.4f}",
                  flush=True)
    r = pd.DataFrame(rows).sort_values(["set", "seed"])
    os.makedirs(OUT, exist_ok=True)
    r.to_csv(os.path.join(OUT, "seed_spread.csv"), index=False, float_format="%.5f")
    t = r.groupby("set").pooled.agg(["mean", "std", "min", "max"]).loc[names]
    print("\nPooled leave-one-anion-out MAE in x(CO2) over seeds", seeds)
    print(t.round(4).to_string())


if __name__ == "__main__":
    main()
