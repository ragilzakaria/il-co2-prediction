"""Do COSMO-RS σ-moments help the model predict new anions, beyond ln K_H and the wB97X-D3 descriptors?

Features: 7 σ-moments per ion (area, volume, M2, M3, M4, H-bond acceptor and donor moments;
scripts/cosmors_sigma_moments.py, calculations/cosmors/sigma_moments.csv), for cation and anion = 14
columns, from quantum chemistry only.

Gradient boosting, leave-one-anion-out (23 folds), curated dataset, as scripts/check_cosmors_feature.py:
  counts, counts + σ, wB97X-D3, wB97X-D3 + σ, wB97X-D3 + ln K_H, wB97X-D3 + ln K_H + σ
Null model: the σ-moment vectors shuffled among the ions of the same kind (each cation gets the 7
moments of another cation, each anion those of another anion), as scripts/check_descriptor_noise.py.

Output: paper/data/sigma_moments_feature_null.csv      pooled and per-anion MAE of every null draw
        paper/data/sigma_moments_feature_summary.csv   per anion: MAE of each feature set, null 5-95 %

Run from the repository root (fits spread over --workers processes, one thread each):
    python scripts/check_sigma_moments_feature.py [--draws 20] [--draws-wb 10] [--workers 7]
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_cosmors_feature import OUT, ROOT, arg  # noqa: E402
from check_cosmors_feature import prepare as prepare_kh  # noqa: E402
from check_descriptor_noise import per_anion_mae  # noqa: E402
from check_validation import SEED  # noqa: E402
from cosmors_sigma_moments import MOMENTS  # noqa: E402

SIGMA = os.path.join(ROOT, "calculations", "cosmors", "sigma_moments.csv")
SIDES = ("cation", "anion")


def sig_cols(suffix=""):
    return [f"{side}_sig_{m}{suffix}" for side in SIDES for m in MOMENTS]


def prepare():
    df, y, sets = prepare_kh()
    sm = pd.read_csv(SIGMA).set_index("ion_abbreviation")[MOMENTS]
    for side in SIDES:
        for m in MOMENTS:
            df[f"{side}_sig_{m}"] = df[side].map(sm[m])
    assert df[sig_cols()].notna().all().all(), "missing σ-moments"
    base = {"counts": sets["counts"], "wb97xd3": sets["wb97xd3"], "wb97xd3_KH": sets["wb97xd3_KH"]}
    return df, y, base, sm


_W = {}


def _init():
    _W["df"], _W["y"], _W["base"], _W["sm"] = prepare()


def _fit(task):
    """task: (base set, None) real without σ; (base set, "real") real with σ; (base set, (perm_c, perm_a))
    with the σ-moment table permuted among cations and among anions."""
    name, how = task
    df, y, base, sm = _W["df"], _W["y"], _W["base"], _W["sm"]
    if how is None:
        return per_anion_mae(df, base[name], y)
    if isinstance(how, str):
        return per_anion_mae(df, base[name] + sig_cols(), y)
    d = df.copy()
    for side, perm in zip(SIDES, how):
        ions = sorted(df[side].unique())
        table = sm.loc[ions]
        shuffled = pd.DataFrame(table.to_numpy()[perm], index=table.index, columns=MOMENTS)
        for m in MOMENTS:
            d[f"{side}_sig_{m}__null"] = d[side].map(shuffled[m])
    return per_anion_mae(d, base[name] + sig_cols("__null"), y)


def main():
    draws = {"counts": arg("--draws", 20), "wb97xd3": arg("--draws-wb", 10), "wb97xd3_KH": arg("--draws-wb", 10)}
    workers = arg("--workers", 7)
    df, y, base, _ = prepare()
    n_ion = {side: df[side].nunique() for side in SIDES}

    rng = np.random.default_rng(SEED)
    tasks = {}
    for b in base:
        tasks[("real", b, "")] = (b, None)
        tasks[("real", b, "sig")] = (b, "real")
    for b, n in draws.items():
        for d in range(n):
            tasks[("null", b, d)] = (b, tuple(rng.permutation(n_ion[s]) for s in SIDES))

    t0, res = time.time(), {}
    with ProcessPoolExecutor(workers, initializer=_init) as ex:
        futs = {ex.submit(_fit, t): key for key, t in tasks.items()}
        for f in as_completed(futs):
            key = futs[f]
            res[key] = f.result()
            print(f"[{time.time() - t0:5.0f}s] {len(res):2d}/{len(tasks)} {key[0]} {key[1]} {key[2]}: "
                  f"pooled {res[key][1]:.4f}", flush=True)

    rows = [dict(base=key[1], draw=key[2], pooled=v[1], **v[0].to_dict())
            for key, v in sorted(res.items(), key=lambda kv: str(kv[0])) if key[0] == "null"]
    null = pd.DataFrame(rows)
    os.makedirs(OUT, exist_ok=True)
    null.to_csv(os.path.join(OUT, "sigma_moments_feature_null.csv"), index=False, float_format="%.5f")

    anions = res[("real", "counts", "")][0].index
    s = pd.DataFrame({"n_points": df.groupby("anion").size()}).loc[anions]
    for b in base:
        s[f"MAE_{b}"] = res[("real", b, "")][0]
        s[f"MAE_{b}_sig"] = res[("real", b, "sig")][0]
        s[f"effect_{b}_sig"] = s[f"MAE_{b}_sig"] - s[f"MAE_{b}"]
        eff = null[null.base == b][anions].sub(s[f"MAE_{b}"], axis=1)
        s[f"null_{b}_p05"], s[f"null_{b}_p95"] = eff.quantile(0.05), eff.quantile(0.95)
        s[f"outside_{b}_sig"] = np.where(s[f"effect_{b}_sig"] < s[f"null_{b}_p05"], "better",
                                         np.where(s[f"effect_{b}_sig"] > s[f"null_{b}_p95"], "worse", ""))
    s.reset_index().rename(columns={"index": "anion"}).to_csv(
        os.path.join(OUT, "sigma_moments_feature_summary.csv"), index=False, float_format="%.4f")

    print("\nPooled leave-one-anion-out MAE in x(CO2):")
    for b in base:
        nb = null[null.base == b].pooled
        real = res[("real", b, "sig")][1]
        print(f"  {b:10s} {res[('real', b, '')][1]:.4f}   + σ {real:.4f}   + shuffled σ: mean {nb.mean():.4f}, "
              f"range {nb.min():.4f}-{nb.max():.4f} ({len(nb)} draws); real better than {(nb > real).mean():.0%}")
        print(f"    better than the noise band: {list(s.index[s[f'outside_{b}_sig'] == 'better'])}")
        print(f"    worse than the noise band:  {list(s.index[s[f'outside_{b}_sig'] == 'worse'])}")


if __name__ == "__main__":
    main()
