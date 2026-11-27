"""Are the per-anion effects of the DFT descriptors real, or noise? A shuffled-descriptor null test.

scripts/per_il_errors.py found that adding the 10 r2SCAN-3c descriptors changes the
leave-one-anion-out error of single anions by up to +-0.03 in x(CO2) (worse for the small oxyanions,
better for [BF4], [Tf2N] and others). Tree models react to any extra ion-level column, so part of
that could be noise. Null model: the same 10 descriptor columns, but the descriptor vectors are
shuffled among the ions (each anion gets the five values of another anion, each cation those of
another cation). The values keep their real distribution, but carry no chemistry. Repeated N times.

For every held-out anion this gives the spread of "effect of adding 10 meaningless ion-level
numbers". A real effect is one outside that spread.

Model: histogram gradient boosting (settings of scripts/compare_descriptor_sets.py), validation:
leave-one-anion-out (23 folds), curated dataset. If data/processed/dft_descriptors_wb97xd3.csv
holds every ion, the wB97X-D3/def2-TZVPD descriptors are tested against the same null.

Output: paper/data/descriptor_noise_null.csv      MAE per anion for every null draw
        paper/data/descriptor_noise_summary.csv   per anion: real effect, null mean, 5-95 % range,
                                                  share of null draws at least as extreme
        images/checks/descriptor_noise_null.png

Run from the repository root (about 45 min with the default 20 draws):
    python scripts/check_descriptor_noise.py [--draws N]
    python scripts/check_descriptor_noise.py --plot-only    # redraw the figure from the summary CSV
    python scripts/check_descriptor_noise.py --add-wb97xd3  # place the wB97X-D3 result of
        scripts/compare_orbital_levels.py (same model and folds) against the saved null draws
"""
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import LeaveOneGroupOut

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_validation import DFT, SEED, STATE, TARGET, load  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
DATA = os.path.join(ROOT, "data", "processed", "model_dataset_curated.csv")
WB97X = os.path.join(ROOT, "data", "processed", "dft_descriptors_wb97xd3.csv")
OUT = os.path.join(ROOT, "paper", "data")
FIG = os.path.join(ROOT, "images", "checks", "descriptor_noise_null.png")
PROPS = ["Energy_Hartree", "HOMO_Eh", "LUMO_Eh", "HOMO_LUMO_Gap_Eh", "Dipole_Moment_Debye"]


def model():
    return HistGradientBoostingRegressor(max_iter=500, learning_rate=0.05, random_state=SEED)


def per_anion_mae(df, cols, y):
    X = df[cols].to_numpy(float)
    oof = np.zeros_like(y)
    for tr, te in LeaveOneGroupOut().split(X, y, df.anion):
        oof[te] = model().fit(X[tr], y[tr]).predict(X[te])
    ae = pd.Series(np.abs(oof - y), index=df.index)
    return ae.groupby(df.anion).mean(), ae.mean()


def shuffled(df, rng):
    """Copy of the DFT columns with the per-ion descriptor vectors permuted among ions of a kind."""
    out = {}
    for side in ("cation", "anion"):
        cols = [f"{side}_{p}" for p in PROPS]
        table = df.groupby(side)[cols].first()
        perm = table.copy()
        perm[:] = table.to_numpy()[rng.permutation(len(table))]
        mapped = df[[side]].join(perm, on=side)
        for c in cols:
            out[c + "__null"] = mapped[c].to_numpy()
    return pd.DataFrame(out, index=df.index)


def add_wb97xd3():
    """Add the wB97X-D3 per-anion effect from paper/data/orbital_level_per_anion.csv (GB rows) to the
    summary, using the saved null draws; no refitting of the null models."""
    s = pd.read_csv(os.path.join(OUT, "descriptor_noise_summary.csv")).set_index("anion")
    null = pd.read_csv(os.path.join(OUT, "descriptor_noise_null.csv"))
    lev = pd.read_csv(os.path.join(OUT, "orbital_level_per_anion.csv"))
    lev = lev[lev.model == "GB"].set_index("anion")
    # same deterministic model and folds: the counts and r2SCAN-3c errors must agree
    assert (lev.counts - s.MAE_counts).abs().max() < 1e-3, "counts baseline differs"
    assert (lev.r2scan3c - s.MAE_counts - s.effect_r2scan3c).abs().max() < 1e-3, "r2SCAN-3c result differs"
    eff_null = null.pivot(index="draw", columns="anion", values="MAE") - s.MAE_counts
    s["effect_wb97xd3"] = lev.wb97xd3 - s.MAE_counts
    for a in s.index:
        e, eff = eff_null[a], s.at[a, "effect_wb97xd3"]
        s.at[a, "p_wb97xd3"] = (np.abs(e - e.mean()) >= abs(eff - e.mean())).mean()
        s.at[a, "outside_null_wb97xd3"] = not (e.quantile(0.05) <= eff <= e.quantile(0.95))
    s = s.reset_index()
    s.to_csv(os.path.join(OUT, "descriptor_noise_summary.csv"), index=False, float_format="%.4f")
    plot(s, ["r2scan3c", "wb97xd3"], null.draw.nunique())
    cols = ["anion", "n_points", "null_mean_effect", "null_p05", "null_p95", "effect_r2scan3c",
            "outside_null_r2scan3c", "effect_wb97xd3", "outside_null_wb97xd3"]
    print(s[cols].round(3).to_string(index=False))
    for l in ("r2scan3c", "wb97xd3"):
        below = (s[f"effect_{l}"] < s.null_p05).sum()
        above = (s[f"effect_{l}"] > s.null_p95).sum()
        print(f"{l}: better than the noise band for {below} anions, worse for {above}, of {len(s)}")


def main():
    if "--add-wb97xd3" in sys.argv:
        add_wb97xd3()
        return
    if "--plot-only" in sys.argv:
        s = pd.read_csv(os.path.join(OUT, "descriptor_noise_summary.csv"))
        n = pd.read_csv(os.path.join(OUT, "descriptor_noise_null.csv")).draw.nunique()
        plot(s, [c[len("effect_"):] for c in s.columns if c.startswith("effect_")], n)
        return
    n_draws = int(sys.argv[sys.argv.index("--draws") + 1]) if "--draws" in sys.argv else 20
    df, _, counts = load(DATA)
    y = df[TARGET].to_numpy(float)
    t0 = time.time()

    real = {}
    base, base_pooled = per_anion_mae(df, STATE + counts, y)
    real["r2scan3c"] = per_anion_mae(df, STATE + DFT + counts, y)
    print(f"[{time.time() - t0:5.0f}s] counts and r2SCAN-3c done", flush=True)
    if os.path.exists(WB97X):
        wb = pd.read_csv(WB97X).set_index("ion_abbreviation")
        if set(df.cation) | set(df.anion) <= set(wb.index):
            for side in ("cation", "anion"):
                for p in PROPS:
                    df[f"{side}_{p}__wb97xd3"] = df[side].map(wb[p])
            real["wb97xd3"] = per_anion_mae(df, STATE + [f"{c}__wb97xd3" for c in DFT] + counts, y)
            print(f"[{time.time() - t0:5.0f}s] wB97X-D3 done", flush=True)
        else:
            print("wB97X-D3 descriptors incomplete: tested only r2SCAN-3c")

    rng = np.random.default_rng(SEED)
    null_rows, null_pooled = [], []
    for k in range(n_draws):
        nd = df.join(shuffled(df, rng))
        mae, pooled = per_anion_mae(nd, STATE + [c + "__null" for c in DFT] + counts, y)
        null_pooled.append(pooled)
        null_rows += [dict(draw=k, anion=a, MAE=v) for a, v in mae.items()]
        print(f"[{time.time() - t0:5.0f}s] null draw {k + 1}/{n_draws}: pooled MAE {pooled:.4f}", flush=True)
    null = pd.DataFrame(null_rows)
    os.makedirs(OUT, exist_ok=True)
    null.to_csv(os.path.join(OUT, "descriptor_noise_null.csv"), index=False, float_format="%.5f")

    eff_null = null.pivot(index="draw", columns="anion", values="MAE") - base   # null effect vs counts
    n_pts = df.anion.value_counts()
    rows = []
    for a in base.index:
        e = eff_null[a]
        r = dict(anion=a, n_points=n_pts[a], MAE_counts=base[a], null_mean_effect=e.mean(),
                 null_p05=e.quantile(0.05), null_p95=e.quantile(0.95))
        for lev, (mae, _) in real.items():
            eff = mae[a] - base[a]
            r[f"effect_{lev}"] = eff
            # share of null draws whose effect is at least as far from the null mean
            r[f"p_{lev}"] = (np.abs(e - e.mean()) >= abs(eff - e.mean())).mean()
            r[f"outside_null_{lev}"] = not (e.quantile(0.05) <= eff <= e.quantile(0.95))
        rows.append(r)
    s = pd.DataFrame(rows).sort_values("effect_r2scan3c")
    s.to_csv(os.path.join(OUT, "descriptor_noise_summary.csv"), index=False, float_format="%.4f")
    plot(s, list(real), n_draws)

    print(f"\npooled MAE: counts {base_pooled:.4f}; " +
          "; ".join(f"{lev} {p:.4f}" for lev, (_, p) in real.items()) +
          f"; shuffled {np.mean(null_pooled):.4f} (5-95 %: {np.quantile(null_pooled, 0.05):.4f}-"
          f"{np.quantile(null_pooled, 0.95):.4f})")
    cols = ["anion", "n_points", "MAE_counts", "null_mean_effect", "null_p05", "null_p95"] + \
           [c for lev in real for c in (f"effect_{lev}", f"p_{lev}", f"outside_null_{lev}")]
    print(s[cols].round(3).to_string(index=False))
    for lev in real:
        print(f"{lev}: outside the null 5-95 % range for {s[f'outside_null_{lev}'].sum()} of {len(s)} anions")
    print(f"Total {time.time() - t0:.0f}s")


def plot(s, levels, n_draws):
    """Per anion: the null 5-95 % range, its mean, and the real effect of each descriptor level."""
    fig, ax = plt.subplots(figsize=(8, 7))
    yy = np.arange(len(s))
    ax.hlines(yy, s.null_p05, s.null_p95, color="#b5b3ad", lw=6, label="shuffled descriptors, 5–95 %")
    ax.plot(s.null_mean_effect, yy, "|", color="#52514e", ms=10, label="shuffled descriptors, mean")
    marks = {"r2scan3c": ("o", "#2a78d6", "r²SCAN-3c descriptors"),
             "wb97xd3": ("s", "#eb6834", "ωB97X-D3/def2-TZVPD descriptors")}
    for lev in levels:
        m, c, lab = marks[lev]
        ax.plot(s[f"effect_{lev}"], yy, m, color=c, ms=7, label=lab, markeredgecolor="#fcfcfb", mew=1.5)
    ax.axvline(0, color="#52514e", lw=0.8)
    ax.set_yticks(yy)
    ax.set_yticklabels([f"{a} (n={n})" for a, n in zip(s.anion, s.n_points)], fontsize=8)
    ax.set_xlabel("change in leave-one-anion-out MAE from adding 10 ion descriptors")
    ax.set_title(f"Real descriptors against {n_draws} shuffled sets (GB)", fontsize=10, loc="left")
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#e6e5e1", lw=0.6)
    ax.set_axisbelow(True)
    fig.tight_layout()
    os.makedirs(os.path.dirname(FIG), exist_ok=True)
    fig.savefig(FIG, dpi=160)
    print(f"wrote {os.path.relpath(FIG, ROOT)}")


if __name__ == "__main__":
    main()
