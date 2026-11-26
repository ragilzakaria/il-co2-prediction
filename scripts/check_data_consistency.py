"""Physical and cross-source consistency checks of the curated CO2-solubility dataset.

1. Monotonicity: within one data set and one isotherm (T within 0.5 K), x_CO2 must not fall as
   the pressure rises. A drop counts as a violation when it exceeds max(2 % of x, 0.002, the
   stated uncertainty).
2. Agreement between publications: each point is compared with the isotherms of other
   publications for the same IL at the same temperature (within 1 K). The other isotherm is
   interpolated linearly in ln x against ln P, only inside its own pressure range (no
   extrapolation). The consensus value is the median over the other publications. Each data set
   gets the median relative deviation of its points from the consensus.

Output: data/processed/checks/consistency_by_set.csv, consistency_by_point.csv

Run from the repository root:
    python scripts/check_data_consistency.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "processed" / "curated_co2_solubility_dataset.csv"
OUT = ROOT / "data" / "processed" / "checks"


def isotherms(g):
    """Split one data set into isotherms: points within 0.5 K of each other, >= 2 points."""
    g = g.sort_values("Temperature_K")
    labels = (g.Temperature_K.diff().fillna(0) > 0.5).cumsum()
    return [iso.sort_values("Pressure_kPa") for _, iso in g.groupby(labels) if len(iso) >= 2]


def interp(iso, P):
    iso = iso.groupby("Pressure_kPa").CO2_solubility.mean()
    p, x = iso.index.values, iso.values
    if len(p) < 2 or not (p.min() <= P <= p.max()):
        return None
    return float(np.exp(np.interp(np.log(P), np.log(p), np.log(x))))


def main():
    df = pd.read_csv(DATA)
    OUT.mkdir(parents=True, exist_ok=True)
    df["mono_violation"] = False
    for sid, g in df.groupby("set_id"):
        for iso in isotherms(g):
            x = iso.CO2_solubility.values
            u = iso.solubility_uncertainty.fillna(0).values
            for k in range(1, len(x)):
                tol = max(0.02 * x[k - 1], 0.002, u[k], u[k - 1])
                if iso.Pressure_kPa.values[k] > iso.Pressure_kPa.values[k - 1] and x[k] < x[k - 1] - tol:
                    df.loc[iso.index[k], "mono_violation"] = True

    # isotherms of every publication, per IL
    by_il = {}
    for (cat, an, cit), g in df.groupby(["cation", "anion", "ref_citation"]):
        by_il.setdefault((cat, an), []).append((cit, isotherms(g)))
    dev, n_other = np.full(len(df), np.nan), np.zeros(len(df), int)
    spread, absdev = np.full(len(df), np.nan), np.full(len(df), np.nan)
    for pos, r in enumerate(df.itertuples()):
        refs = []
        for cit, isos in by_il[(r.cation, r.anion)]:
            if cit == r.ref_citation:
                continue
            vals = [interp(iso, r.Pressure_kPa) for iso in isos
                    if abs(iso.Temperature_K.mean() - r.Temperature_K) <= 1.0]
            vals = [v for v in vals if v]
            if vals:
                refs.append(np.median(vals))
        if refs:
            dev[pos] = r.CO2_solubility / np.median(refs) - 1
            absdev[pos] = r.CO2_solubility - np.median(refs)
            n_other[pos] = len(refs)
            spread[pos] = max(refs) / min(refs) - 1 if len(refs) >= 2 else np.nan
    df["rel_dev_from_others"], df["n_other_publications"] = dev, n_other
    df["abs_dev_from_others"], df["spread_of_others"] = absdev, spread
    # a point disagrees when it is > 30 % and > 0.01 in mole fraction from the consensus of >= 2
    # other publications that agree with each other within 15 %
    df["disagrees"] = ((df.n_other_publications >= 2) & (df.spread_of_others <= 0.15)
                       & (df.rel_dev_from_others.abs() > 0.30) & (df.abs_dev_from_others.abs() > 0.01))
    df.to_csv(OUT / "consistency_by_point.csv", index=False)

    rows = []
    for sid, g in df.groupby("set_id"):
        c = g[g.n_other_publications >= 2]
        rows.append({"set_id": sid, "cation": g.cation.iloc[0], "anion": g.anion.iloc[0],
                     "ref_citation": g.ref_citation.iloc[0], "n_points": len(g),
                     "mono_violations": int(g.mono_violation.sum()),
                     "n_compared": len(c), "median_rel_dev": c.rel_dev_from_others.median() if len(c) else np.nan,
                     "median_spread_of_others": c.spread_of_others.median() if len(c) else np.nan,
                     "n_disagree": int(g.disagrees.sum()),
                     "x_range": f"{g.CO2_solubility.min():.3g}-{g.CO2_solubility.max():.3g}",
                     "T_range": f"{g.Temperature_K.min():.0f}-{g.Temperature_K.max():.0f}",
                     "P_range_kPa": f"{g.Pressure_kPa.min():.3g}-{g.Pressure_kPa.max():.3g}"})
    s = pd.DataFrame(rows)
    # verdict per data set
    s["verdict"] = "not comparable"
    comp = s.n_compared >= 3
    s.loc[comp, "verdict"] = "agrees"
    s.loc[comp & (s.median_rel_dev.abs() > 0.30) & (s.median_spread_of_others > 0.15), "verdict"] = "conflicting sources"
    s.loc[comp & (s.median_rel_dev.abs() > 0.30) & (s.median_spread_of_others <= 0.15), "verdict"] = "outlier"
    s.to_csv(OUT / "consistency_by_set.csv", index=False)
    print(f"points compared with >= 2 other publications: {(df.n_other_publications >= 2).sum()} of {len(df)}")
    d = df.loc[df.n_other_publications >= 2, "rel_dev_from_others"].abs()
    print(f"  |relative deviation| median {d.median():.3f}, 90th pct {d.quantile(0.9):.3f}, >30 %: {(d > 0.3).sum()}")
    print(f"monotonicity violations: {df.mono_violation.sum()} points in {(s.mono_violations > 0).sum()} sets")
    print(f"sets compared (>= 3 points vs >= 2 publications): {(s.n_compared >= 3).sum()} of {len(s)}")
    print("verdicts:", s.verdict.value_counts().to_dict())
    print(f"points that disagree with a coherent consensus: {df.disagrees.sum()}")


if __name__ == "__main__":
    main()
