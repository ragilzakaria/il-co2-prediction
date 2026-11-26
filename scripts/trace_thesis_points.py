"""Trace each point of the thesis training set to a publication, where possible.

The thesis set (data/thesis/final_model_dataset_v3.csv, 11,256 points) comes from Tables S5 and S6
of the Supporting Information of Laakso et al. 2025 (doi:10.1016/j.ces.2025.121226), which took the
experimental data from Song et al. 2020 (doi:10.1016/j.ces.2020.115752); those tables carry no
per-point sources. This script matches every thesis point to a point of the referenced ILThermo
dataset (same IL; T within 0.25 K; P within 1.5 % or 1.5 kPa; x_CO2 within 2 % or 0.003) and records
the citation and verified DOI of the match. It also flags exact duplicate rows.

Output: data/thesis/thesis_point_references.csv (one row per thesis point, in file order)

Run from the repository root:
    python scripts/trace_thesis_points.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
THESIS = ROOT / "data" / "thesis" / "final_model_dataset_v3.csv"
REF = ROOT / "data" / "processed" / "referenced_co2_solubility_dataset.csv"
OUT = ROOT / "data" / "thesis" / "thesis_point_references.csv"
COLS = ["cation", "anion", "Temperature_K", "Pressure_kPa", "CO2_solubility"]


def main():
    th = pd.read_csv(THESIS, usecols=COLS)
    ref = pd.read_csv(REF)
    groups = {k: g for k, g in ref.groupby(["cation", "anion"])}
    cits, dois, n_cand = [], [], []
    for r in th.itertuples(index=False):
        g = groups.get((r.cation, r.anion))
        m = None if g is None else g[
            (abs(g.Temperature_K - r.Temperature_K) <= 0.25)
            & (abs(g.Pressure_kPa - r.Pressure_kPa) <= max(0.015 * r.Pressure_kPa, 1.5))
            & (abs(g.CO2_solubility - r.CO2_solubility) <= max(0.003, 0.02 * r.CO2_solubility))]
        if m is None or m.empty:
            cits.append(""); dois.append(""); n_cand.append(0)
        else:
            # closest match in solubility
            best = m.iloc[(m.CO2_solubility - r.CO2_solubility).abs().argmin()]
            cits.append(best.ref_citation); dois.append(best.ref_doi); n_cand.append(m.ref_doi.nunique())
    out = th.assign(
        dataset_source="Song et al. 2020, doi:10.1016/j.ces.2020.115752 (via Laakso et al. 2025 SI, doi:10.1016/j.ces.2025.121226)",
        matched_citation=cits, matched_doi=dois, n_candidate_publications=n_cand,
        exact_duplicate_of_earlier_row=th.duplicated(COLS))
    out.to_csv(OUT, index=False)
    matched = out.matched_doi != ""
    il = out.groupby(["cation", "anion"]).apply(lambda g: (g.matched_doi != "").mean(), include_groups=False)
    print(f"{len(out)} thesis points; matched to a referenced publication: {matched.sum()} ({matched.mean():.1%})")
    print(f"ILs: {len(il)}; fully matched {int((il == 1).sum())}, partly {int(((il > 0) & (il < 1)).sum())}, "
          f"not at all {int((il == 0).sum())}")
    print(f"exact duplicate rows: {out.exact_duplicate_of_earlier_row.sum()}")
    print(f"distinct publications matched: {out.loc[matched, 'matched_doi'].nunique()}")
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
