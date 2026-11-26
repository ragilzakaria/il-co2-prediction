"""Which ionic liquids and ion families does the model predict worst? Curated dataset.

Two group-wise validations on data/processed/model_dataset_curated.csv (12,026 points, 110 ILs):

  leave-one-IL-out      110 folds: each IL is predicted by a model trained on the other 109
                        (its cation and anion may occur in other ILs of the training set)
  leave-one-anion-out    23 folds: each anion with all its ILs is held out (a new anion)

Models: Random Forest and histogram gradient boosting (settings as in
scripts/compare_descriptor_sets.py), each with two feature sets:
  counts   T, P + elemental atom counts of both ions (corrected structures)
  dft      T, P + the 10 corrected r2SCAN-3c descriptors + atom counts

For every IL the output also lists what may explain a large error: number of publications and data
sets, the share of points that disagree with other publications (scripts/check_data_consistency.py),
and how many other ILs share its cation and its anion.

Output
  paper/data/per_il_errors_curated.csv        one row per IL and model/feature set
  paper/data/per_family_errors_curated.csv    one row per anion family / cation family and model
  paper/data/per_anion_errors_curated.csv     leave-one-anion-out error per anion, counts against dft
  images/checks/per_il_errors.png             worst ILs, anion-family errors, effect of the DFT
                                              descriptors per held-out anion

Run from the repository root (about 40 min on a laptop; --plot-only redraws the figure from the CSVs):
    python scripts/per_il_errors.py [--plot-only]
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
from sklearn.model_selection import LeaveOneGroupOut

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_validation import DFT, SEED, STATE, TARGET, load  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
DATA = os.path.join(ROOT, "data", "processed", "model_dataset_curated.csv")
CONS = os.path.join(ROOT, "data", "processed", "checks", "consistency_by_point.csv")
OUT = os.path.join(ROOT, "paper", "data")
FIG = os.path.join(ROOT, "images", "checks", "per_il_errors.png")

ANION_FAMILY = {
    "sulfonyl amide/methide": ["[Tf2N]", "[BETA]", "[FSA]", "[methide]"],
    "fluoro-borate/phosphate": ["[BF4]", "[PF6]"],
    "cyano/thiocyanate": ["[C(CN)3]", "[DCA]", "[SCN]", "[TCB]"],
    "sulfonate": ["[TfO]", "[MeSO3]", "[NfO]", "[TOS]"],
    "sulfate": ["[MeSO4]", "[EtSO4]", "[MDEGSO4]", "[HSO4]"],
    "phosphate": ["[DEPO4]", "[DMPO4]"],
    "carboxylate": ["[TFA]", "[C3F7CO2]"],
    "nitrate": ["[NO3]"],
}
CATION_FAMILY = {
    "imidazolium, alkyl": ["[MMIM]", "[EMIM]", "[PMIM]", "[BMIM]", "[C5MIM]", "[HMIM]", "[C7MIM]", "[OMIM]",
                           "[NMIM]", "[DMIM]", "[BMMIM]", "[BBIM]"],
    "imidazolium, functionalised": ["[HOC2MIM]", "[HOC3MIM]", "[C8F13MIM]", "[aMIM]"],
    "pyridinium": ["[C4PY]", "[C4MPY]", "[PMPY]", "[HMPY]", "[HOEPY]"],
    "pyrrolidinium/piperidinium": ["[C3MPYR]", "[BMPYR]", "[C5MPYR]", "[HMPYR]", "[C7MPYR]", "[OMPYR]",
                                   "[COCMPYR]", "[PMPip]"],
    "ammonium": ["[N4,1,1,1]", "[N6,1,1,1]", "[N2,1,1,3]", "[N1,4,4,4]", "[deme]", "[thma]", "[N1,1,1,2-OH]"],
    "phosphonium/sulfonium": ["[P1,4,4,4]", "[P2,4,4,4]", "[S222]"],
}
# categorical palette, fixed order (dataviz reference palette, light mode)
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]


def family_of(mapping):
    return {ion: fam for fam, ions in mapping.items() for ion in ions}


def models():
    return {
        "RF": lambda: RandomForestRegressor(n_estimators=300, random_state=SEED, n_jobs=-1),
        "GB": lambda: HistGradientBoostingRegressor(max_iter=500, learning_rate=0.05, random_state=SEED),
    }


def oof_predictions(fn, X, y, groups):
    oof = np.zeros_like(y)
    for tr, te in LeaveOneGroupOut().split(X, y, groups):
        oof[te] = fn().fit(X[tr], y[tr]).predict(X[te])
    return oof


def main():
    if "--plot-only" in sys.argv:
        plot(pd.read_csv(os.path.join(OUT, "per_il_errors_curated.csv")),
             pd.read_csv(os.path.join(OUT, "per_family_errors_curated.csv")))
        return
    df, _, counts = load(DATA)
    an_fam, cat_fam = family_of(ANION_FAMILY), family_of(CATION_FAMILY)
    missing = (set(df.anion) - set(an_fam)) | (set(df.cation) - set(cat_fam))
    assert not missing, f"ions without a family: {missing}"
    df["anion_family"], df["cation_family"] = df.anion.map(an_fam), df.cation.map(cat_fam)

    cons = pd.read_csv(CONS, usecols=["set_id", "Temperature_K", "Pressure_kPa", TARGET, "disagrees",
                                      "n_other_publications"])
    key = ["set_id", "Temperature_K", "Pressure_kPa", TARGET]
    assert not cons.duplicated(key).any()
    df = df.merge(cons, on=key, how="left", validate="one_to_one")
    assert df.disagrees.notna().all(), "every point must have a consistency check"

    y = df[TARGET].to_numpy(float)
    feats = {"counts": STATE + counts, "dft": STATE + DFT + counts}
    schemes = {"leave_one_IL_out": df.IL.to_numpy(), "leave_one_anion_out": df.anion.to_numpy()}
    t0 = time.time()
    for mname, fn in models().items():
        for fname, cols in feats.items():
            X = df[cols].to_numpy(float)
            for sname, groups in schemes.items():
                col = f"pred_{mname}_{fname}_{sname}"
                df[col] = oof_predictions(fn, X, y, groups)
                mae = np.abs(df[col] - y).mean()
                print(f"[{time.time() - t0:5.0f}s] {mname} {fname:6s} {sname:20s} pooled MAE {mae:.4f}", flush=True)

    # ---- per IL
    ils = df.groupby("IL")
    base = pd.DataFrame({
        "cation": ils.cation.first(), "anion": ils.anion.first(),
        "cation_family": ils.cation_family.first(), "anion_family": ils.anion_family.first(),
        "n_points": ils.size(), "n_publications": ils.ref_citation.nunique(), "n_sets": ils.set_id.nunique(),
        "x_mean": ils[TARGET].mean(), "x_max": ils[TARGET].max(),
        "T_range": ils.Temperature_K.agg(lambda s: f"{s.min():.0f}-{s.max():.0f}"),
        "P_max_kPa": ils.Pressure_kPa.max(),
        "share_disagreeing_points": ils.disagrees.mean(),
        "share_points_comparable": ils.n_other_publications.agg(lambda s: (s >= 2).mean()),
    })
    il_of = df.drop_duplicates("IL")
    base["n_other_ILs_same_cation"] = base.cation.map(il_of.cation.value_counts()) - 1
    base["n_other_ILs_same_anion"] = base.anion.map(il_of.anion.value_counts()) - 1
    rows = []
    for col in [c for c in df.columns if c.startswith("pred_")]:
        _, mname, fname, sname = col.split("_", 3)
        e = df[col] - y
        g = pd.DataFrame({"IL": df.IL, "e": e, "ae": e.abs()}).groupby("IL")
        r = base.copy()
        r.insert(0, "model", mname); r.insert(1, "features", fname); r.insert(2, "scheme", sname)
        r["MAE"], r["bias"] = g.ae.mean(), g.e.mean()
        r["rel_MAE"] = r.MAE / r.x_mean
        rows.append(r.reset_index())
    per_il = pd.concat(rows, ignore_index=True)
    os.makedirs(OUT, exist_ok=True)
    per_il.to_csv(os.path.join(OUT, "per_il_errors_curated.csv"), index=False, float_format="%.4f")

    # ---- per family (mean over the family's points)
    fam_rows = []
    for col in [c for c in df.columns if c.startswith("pred_")]:
        _, mname, fname, sname = col.split("_", 3)
        ae = (df[col] - y).abs()
        for level in ("anion_family", "cation_family"):
            for fam, idx in df.groupby(level).groups.items():
                fam_rows.append(dict(model=mname, features=fname, scheme=sname, level=level, family=fam,
                                     n_ILs=df.loc[idx, "IL"].nunique(), n_points=len(idx),
                                     MAE=ae[idx].mean(), x_mean=y[idx].mean(),
                                     rel_MAE=ae[idx].mean() / y[idx].mean()))
    per_fam = pd.DataFrame(fam_rows)
    per_fam.to_csv(os.path.join(OUT, "per_family_errors_curated.csv"), index=False, float_format="%.4f")
    sel, pf = plot(per_il, per_fam)

    # ---- console summary
    print("\nPooled MAE by model, features and scheme:")
    for col in [c for c in df.columns if c.startswith("pred_")]:
        print(f"  {col[5:]:40s} {np.abs(df[col] - y).mean():.4f}")
    show = ["IL", "n_points", "n_publications", "x_mean", "MAE", "rel_MAE", "bias", "share_disagreeing_points",
            "n_other_ILs_same_cation", "n_other_ILs_same_anion"]
    print("\nWorst 15 ILs, GB dft, leave-one-IL-out:")
    print(sel.sort_values("MAE", ascending=False)[show].head(15).round(3).to_string(index=False))
    print("\nAnion families, GB dft:")
    print(pf.pivot(index="family", columns="scheme", values="MAE").round(4).to_string())
    cf = per_fam[(per_fam.model == "GB") & (per_fam.features == "dft") & (per_fam.level == "cation_family")]
    print("\nCation families, GB dft:")
    print(cf.pivot(index="family", columns="scheme", values="MAE").round(4).to_string())
    print(f"\nTotal {time.time() - t0:.0f}s")


def per_anion(per_il):
    """Leave-one-anion-out MAE per held-out anion (point-weighted over its ILs), counts against dft."""
    lao = per_il[per_il.scheme == "leave_one_anion_out"].assign(w=lambda d: d.MAE * d.n_points)
    g = lao.groupby(["model", "features", "anion"]).agg(w=("w", "sum"), n_points=("n_points", "sum"),
                                                         n_ILs=("IL", "size"), family=("anion_family", "first"))
    g["MAE"] = g.w / g.n_points
    t = g.MAE.unstack("features").reset_index()
    t = t.merge(g.reset_index().drop_duplicates(["model", "anion"])[["model", "anion", "family", "n_ILs",
                                                                   "n_points"]], on=["model", "anion"])
    t["dft_minus_counts"] = t.dft - t.counts
    return t.rename(columns={"counts": "MAE_counts", "dft": "MAE_dft"})


def plot(per_il, per_fam):
    pa = per_anion(per_il)
    pa.to_csv(os.path.join(OUT, "per_anion_errors_curated.csv"), index=False, float_format="%.4f")
    sel = per_il[(per_il.model == "GB") & (per_il.features == "dft") & (per_il.scheme == "leave_one_IL_out")]
    worst = sel.sort_values("MAE", ascending=False).head(20).iloc[::-1]
    fams = list(ANION_FAMILY)
    color = {f: COLORS[i] for i, f in enumerate(fams)}
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(17, 6.4), gridspec_kw={"width_ratios": [1.15, 1, 0.9]})
    a1.barh(range(len(worst)), worst.MAE, color=[color[f] for f in worst.anion_family], height=0.7)
    a1.set_yticks(range(len(worst)))
    a1.set_yticklabels([f"{il}  (n={n}, {p} pub.)" for il, n, p in
                        zip(worst.IL, worst.n_points, worst.n_publications)], fontsize=8)
    a1.set_xlabel("MAE in x(CO$_2$), leave-one-IL-out")
    a1.set_title("(a) 20 worst-predicted ILs (GB, DFT + counts)", fontsize=10, loc="left")
    handles = [plt.Rectangle((0, 0), 1, 1, color=color[f]) for f in fams]
    a1.legend(handles, fams, fontsize=7, frameon=False,
              loc="lower right", title="anion family (a, c)", title_fontsize=7)

    pf = per_fam[(per_fam.model == "GB") & (per_fam.features == "dft") & (per_fam.level == "anion_family")]
    piv = pf.pivot(index="family", columns="scheme", values="MAE").loc[fams[::-1]]
    nil = pf.drop_duplicates("family").set_index("family").n_ILs
    yy = np.arange(len(piv))
    a2.barh(yy + 0.18, piv.leave_one_IL_out, height=0.34, color="#9b9a95", label="leave-one-IL-out")
    a2.barh(yy - 0.18, piv.leave_one_anion_out, height=0.34, color="#3d3c39", label="leave-one-anion-out")
    a2.set_yticks(yy)
    a2.set_yticklabels([f"{f} ({nil[f]} ILs)" for f in piv.index], fontsize=8)
    a2.set_xlabel("MAE in x(CO$_2$)")
    a2.set_title("(b) error by anion family (GB, DFT + counts)", fontsize=10, loc="left")
    a2.legend(fontsize=8, frameon=False, loc="upper right")

    d = pa[pa.model == "GB"].sort_values("dft_minus_counts", ascending=False)
    yy = np.arange(len(d))
    a3.barh(yy, d.dft_minus_counts, color=[color[f] for f in d.family], height=0.7)
    a3.axvline(0, color="#52514e", lw=0.8)
    a3.set_yticks(yy)
    a3.set_yticklabels([f"{a} (n={n})" for a, n in zip(d.anion, d.n_points)], fontsize=8)
    a3.set_xlabel("change in MAE from adding DFT descriptors")
    a3.set_title("(c) DFT effect per held-out anion (GB)", fontsize=10, loc="left")
    a3.text(0.98, 0.02, "left of 0: DFT helps", transform=a3.transAxes, ha="right", fontsize=8, color="#52514e")
    for ax in (a1, a2, a3):
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="x", color="#e6e5e1", lw=0.6)
        ax.set_axisbelow(True)
    fig.tight_layout()
    os.makedirs(os.path.dirname(FIG), exist_ok=True)
    fig.savefig(FIG, dpi=160)
    print(f"wrote {os.path.relpath(FIG, ROOT)}")
    return sel, pf


if __name__ == "__main__":
    main()
