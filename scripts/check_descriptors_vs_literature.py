"""Check our single-ion DFT descriptors against published computed values.

Compares, ion by ion:
  * -E(HOMO) with published ionisation energies of the ions in an IL-like continuum
    (B3LYP-D3(BJ)/6-311+G(d,p), SMD; PCCP 2025, 27, 9068, ESI Tables S12-S13);
  * the HOMO-LUMO gap with the PBE-D3(BJ)/DZP gaps of the ionic-liquids database
    (Stavarache et al. 2026, Zenodo 10.5281/zenodo.20539350).

Absolute values are not expected to agree (different functionals, gas phase versus
continuum), so the check is on ordering (Spearman rank correlation) and on linearity
(Pearson r). Two descriptor sets are checked:
  current  data/processed/dft_descriptors_r2scan3c.csv  (corrected parser)
  thesis   data/processed/orca_descriptors_individual.csv (values used to train the thesis model)

Ions whose stored outputs were calculated as the wrong species are excluded.

Writes data/processed/checks/descriptors_vs_literature.csv and
images/checks/descriptors_vs_literature.png, and prints a summary.

Run from the repository root:
    python scripts/check_descriptors_vs_literature.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CURRENT = ROOT / "data" / "processed" / "dft_descriptors_r2scan3c.csv"
THESIS = ROOT / "data" / "processed" / "orca_descriptors_individual.csv"
IONS = ROOT / "data" / "ions" / "ions_with_smiles.csv"
REFS = ROOT / "data" / "literature_reference" / "reference_values.csv"
OUT_CSV = ROOT / "data" / "processed" / "checks" / "descriptors_vs_literature.csv"
OUT_PNG = ROOT / "images" / "checks" / "descriptors_vs_literature.png"

HARTREE_EV = 27.211386
# Outputs in calculations/ions/ for these were computed as the wrong species (see commit 0953722).
WRONG_SPECIES = {"[DCA]", "[DBPO4]", "[C12PhSO3]", "[MDEA]", "[HMMIM]", "[P4,4,4,4]"}


def load_ours(path, label):
    d = pd.read_csv(path)
    d = d[~d["ion_abbreviation"].isin(WRONG_SPECIES)].copy()
    d["neg_HOMO_eV"] = -d["HOMO_Eh"] * HARTREE_EV
    d["gap_eV"] = (d["LUMO_Eh"] - d["HOMO_Eh"]) * HARTREE_EV
    d["set"] = label
    return d[["ion_abbreviation", "set", "neg_HOMO_eV", "gap_eV"]].rename(columns={"ion_abbreviation": "ion"})


def stats(x, y):
    n = len(x)
    if n < 4:
        return dict(n=n, spearman=float("nan"), pearson=float("nan"))
    return dict(n=n, spearman=x.rank().corr(y.rank()), pearson=x.corr(y))


def main():
    types = pd.read_csv(IONS).set_index("ion_abbreviation")["ion_type"]
    refs = pd.read_csv(REFS)
    ie = refs[refs["quantity"].str.startswith("ionisation energy")][["ion", "value"]]
    ie = ie.rename(columns={"value": "ref_IE_eV"})
    gap = refs[(refs["category"] == "QM single ion") & (refs["quantity"] == "HOMO-LUMO gap")][["ion", "value"]]
    gap = gap.rename(columns={"value": "ref_gap_PBE_eV"})

    ours = pd.concat([load_ours(CURRENT, "current"), load_ours(THESIS, "thesis")])
    ours["type"] = ours["ion"].map(types)
    table = ours.merge(ie, on="ion", how="left").merge(gap, on="ion", how="left")
    table = table[table["ref_IE_eV"].notna() | table["ref_gap_PBE_eV"].notna()]
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    table.sort_values(["set", "type", "ref_IE_eV"]).to_csv(OUT_CSV, index=False, float_format="%.4f")

    print("Check of single-ion descriptors against literature (rank = Spearman, r = Pearson)\n")
    rows = []
    for s in ("current", "thesis"):
        for t in ("anion", "cation"):
            sub = table[(table["set"] == s) & (table["type"] == t)]
            a = sub.dropna(subset=["ref_IE_eV"])
            b = sub.dropna(subset=["ref_gap_PBE_eV"])
            r1, r2 = stats(a["neg_HOMO_eV"], a["ref_IE_eV"]), stats(b["gap_eV"], b["ref_gap_PBE_eV"])
            rows.append((s, t, r1, r2))
            print(f"{s:8s} {t:7s}  -HOMO vs IE:  n={r1['n']:2d}  rank={r1['spearman']:+.2f}  r={r1['pearson']:+.2f}"
                  f"   |  gap vs PBE gap:  n={r2['n']:2d}  rank={r2['spearman']:+.2f}  r={r2['pearson']:+.2f}")

    # Largest rank disagreements for the current anions
    a = table[(table["set"] == "current") & (table["type"] == "anion")].dropna(subset=["ref_IE_eV"]).copy()
    a["rank_ours"] = a["neg_HOMO_eV"].rank()
    a["rank_ref"] = a["ref_IE_eV"].rank()
    a["rank_shift"] = a["rank_ours"] - a["rank_ref"]
    print("\nCurrent anions, ordered by published IE:")
    print(a.sort_values("ref_IE_eV")[["ion", "ref_IE_eV", "neg_HOMO_eV", "rank_ref", "rank_ours", "rank_shift"]]
          .to_string(index=False, float_format=lambda v: f"{v:.2f}"))

    # Figure: -HOMO vs published IE, anions and cations, current and thesis sets
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for ax, t in zip(axes, ("anion", "cation")):
        for s, m, c in (("thesis", "x", "#B5651D"), ("current", "o", "#1F6F78")):
            sub = table[(table["set"] == s) & (table["type"] == t)].dropna(subset=["ref_IE_eV"])
            ax.scatter(sub["ref_IE_eV"], sub["neg_HOMO_eV"], marker=m, color=c, label=s, s=36)
            if s == "current":
                for _, r in sub.iterrows():
                    ax.annotate(r["ion"].strip("[]"), (r["ref_IE_eV"], r["neg_HOMO_eV"]), fontsize=7,
                                xytext=(3, 3), textcoords="offset points")
        ax.set_xlabel("published IE in IL, B3LYP-D3/SMD (eV)")
        ax.set_ylabel("our -E(HOMO), gas phase (eV)")
        ax.set_title(f"{t}s")
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=160)
    print(f"\nWrote {OUT_CSV.relative_to(ROOT)} and {OUT_PNG.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
