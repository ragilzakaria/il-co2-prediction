"""Compare our ion-pair interaction energies with the published CCSD(T)-level references.

Reads calculations/benchmark_ion_pairs/manifest.csv and outputs/*.out (from submit.sh), computes

  r2SCAN-3c          E_int = E(pair) - E(cat) - E(an)            (gCP and D4 are built in)
  wB97X-D3           the same, without counterpoise
  wB97X-D3 + CP      E(pair) - E(cat in pair basis) - E(an in pair basis)

and reports the deviation from the reference (method - reference, kJ/mol): mean (MD), mean
absolute (MAD) and largest (MAXD), overall, per reference set and per anion.
Writes calculations/benchmark_ion_pairs/results.csv.

Run from the repository root:
    python scripts/analyse_ionpair_benchmark.py
"""
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "calculations" / "benchmark_ion_pairs"
HARTREE_KJ = 2625.49964


def energy(name):
    f = BENCH / "outputs" / f"{name}.out"
    if not f.exists():
        return None
    text = f.read_text(errors="ignore")
    if "ORCA TERMINATED NORMALLY" not in text:
        return None
    e = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", text)
    return float(e[-1]) if e else None


def main():
    man = pd.read_csv(BENCH / "manifest.csv")
    rows = []
    for _, m in man.iterrows():
        e = {f"{meth}_{part}": energy(f"{m.id}__{meth}__{part}")
             for meth, parts in (("r2scan3c", ["pair", "cat", "an"]), ("wb97xd3", ["pair", "cat", "an", "catG", "anG"]))
             for part in parts}
        r = m.to_dict()

        def eint(p, c, a):
            return None if None in (e[p], e[c], e[a]) else (e[p] - e[c] - e[a]) * HARTREE_KJ

        r["r2SCAN-3c"] = eint("r2scan3c_pair", "r2scan3c_cat", "r2scan3c_an")
        r["wB97X-D3"] = eint("wb97xd3_pair", "wb97xd3_cat", "wb97xd3_an")
        r["wB97X-D3+CP"] = eint("wb97xd3_pair", "wb97xd3_catG", "wb97xd3_anG")
        rows.append(r)
    df = pd.DataFrame(rows)
    methods = ["r2SCAN-3c", "wB97X-D3", "wB97X-D3+CP"]
    for meth in methods:
        df[f"dev_{meth}"] = df[meth] - df["E_ref_kJmol"]
    df.to_csv(BENCH / "results.csv", index=False, float_format="%.2f")

    def summary(sub, label):
        out = [f"{label:28s} n"]
        for meth in methods:
            d = sub[f"dev_{meth}"].dropna()
            if len(d):
                out.append(f"{meth}: n={len(d):2d} MD={d.mean():+6.1f} MAD={d.abs().mean():5.1f} MAXD={d.abs().max():5.1f}")
        print("  ".join(out))

    print("Deviation from the CCSD(T)-level reference, kJ/mol (method - reference)\n")
    summary(df, "all")
    for src, sub in df.groupby("source"):
        summary(sub, src)
    print()
    for an, sub in df.groupby("anion"):
        summary(sub, an)
    print(f"\nWrote {(BENCH / 'results.csv').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
