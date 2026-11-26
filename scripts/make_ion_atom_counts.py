"""Elemental atom counts of every ion, from its optimised ORCA structure.

data/ions/ions_with_atomic_counts.csv was counted from the thesis SMILES, six of which describe
the wrong species ([DCA] with an extra H, [HMMIM] with PF6, [P4,4,4,4] with Br, and [MDEA],
[DBPO4] and [C12PhSO3] off by one H). This script counts the atoms of the corrected structures in
calculations/ions/ instead. [Cl] (no calculation) is a single Cl atom.

Output: data/ions/ion_atom_counts_r2scan3c.csv, same columns as the thesis file.

Run from the repository root:
    python scripts/make_ion_atom_counts.py
"""
from collections import Counter
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT / "data" / "ions" / "ions_with_atomic_counts.csv"
OUT = ROOT / "data" / "ions" / "ion_atom_counts_r2scan3c.csv"


def safe_name(abbr):
    return (abbr.strip("[]").replace("+", "_p").replace("-", "_m")
            .replace("(", "_").replace(")", "_"))


def main():
    old = pd.read_csv(OLD)
    num_cols = [c for c in old.columns if c.startswith("num_")]
    new = old.copy()
    changed = []
    for i, r in old.iterrows():
        f = ROOT / "calculations" / "ions" / f"{safe_name(r.ion_abbreviation)}.out"
        if r.ion_abbreviation == "[Cl]":
            counts = Counter({"Cl": 1})
        else:
            block = f.read_text(errors="ignore").split("CARTESIAN COORDINATES (ANGSTROEM)")[-1].split("\n\n")[0]
            counts = Counter(l.split()[0] for l in block.splitlines()[2:] if len(l.split()) == 4)
        unknown = set(counts) - {c[4:] for c in num_cols}
        assert not unknown, (r.ion_abbreviation, unknown)
        for c in num_cols:
            new.at[i, c] = counts.get(c[4:], 0)
        if (new.loc[i, num_cols].values != old.loc[i, num_cols].values).any():
            changed.append(r.ion_abbreviation)
    new = new.drop(columns=[c for c in ["smiles"] if c in new.columns])
    new.to_csv(OUT, index=False)
    print(f"wrote {OUT.relative_to(ROOT)}; counts changed for {changed}")


if __name__ == "__main__":
    main()
