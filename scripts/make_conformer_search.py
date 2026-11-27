"""Write CREST conformer-search inputs (GFN2-xTB) for every computed ion. Roadmap step 3.1.

Every ion so far is a single structure built from its SMILES and optimised with r2SCAN-3c. For the
flexible ions (alkyl chains, ether and hydroxy groups) that is probably not the lowest-energy
conformer. This script writes the r2SCAN-3c geometry of each ion (calculations/ions/<ion>.out) as
the start structure of a CREST search (iMTD-GC workflow, GFN2-xTB, gas phase), plus a manifest with
charge and size. [Cl] (one atom) is skipped.

Why CREST and not ORCA's GOAT: GOAT on one core needed about 5 h for [BMIM] (25 atoms), and the
parallel ORCA it would need on Mahameru was not yet tested; CREST runs in parallel with OpenMP.
Program: CREST 3.0.2, official static build, installed in ~/software/crest-3.0.2 on Mahameru.
Cite: Pracht, Bohle, Grimme, Phys. Chem. Chem. Phys. 22, 7169 (2020), doi:10.1039/C9CP06869D;
Pracht et al., J. Chem. Phys. 160, 114110 (2024), doi:10.1063/5.0197592.

The GFN2-xTB ranking is only a first filter: step 3.2 re-optimises the lowest conformers with
r2SCAN-3c and keeps the lowest at that level. Isolated ions with long chains may fold in the gas
phase; this is checked when the results come back.

Output: calculations/conformers_crest/inputs/<ion>.xyz, calculations/conformers_crest/manifest.csv
        (file, charge, n_atoms, crest_flags), calculations/conformers_crest/jobs.txt (for submit.sh)

CREST settings: the default iMTD-GC search for ions up to 50 atoms; the --quick search (fewer and
shorter metadynamics runs) for the six larger ions, whose default search would not fit the 23 h
limit of the short partition.

Run from the repository root:
    python scripts/make_conformer_search.py
"""
import re
from pathlib import Path

import pandas as pd

from make_ion_singlepoints import final_geometry, safe_name

ROOT = Path(__file__).resolve().parents[1]
OPT_DIR = ROOT / "calculations" / "ions"
DEST = ROOT / "calculations" / "conformers_crest"
IONS = ROOT / "data" / "ions" / "ions_with_smiles.csv"
QUICK_ABOVE = 50   # atoms


def main():
    ions = pd.read_csv(IONS).drop_duplicates("ion_abbreviation")
    (DEST / "inputs").mkdir(parents=True, exist_ok=True)
    rows, skipped = [], []
    for _, ion in ions.iterrows():
        abbr, name = ion.ion_abbreviation, safe_name(ion.ion_abbreviation)
        out = OPT_DIR / f"{name}.out"
        text = out.read_text(errors="ignore") if out.exists() else ""
        if "THE OPTIMIZATION HAS CONVERGED" not in text:
            skipped.append(abbr)
            continue
        atoms = final_geometry(text)
        n = int(re.search(r"Number of atoms\s*\.+\s*(\d+)", text).group(1))
        assert len(atoms) == n, (abbr, len(atoms), n)
        charge = 1 if ion.ion_type == "cation" else -1
        xyz = [str(n), f"{abbr} charge {charge:+d}; r2SCAN-3c geometry of calculations/ions/{name}.out"]
        xyz += [f"{e:<2} {float(x):14.8f} {float(y):14.8f} {float(z):14.8f}" for e, x, y, z in atoms]
        (DEST / "inputs" / f"{name}.xyz").write_text("\n".join(xyz) + "\n")
        rows.append({"ion_abbreviation": abbr, "file": name, "ion_type": ion.ion_type,
                     "charge": charge, "n_atoms": n,
                     "crest_flags": f"--gfn2 --chrg {charge}" + (" --quick" if n > QUICK_ABOVE else "")})
    man = pd.DataFrame(rows).sort_values("n_atoms", ascending=False)
    man.to_csv(DEST / "manifest.csv", index=False)
    # read by submit.sh: one line per ion, largest first, "<file>\t<crest flags>"
    (DEST / "jobs.txt").write_text("".join(f"{r.file}\t{r.crest_flags}\n" for r in man.itertuples()))
    print(f"{len(man)} CREST inputs in {DEST.relative_to(ROOT)}; skipped: {skipped}")
    print("quick mode:", man[man.n_atoms > QUICK_ABOVE][["ion_abbreviation", "n_atoms"]].values.tolist())


if __name__ == "__main__":
    main()
