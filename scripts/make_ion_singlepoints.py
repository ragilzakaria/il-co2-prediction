"""Write wB97X-D3/def2-TZVPD single-point inputs on the r2SCAN-3c ion geometries.

Roadmap steps 2.2 and 2.3 (docs/QM_PROTOCOL.md): one single point per ion, on the final geometry
of its r2SCAN-3c optimisation in calculations/ions/<ion>.out. Also writes copies of a few ions
with %pal for the parallel-ORCA timing test.

Output: calculations/ions_wb97xd3/inputs/<ion>.inp, <ion>__pal4.inp, <ion>__pal16.inp,
        calculations/ions_wb97xd3/manifest.csv

Run from the repository root:
    python scripts/make_ion_singlepoints.py
"""
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OPT_DIR = ROOT / "calculations" / "ions"
DEST = ROOT / "calculations" / "ions_wb97xd3"
IONS = ROOT / "data" / "ions" / "ions_with_smiles.csv"

METHOD = "! wB97X-D3 def2-TZVPD def2/J RIJCOSX TightSCF"
MAXCORE = 2500
# Parallel test (step 2.2): a small anion, a small cation and a mid-sized cation.
PAL_TEST = {"[Tf2N]": (4, 16), "[BMIM]": (4, 16), "[P4,4,4,4]": (4, 16)}


def safe_name(abbr):
    # Same rule as scripts/extract_descriptors.py.
    return (abbr.strip("[]").replace("+", "_p").replace("-", "_m")
            .replace("(", "_").replace(")", "_"))


def final_geometry(text):
    block = text.split("CARTESIAN COORDINATES (ANGSTROEM)")[-1].split("\n\n")[0]
    rows = [l.split() for l in block.splitlines()[2:]]
    return [r for r in rows if len(r) == 4]


def write_input(path, title, charge, atoms, nprocs=None):
    lines = [f"# {title}", METHOD, f"%maxcore {MAXCORE}"]
    if nprocs:
        lines.append(f"%pal nprocs {nprocs} end")
    lines.append(f"* xyz {charge} 1")
    lines += [f"{e:<2} {float(x):14.8f} {float(y):14.8f} {float(z):14.8f}" for e, x, y, z in atoms]
    lines.append("*")
    path.write_text("\n".join(lines) + "\n")


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
        title = f"{abbr}: wB97X-D3/def2-TZVPD single point on the r2SCAN-3c geometry of calculations/ions/{name}.out"
        write_input(DEST / "inputs" / f"{name}.inp", title, charge, atoms)
        rows.append({"ion_abbreviation": abbr, "file": name, "ion_type": ion.ion_type,
                     "charge": charge, "n_atoms": n, "nprocs": 1})
        for p in PAL_TEST.get(abbr, ()):
            write_input(DEST / "inputs" / f"{name}__pal{p}.inp", title + f" (parallel test, {p} cores)",
                        charge, atoms, nprocs=p)
            rows.append({"ion_abbreviation": abbr, "file": f"{name}__pal{p}", "ion_type": ion.ion_type,
                         "charge": charge, "n_atoms": n, "nprocs": p})
    man = pd.DataFrame(rows)
    man.to_csv(DEST / "manifest.csv", index=False)
    print(f"{(man.nprocs == 1).sum()} ions, {(man.nprocs > 1).sum()} parallel-test inputs in {DEST.relative_to(ROOT)}")
    print(f"skipped (no converged optimisation): {skipped}")
    print(man[man.nprocs == 1].n_atoms.describe()[["min", "50%", "max"]].to_dict())


if __name__ == "__main__":
    main()
