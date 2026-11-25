"""Rebuild structures and ORCA inputs for ions whose SMILES were wrong.

Six SMILES in data/ions/ions_with_smiles.csv described the wrong species (a neutral
molecule, or the ion together with its counter-ion). Notebook 06 then chose the
multiplicity from the electron count, so ORCA ran them as doublet radicals.

This script:
  1. checks every SMILES in data/ions/ions_with_smiles.csv (one fragment, formal
     charge +1 for cations / -1 for anions, even electron count) and stops if
     any fails;
  2. rebuilds the 3D structure of the fixed ions with the original protocol
     (RDKit ETKDGv3 embedding + MMFF optimisation, now with a fixed seed);
  3. writes calculations/structures/<name>.xyz and calculations/ions/<name>.inp with the
     correct charge and multiplicity 1, coordinates embedded in the input.

Run from the repository root (needs rdkit):
    python scripts/rebuild_fixed_ions.py
"""
import csv
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

ROOT = Path(__file__).resolve().parents[1]
IONS_CSV = ROOT / "data" / "ions" / "ions_with_smiles.csv"
STRUCT_DIR = ROOT / "calculations" / "structures"
INP_DIR = ROOT / "calculations" / "ions"

FIXED = ["[DCA]", "[DBPO4]", "[C12PhSO3]", "[MDEA]", "[HMMIM]", "[P4,4,4,4]"]
SEED = 42


def safe_name(abbr: str) -> str:
    # Same rule as notebooks 05/06.
    return (abbr.strip("[]").replace("+", "_p").replace("-", "_m")
            .replace("(", "_").replace(")", "_"))


def check(row) -> Chem.Mol:
    mol = Chem.MolFromSmiles(row["smiles"])
    if mol is None:
        raise ValueError(f"{row['ion_abbreviation']}: SMILES does not parse")
    want = 1 if row["ion_type"] == "cation" else -1
    charge = Chem.GetFormalCharge(mol)
    frags = len(Chem.GetMolFrags(mol))
    electrons = sum(a.GetAtomicNum() for a in Chem.AddHs(mol).GetAtoms()) - charge
    problems = []
    if frags != 1:
        problems.append(f"{frags} fragments")
    if charge != want:
        problems.append(f"charge {charge:+d}, expected {want:+d}")
    if electrons % 2:
        problems.append(f"odd electron count ({electrons})")
    if problems:
        raise ValueError(f"{row['ion_abbreviation']}: " + "; ".join(problems))
    return mol


def main() -> None:
    rows = list(csv.DictReader(open(IONS_CSV, encoding="utf-8")))
    mols = {r["ion_abbreviation"]: (r, check(r)) for r in rows}
    print(f"All {len(rows)} SMILES pass the charge/fragment/electron checks.")

    for abbr in FIXED:
        row, mol = mols[abbr]
        charge = Chem.GetFormalCharge(mol)
        molh = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = SEED
        if AllChem.EmbedMolecule(molh, params) != 0:
            raise RuntimeError(f"{abbr}: embedding failed")
        AllChem.MMFFOptimizeMolecule(molh, maxIters=2000)

        name = safe_name(abbr)
        xyz = Chem.MolToXYZBlock(molh)
        (STRUCT_DIR / f"{name}.xyz").write_text(xyz)
        coords = "".join(xyz.splitlines(keepends=True)[2:])
        (INP_DIR / f"{name}.inp").write_text(
            "! r2SCAN-3c Opt Freq MiniPrint\n"
            f"* xyz {charge} 1\n{coords}*\n"
        )
        print(f"{abbr:12} {CalcMolFormula(molh):14} charge {charge:+d}, "
              f"mult 1, {molh.GetNumAtoms()} atoms -> {name}.inp")


if __name__ == "__main__":
    main()
