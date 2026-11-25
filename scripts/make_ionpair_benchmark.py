"""Write ORCA inputs to benchmark our DFT methods on published ion-pair interaction energies.

Reference sets (geometries in data/literature_reference/geometries/):
  karu2018  DLPNO-CCSD(T)/CBS interaction energies on B2PLYP-D3/def2-TZVPP geometries
            (Karu et al. 2018, doi:10.1002/qua.25582); [BMIM], [BMPYR], [C4PY] x 10 anions
  il2013    CCSD(T)/CBS, counterpoise corrected (Zahn et al. 2013, doi:10.1039/C3CP51682B);
            the most stable structure of each pair whose ions are both in our ion list

For every structure, single points at the published geometry (no optimisation):
  r2scan3c  ! r2SCAN-3c                                 pair, cation, anion
  wb97xd3   ! wB97X-D3 def2-TZVPD (RIJCOSX)             pair, cation, anion, and both ions in the
                                                         full pair basis (ghost atoms) for counterpoise

E_int = E(pair) - E(cation) - E(anion), ions frozen at the pair geometry, as in both references.

Writes calculations/benchmark_ion_pairs/inputs/*.inp and calculations/benchmark_ion_pairs/manifest.csv.
Run from the repository root (needs rdkit):
    python scripts/make_ionpair_benchmark.py
"""
import csv
import json
import re
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from ionmatch import ion_table, skeleton  # noqa: E402

GEOM = ROOT / "data" / "literature_reference" / "geometries"
KARU_REF = ROOT / "data" / "literature_reference" / "raw" / "karu2018_pairs.csv"
IONS = ROOT / "data" / "ions" / "ions_with_smiles.csv"
OUT = ROOT / "calculations" / "benchmark_ion_pairs"

METHODS = {
    "r2scan3c": ("! r2SCAN-3c TightSCF", ["pair", "cat", "an"]),
    "wb97xd3": ("! wB97X-D3 def2-TZVPD def2/J RIJCOSX TightSCF", ["pair", "cat", "an", "catG", "anG"]),
}
MAXCORE = 2500


def read_xyz_block(lines):
    atoms = []
    for line in lines:
        el, x, y, z = line.split()[:4]
        atoms.append((el.strip().capitalize(), float(x), float(y), float(z)))
    return atoms


def split_ions(atoms, table, types):
    """Return (cation_name, cation_atom_indices, anion_name, anion_atom_indices) or None."""
    block = "\n".join([str(len(atoms)), ""] + [f"{a} {x:.6f} {y:.6f} {z:.6f}" for a, x, y, z in atoms])
    mol = Chem.MolFromXYZBlock(block)
    rdDetermineBonds.DetermineConnectivity(mol)
    frags = Chem.GetMolFrags(mol, asMols=False)
    if len(frags) != 2:
        return None
    found = {}
    for idx in frags:
        sub = Chem.RWMol(Chem.Mol(mol))
        for i in sorted(set(range(mol.GetNumAtoms())) - set(idx), reverse=True):
            sub.RemoveAtom(i)
        sub = sub.GetMol()
        sub.UpdatePropertyCache(strict=False)
        name = table.get(skeleton(sub))
        if name is None:
            return None
        found[types[name]] = (name, list(idx))
    if set(found) != {"cation", "anion"}:
        return None
    return found["cation"][0], found["cation"][1], found["anion"][0], found["anion"][1]


def write_inputs(sid, atoms, cat_idx, an_idx):
    for method, (keywords, parts) in METHODS.items():
        for part in parts:
            if part == "pair":
                charge, use = 0, [(a, False) for a in range(len(atoms))]
            elif part in ("cat", "catG"):
                charge = 1
                use = [(a, False) for a in cat_idx] + ([(a, True) for a in an_idx] if part == "catG" else [])
            else:
                charge = -1
                use = [(a, False) for a in an_idx] + ([(a, True) for a in cat_idx] if part == "anG" else [])
            lines = [keywords, f"%maxcore {MAXCORE}", f"* xyz {charge} 1"]
            for a, ghost in use:
                el, x, y, z = atoms[a]
                lines.append(f"{el + (':' if ghost else ''):3s} {x:14.8f} {y:14.8f} {z:14.8f}")
            lines.append("*")
            (OUT / "inputs" / f"{sid}__{method}__{part}.inp").write_text("\n".join(lines) + "\n")


def main():
    table = ion_table()
    types = {r["ion_abbreviation"]: r["ion_type"] for r in csv.DictReader(open(IONS))}
    (OUT / "inputs").mkdir(parents=True, exist_ok=True)
    manifest = []

    karu_ref = {r["label"]: float(r["E_int_DLPNO-CCSD(T)/CBS_kJmol"]) for r in csv.DictReader(open(KARU_REF))}
    for f in sorted((GEOM / "karu2018").glob("*.xyz")):
        lines = f.read_text().splitlines()
        atoms = read_xyz_block(lines[2:2 + int(lines[0])])
        res = split_ions(atoms, table, types)
        if res is None:
            print(f"[skip] {f.name}: could not identify both ions")
            continue
        cat, ci, an, ai = res
        sid = f"karu_{f.stem}"
        write_inputs(sid, atoms, ci, ai)
        manifest.append(dict(id=sid, source="karu2018", cation=cat, anion=an, n_atoms=len(atoms),
                             E_ref_kJmol=karu_ref[f.stem], ref_method="DLPNO-CCSD(T)/CBS"))

    best = {}
    lines = (GEOM / "il2013" / "IL-2013-set.xyz").read_text().splitlines()
    i = k = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i].split()[0])
        k += 1
        e_ref = float(re.search(r"(-?\d+\.\d+)", lines[i + 1]).group(1))
        atoms = read_xyz_block(lines[i + 2:i + 2 + n])
        i += 2 + n
        res = split_ions(atoms, table, types)
        if res is None:
            continue
        key = (res[0], res[2])
        if key not in best or e_ref < best[key][0]:
            best[key] = (e_ref, k, atoms, res)
    for (cat, an), (e_ref, k, atoms, res) in sorted(best.items()):
        sid = f"il2013_{k:03d}"
        write_inputs(sid, atoms, res[1], res[3])
        manifest.append(dict(id=sid, source="il2013", cation=cat, anion=an, n_atoms=len(atoms),
                             E_ref_kJmol=e_ref, ref_method="CCSD(T)/CBS, counterpoise corrected"))

    with open(OUT / "manifest.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(manifest[0]))
        w.writeheader()
        w.writerows(manifest)
    n_inp = len(list((OUT / "inputs").glob("*.inp")))
    print(f"{len(manifest)} structures ({sum(m['source'] == 'karu2018' for m in manifest)} Karu, "
          f"{sum(m['source'] == 'il2013' for m in manifest)} IL-2013), {n_inp} inputs; "
          f"atoms {min(m['n_atoms'] for m in manifest)}-{max(m['n_atoms'] for m in manifest)}")


if __name__ == "__main__":
    main()
