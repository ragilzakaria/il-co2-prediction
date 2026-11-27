"""Roadmap step 3.2: r2SCAN-3c re-optimisation of the lowest CREST conformers of every ion.

Input: the CREST results of step 3.1, copied back from Mahameru to
calculations/conformers_crest/outputs/<ion>.conformers.xyz (GFN2-xTB energies in Eh on the comment
lines), and the current r2SCAN-3c structure of each ion in calculations/ions/<ion>.out.

Per ion the candidates are
  c0        the current r2SCAN-3c structure (so the result can never be worse than what we have)
  c1..c5    the lowest CREST conformers within 3 kcal/mol (GFN2-xTB), sorted by energy; a conformer is
            skipped as a duplicate when it is within 0.1 kcal/mol of one already chosen and its
            heavy atoms match it within 0.3 A RMSD (mirror images count as matches)
The GFN2-xTB ranking is not reliable to a few kcal/mol, which is why several are re-optimised.

Each candidate gets an r2SCAN-3c optimisation (stage 1 of calculations/conformers_r2scan3c/submit.sh).
The job then keeps the lowest converged candidate per ion and runs its frequencies (stage 2), with
the settings of docs/QM_PROTOCOL.md.

Output: calculations/conformers_r2scan3c/inputs/<ion>__c<k>.inp     optimisations
        calculations/conformers_r2scan3c/freq_templates/<ion>.freq.inp  frequency input; the job
                                                    fills in the geometry of the lowest candidate
        calculations/conformers_r2scan3c/candidates.csv   one row per candidate, with its xTB energy
        calculations/conformers_r2scan3c/ions.txt         "<ion>\t<charge>\t<n_atoms>" for the job

Run from the repository root, after the CREST outputs are copied back:
    python scripts/make_conformer_reopt.py [--crest-dir DIR]
"""
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

from make_ion_singlepoints import final_geometry

ROOT = Path(__file__).resolve().parents[1]
CREST = ROOT / "calculations" / "conformers_crest"
OPT_DIR = ROOT / "calculations" / "ions"
DEST = ROOT / "calculations" / "conformers_r2scan3c"

E_WINDOW = 3.0      # kcal/mol
MAX_CREST = 5
DUP_E = 0.1         # kcal/mol
DUP_RMSD = 0.3      # Angstrom, heavy atoms
HARTREE_KCAL = 627.509
OPT = "! r2SCAN-3c Opt MiniPrint"
FREQ = "! r2SCAN-3c Freq MiniPrint"
MAXCORE = 3000
FREQ_NPROCS = 4


def read_multi_xyz(path):
    """Frames of a multi-structure xyz file: list of (energy, [(el, x, y, z), ...])."""
    lines = path.read_text().splitlines()
    frames, i = [], 0
    while i < len(lines) and lines[i].strip():
        n = int(lines[i])
        energy = float(lines[i + 1].split()[0])
        atoms = [(l.split()[0], *map(float, l.split()[1:4])) for l in lines[i + 2:i + 2 + n]]
        frames.append((energy, atoms))
        i += n + 2
    return frames


def heavy_rmsd(a, b):
    """Heavy-atom RMSD after optimal superposition, mirror images allowed (same atom order)."""
    P = np.array([xyz[1:] for xyz in a if xyz[0] != "H"])
    Q = np.array([xyz[1:] for xyz in b if xyz[0] != "H"])
    P, Q = P - P.mean(0), Q - Q.mean(0)
    best = np.inf
    for mirror in (1, -1):
        Qm = Q * np.array([1, 1, mirror])
        U, S, Vt = np.linalg.svd(P.T @ Qm)
        d = np.sign(np.linalg.det(U @ Vt))
        R = U @ np.diag([1, 1, d]) @ Vt
        best = min(best, np.sqrt(((P @ R - Qm) ** 2).sum(1).mean()))
    return best


def pick(frames):
    frames = sorted(frames, key=lambda f: f[0])
    e0 = frames[0][0]
    chosen = []
    for e, atoms in frames:
        rel = (e - e0) * HARTREE_KCAL
        if rel > E_WINDOW or len(chosen) == MAX_CREST:
            break
        if any(abs(e - ce) * HARTREE_KCAL < DUP_E and heavy_rmsd(atoms, ca) < DUP_RMSD for ce, ca, _ in chosen):
            continue
        chosen.append((e, atoms, rel))
    return chosen


def geometry_block(charge, atoms):
    return [f"* xyz {charge} 1"] + [f"{el:<2} {x:14.8f} {y:14.8f} {z:14.8f}" for el, x, y, z in atoms] + ["*"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--crest-dir", type=Path, default=CREST / "outputs")
    args = ap.parse_args()
    man = pd.read_csv(CREST / "manifest.csv")
    for sub in ("inputs", "freq_templates"):
        (DEST / sub).mkdir(parents=True, exist_ok=True)
    rows, ions, missing = [], [], []
    for r in man.itertuples():
        f = args.crest_dir / f"{r.file}.conformers.xyz"
        if not f.exists():
            missing.append(r.ion_abbreviation)
            continue
        text = (OPT_DIR / f"{r.file}.out").read_text(errors="ignore")
        current = [(e, float(x), float(y), float(z)) for e, x, y, z in final_geometry(text)]
        frames = read_multi_xyz(f)
        assert all([a[0] for a in fr] == [a[0] for a in current] for _, fr in frames), \
            f"{r.file}: CREST changed the atom order or composition"
        cands = [("c0", None, None, current)] + [(f"c{k}", e, rel, atoms)
                                                  for k, (e, atoms, rel) in enumerate(pick(frames), 1)]
        for tag, e, rel, atoms in cands:
            name = f"{r.file}__{tag}"
            src = ("current r2SCAN-3c structure, calculations/ions/" if tag == "c0"
                   else f"CREST conformer {tag[1:]}, GFN2-xTB {rel:.2f} kcal/mol above the lowest, from ")
            lines = [f"# {r.ion_abbreviation} {tag}: r2SCAN-3c optimisation; start: {src}"
                     f"{r.file}{'.out' if tag == 'c0' else '.conformers.xyz'}",
                     OPT, f"%maxcore {MAXCORE}"] + geometry_block(r.charge, atoms)
            (DEST / "inputs" / f"{name}.inp").write_text("\n".join(lines) + "\n")
            rows.append({"ion_abbreviation": r.ion_abbreviation, "file": r.file, "candidate": tag,
                         "xtb_energy_Eh": e, "xtb_rel_kcal": rel, "n_atoms": r.n_atoms, "charge": r.charge})
        (DEST / "freq_templates" / f"{r.file}.freq.inp").write_text("\n".join([
            f"# {r.ion_abbreviation}: r2SCAN-3c frequencies of the lowest re-optimised conformer "
            f"(selected by calculations/conformers_r2scan3c/submit.sh)",
            FREQ, f"%maxcore {MAXCORE}", f"%pal nprocs {FREQ_NPROCS} end",
            f"* xyzfile {r.charge} 1 {r.file}.best.xyz"]) + "\n")
        ions.append(f"{r.file}\t{r.charge}\t{r.n_atoms}\n")
    pd.DataFrame(rows).to_csv(DEST / "candidates.csv", index=False)
    (DEST / "ions.txt").write_text("".join(ions))
    c = pd.DataFrame(rows)
    print(f"{len(ions)} ions, {len(c)} optimisation inputs in {DEST.relative_to(ROOT)}")
    if len(c):
        print("candidates per ion:", c.groupby("file").size().value_counts().sort_index().to_dict())
    if missing:
        print(f"no CREST result yet for {len(missing)} ions: {missing[:10]}{' ...' if len(missing) > 10 else ''}")


if __name__ == "__main__":
    main()
