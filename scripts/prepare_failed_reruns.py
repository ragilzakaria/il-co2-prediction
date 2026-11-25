"""Write new ORCA inputs for the three ions whose r2SCAN-3c runs did not give a minimum.

  [PF6]-      the optimisation hit the 50-cycle limit and one P-F bond broke (2.87 A,
              PF5 + F-). Restart from an ideal octahedron, P-F = 1.60 A.
  [HSO4]-     ended at a saddle point (one imaginary mode, -101 cm-1, the O-H torsion).
              Restart from the old final geometry with the H rotated by 60 degrees about
              the S-O(H) bond.
  [C3F7CO2]-  ended with one tiny imaginary mode (-7 cm-1, an almost free CF3 rotation).
              TightOpt alone kept it (-7.2 cm-1), so restart from the old final geometry
              with TightOpt and the finer DefGrid3 integration grid.

The old outputs are moved to calculations/ions/old_failed/ (by git mv, outside this
script). The new inputs embed their coordinates, as the inputs written by
rebuild_fixed_ions.py do, and keep the protocol: ! r2SCAN-3c Opt Freq.

Run from the repository root:
    python scripts/prepare_failed_reruns.py
"""
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
IONS = ROOT / "calculations" / "ions"
OLD = IONS / "old_failed"
STRUCT = ROOT / "calculations" / "structures"


def last_geometry(out_path):
    text = out_path.read_text()
    blocks = re.findall(r"CARTESIAN COORDINATES \(ANGSTROEM\)\n-+\n(.*?)\n\n", text, re.S)
    atoms = []
    for line in blocks[-1].strip().splitlines():
        el, x, y, z = line.split()[:4]
        atoms.append([el, np.array([float(x), float(y), float(z)])])
    return atoms


def write_input(name, keywords, atoms, charge=-1, comment=""):
    lines = [f"! {keywords}"]
    if comment:
        lines.insert(0, f"# {comment}")
    lines.append(f"* xyz {charge} 1")
    lines += [f"{el:2s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}" for el, p in atoms]
    lines.append("*")
    (IONS / f"{name}.inp").write_text("\n".join(lines) + "\n")
    xyz = [str(len(atoms)), comment or name] + [f"{el:2s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}" for el, p in atoms]
    (STRUCT / f"{name}.xyz").write_text("\n".join(xyz) + "\n")
    print(f"wrote calculations/ions/{name}.inp and calculations/structures/{name}.xyz")


def rotate_about(p, a, b, angle_deg):
    """Rotate point p about the axis through a and b."""
    k = (b - a) / np.linalg.norm(b - a)
    v = p - a
    t = np.radians(angle_deg)
    return a + v * np.cos(t) + np.cross(k, v) * np.sin(t) + k * np.dot(k, v) * (1 - np.cos(t))


def main():
    # PF6-: ideal octahedron
    r = 1.60
    pf6 = [["P", np.zeros(3)]] + [["F", s * r * e] for e in np.eye(3) for s in (1, -1)]
    write_input("PF6", "r2SCAN-3c Opt Freq MiniPrint", pf6,
                comment="octahedral start (P-F 1.60 A); old run dissociated, see old_failed/PF6.out")

    # HSO4-: rotate the acidic H about the S-O(H) bond
    atoms = last_geometry(OLD / "HSO4.out")
    pos = {i: a for i, a in enumerate(atoms)}
    s = next(i for i, a in pos.items() if a[0] == "S")
    h = next(i for i, a in pos.items() if a[0] == "H")
    o = min((i for i, a in pos.items() if a[0] == "O"), key=lambda i: np.linalg.norm(pos[i][1] - pos[h][1]))
    atoms[h][1] = rotate_about(atoms[h][1], atoms[s][1], atoms[o][1], 60.0)
    write_input("HSO4", "r2SCAN-3c Opt Freq MiniPrint", atoms,
                comment="old final geometry with H rotated 60 deg about S-O(H); old run was a saddle point")

    # C3F7CO2-: old final geometry, tighter optimisation
    atoms = last_geometry(OLD / "C3F7CO2.out")
    write_input("C3F7CO2", "r2SCAN-3c TightOpt Freq DefGrid3 MiniPrint", atoms,
                comment="old final geometry, TightOpt + DefGrid3; old run had one -7 cm-1 mode")


if __name__ == "__main__":
    main()
