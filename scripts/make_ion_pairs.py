"""Roadmap step 3.5: starting structures for the ion pairs and anion-CO2 complexes.

Systems
  pairs        the 110 cation-anion pairs of the curated dataset (charge 0)
  CO2          each of the 23 anions with one CO2 molecule (charge -1)

For each system the guest (anion, or CO2) is placed around the host (cation, or anion) at many
positions, in the spirit of the systematic placements used for the IL-2013 set:
  - along 12 directions spread evenly around the host's centre;
  - next to every hydrogen near the charge centre of a cation (H within 3.5 A of N, P or S, and
    all N-H and O-H hydrogens), in the direction of its bond;
  - next to every heavy atom of an anion host (CO2 complexes), pointing outwards.
The guest gets a random orientation (fixed seed) and is moved out along the direction until no
contact is shorter than 2.0 A (to H) or 2.6 A (heavy atoms). Starts whose guest centres lie within
1.0 A of each other are merged; at most 30 starts per system.

Monomer geometries: --geometries current (default) uses calculations/ions/<ion>.out;
--geometries reopt uses the lowest conformers of roadmap step 3.2,
calculations/conformers_r2scan3c/outputs/<ion>.best.xyz. CO2: r2SCAN-3c is not needed for the
start; a linear O=C=O with C=O 1.16 A is used.

The three stages of calculations/ion_pairs/submit.sh then optimise every start with GFN2-xTB,
re-optimise the lowest with r2SCAN-3c and compute the wB97X-D3/def2-TZVPD interaction energy.

Output: calculations/ion_pairs/starts/<system>.xyz   all starts of a system (multi-structure xyz)
        calculations/ion_pairs/systems.tsv           system, kind, host, guest, charge, n_host_atoms,
                                                     n_atoms, n_starts (read by the job script)

Run from the repository root:
    python scripts/make_ion_pairs.py [--geometries current|reopt] [--only SYSTEM ...]
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from make_ion_singlepoints import final_geometry, safe_name

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "calculations" / "ion_pairs"
DATA = ROOT / "data" / "processed" / "curated_co2_solubility_dataset.csv"
MANIFEST = ROOT / "calculations" / "conformers_crest" / "manifest.csv"
SEED = 42
N_DIRECTIONS = 12
MAX_STARTS = 30
CONTACT_H, CONTACT_HEAVY = 2.0, 2.6
CO2 = [("C", 0.0, 0.0, 0.0), ("O", 0.0, 0.0, 1.16), ("O", 0.0, 0.0, -1.16)]
COVALENT = {"H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07, "S": 1.05, "Cl": 1.02}


def load_geometry(abbr, source):
    name = safe_name(abbr)
    if source == "reopt":
        lines = (ROOT / "calculations" / "conformers_r2scan3c" / "outputs" / f"{name}.best.xyz").read_text().splitlines()
        n = int(lines[0])
        return [(l.split()[0], *map(float, l.split()[1:4])) for l in lines[2:2 + n]]
    text = (ROOT / "calculations" / "ions" / f"{name}.out").read_text(errors="ignore")
    return [(e, float(x), float(y), float(z)) for e, x, y, z in final_geometry(text)]


def xyz(atoms):
    return np.array([a[1:] for a in atoms], float)


def fibonacci(n):
    i = np.arange(n) + 0.5
    phi, theta = np.arccos(1 - 2 * i / n), np.pi * (1 + 5 ** 0.5) * i
    return np.c_[np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)]


def random_rotation(rng):
    q = rng.normal(size=4)
    a, b, c, d = q / np.linalg.norm(q)
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d]])


def bonded_parent(atoms, i):
    P = xyz(atoms)
    d = np.linalg.norm(P - P[i], axis=1)
    d[i] = np.inf
    j = int(np.argmin(d))
    return j


def sites(host, kind):
    """(anchor point, unit direction) pairs where the guest is placed."""
    P = xyz(host)
    centre = P.mean(0)
    out = [(centre, d) for d in fibonacci(N_DIRECTIONS)]
    els = [a[0] for a in host]
    if kind == "pair":   # cation host: hydrogens near the charge centre, and acidic hydrogens
        charged = [i for i, e in enumerate(els) if e in ("N", "P", "S")]
        for i, e in enumerate(els):
            if e != "H":
                continue
            j = bonded_parent(host, i)
            near = any(np.linalg.norm(P[i] - P[k]) <= 3.5 for k in charged)
            if near or els[j] in ("N", "O"):
                v = P[i] - P[j]
                out.append((P[i], v / np.linalg.norm(v)))
    else:                # anion host: every heavy atom, outwards
        for i, e in enumerate(els):
            if e != "H":
                v = P[i] - centre
                n = np.linalg.norm(v)
                if n > 1e-6:
                    out.append((P[i], v / n))
    return out


def place(host, guest, anchor, direction, rng):
    H, G = xyz(host), xyz(guest)
    G = (G - G.mean(0)) @ random_rotation(rng).T
    hH = np.array([a[0] == "H" for a in host])
    gH = np.array([a[0] == "H" for a in guest])
    limit = np.where(hH[:, None] | gH[None, :], CONTACT_H, CONTACT_HEAVY)
    for t in np.arange(0.0, 30.0, 0.1):
        c = anchor + direction * t
        D = np.linalg.norm(H[:, None, :] - (G + c)[None, :, :], axis=2)
        if (D >= limit).all():
            return G + c, c
    raise RuntimeError("could not place the guest")


def topology_ok(atoms, n_host):
    """No covalent contact between host and guest (sum of covalent radii + 0.45 A)."""
    P = xyz(atoms)
    r = np.array([COVALENT[a[0]] for a in atoms])
    D = np.linalg.norm(P[:n_host, None] - P[None, n_host:], axis=2)
    return (D > r[:n_host, None] + r[None, n_host:] + 0.45).all()


def write_starts(path, title, host, guest, kind, rng):
    starts, centres = [], []
    for anchor, d in sites(host, kind):
        G, c = place(host, guest, anchor, d, rng)
        if any(np.linalg.norm(c - o) < 1.0 for o in centres):
            continue
        atoms = host + [(guest[k][0], *G[k]) for k in range(len(guest))]
        assert topology_ok(atoms, len(host))
        starts.append(atoms)
        centres.append(c)
    if len(starts) > MAX_STARTS:       # keep an even spread: every k-th start
        idx = np.linspace(0, len(starts) - 1, MAX_STARTS).round().astype(int)
        starts = [starts[i] for i in idx]
    lines = []
    for k, atoms in enumerate(starts):
        lines += [str(len(atoms)), f"{title} start {k}"]
        lines += [f"{e:<2} {x:14.8f} {y:14.8f} {z:14.8f}" for e, x, y, z in atoms]
    path.write_text("\n".join(lines) + "\n")
    return len(starts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--geometries", choices=["current", "reopt"], default="current")
    ap.add_argument("--only", nargs="*", help="system names to write (default: all)")
    args = ap.parse_args()
    man = pd.read_csv(MANIFEST).set_index("ion_abbreviation")
    data = pd.read_csv(DATA, usecols=["cation", "anion"]).drop_duplicates()
    systems = [("pair", f"{safe_name(c)}__{safe_name(a)}", c, a, 0) for c, a in data.itertuples(index=False)]
    systems += [("co2", f"{safe_name(a)}__CO2", a, "CO2", -1) for a in sorted(data.anion.unique())]
    if args.only:
        systems = [s for s in systems if s[1] in args.only]
    (DEST / "starts").mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    rows, skipped = [], []
    cache = {}
    for kind, name, host_abbr, guest_abbr, charge in systems:
        try:
            host = cache.setdefault(host_abbr, load_geometry(host_abbr, args.geometries))
            guest = CO2 if guest_abbr == "CO2" else cache.setdefault(guest_abbr, load_geometry(guest_abbr, args.geometries))
        except FileNotFoundError:
            skipped.append(name)
            continue
        if kind == "pair":
            assert man.loc[host_abbr, "charge"] == 1 and man.loc[guest_abbr, "charge"] == -1, name
        n = write_starts(DEST / "starts" / f"{name}.xyz", f"{host_abbr}{guest_abbr}", host, guest, kind, rng)
        rows.append(dict(system=name, kind=kind, host=host_abbr, guest=guest_abbr, charge=charge,
                         n_host_atoms=len(host), n_atoms=len(host) + len(guest), n_starts=n))
    tab = pd.DataFrame(rows)
    tab.to_csv(DEST / "systems.tsv", sep="\t", index=False)
    print(f"{len(tab)} systems ({(tab.kind == 'pair').sum()} ion pairs, {(tab.kind == 'co2').sum()} anion-CO2), "
          f"{tab.n_starts.sum()} starts; geometries: {args.geometries}")
    print(f"starts per system: min {tab.n_starts.min()}, median {tab.n_starts.median():.0f}, max {tab.n_starts.max()}")
    if skipped:
        print(f"skipped (no geometry): {skipped}")


if __name__ == "__main__":
    main()
