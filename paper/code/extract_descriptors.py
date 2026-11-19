#!/usr/bin/env python3
"""
Re-extract per-ion DFT descriptors from the ORCA outputs in calculations/ions/.

The outputs are r2SCAN-3c geometry optimisations followed by frequency
calculations (ORCA 6.0). For each ion this script reads, from the *final*
(optimised) structure:

  Energy_Hartree       last "FINAL SINGLE POINT ENERGY"
  HOMO_Eh / LUMO_Eh    highest occupied / lowest virtual orbital of the last
                       "ORBITAL ENERGIES" block (that block only)
  HOMO_LUMO_Gap_Eh     LUMO - HOMO
  Dipole_Moment_Debye  last "Magnitude (Debye)" (total dipole magnitude)

and records whether the optimisation converged and how many imaginary modes the
frequency calculation found. Writes paper/data/dft_descriptors_r2scan3c.csv.

Run from the repository root:
    python paper/code/extract_descriptors.py
"""

import os
import re

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "../.."))
OUT_DIR = os.path.join(ROOT, "calculations", "ions")
IONS = os.path.join(ROOT, "data", "ions_with_smiles.csv")
DEST = os.path.abspath(os.path.join(HERE, "../data/dft_descriptors_r2scan3c.csv"))


def safe_name(abbr):
    # Same file-naming rule as notebooks/06-inp-generation and 07-descriptor-extraction.
    return (abbr.strip("[]").replace("+", "_p").replace("-", "_m")
            .replace("(", "_").replace(")", "_"))


def last_orbital_block(text):
    block = text.split("ORBITAL ENERGIES")[-1].splitlines()
    orbitals, started = [], False
    for line in block[1:]:
        parts = line.split()
        if len(parts) == 4 and parts[0].isdigit():
            orbitals.append((float(parts[1]), float(parts[2])))
            started = True
        elif started:
            break
    homo = max(e for occ, e in orbitals if occ > 0)
    lumo = min(e for occ, e in orbitals if occ == 0)
    return homo, lumo


AU_TO_DEBYE = 2.541746


def dipole_from_components(text):
    # Fallback for runs that stopped before the final property block.
    comps = re.findall(r"Total Dipole Moment\s*:\s*(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", text)
    if not comps:
        return None
    x, y, z = map(float, comps[-1])
    return AU_TO_DEBYE * (x * x + y * y + z * z) ** 0.5


def parse(path):
    text = open(path, errors="ignore").read()
    if "ORCA TERMINATED NORMALLY" not in text:
        return None
    method = re.search(r"\|\s*1>\s*!\s*(.+)", text).group(1).strip()
    energies = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", text)
    dipoles = re.findall(r"Magnitude \(Debye\)\s*:\s*(\d+\.\d+)", text)
    homo, lumo = last_orbital_block(text)
    return {
        "method": method,
        "Energy_Hartree": float(energies[-1]),
        "HOMO_Eh": homo,
        "LUMO_Eh": lumo,
        "HOMO_LUMO_Gap_Eh": lumo - homo,
        "Dipole_Moment_Debye": float(dipoles[-1]) if dipoles else dipole_from_components(text),
        "opt_converged": ("OPTIMIZATION RUN DONE" in text) or ("HAS CONVERGED" in text),
        "n_imaginary": len(re.findall(r"\*\*\*imaginary mode\*\*\*", text)),
        "n_opt_cycles": len(energies),
    }


def main():
    ions = pd.read_csv(IONS)["ion_abbreviation"].dropna().unique()
    rows, missing = [], []
    for abbr in ions:
        path = os.path.join(OUT_DIR, safe_name(abbr) + ".out")
        rec = parse(path) if os.path.exists(path) else None
        if rec is None:
            missing.append(abbr)
            continue
        rows.append({"ion_abbreviation": abbr, **rec})
    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(DEST), exist_ok=True)
    df.to_csv(DEST, index=False)
    print(f"[descriptors] {len(df)} ions written to {DEST}")
    print(f"[descriptors] methods: {df['method'].value_counts().to_dict()}")
    print(f"[descriptors] not converged: {df.loc[~df.opt_converged, 'ion_abbreviation'].tolist()}")
    print(f"[descriptors] with imaginary modes: "
          f"{df.loc[df.n_imaginary > 0, 'ion_abbreviation'].tolist()}")
    print(f"[descriptors] no usable output: {missing}")


if __name__ == "__main__":
    main()
