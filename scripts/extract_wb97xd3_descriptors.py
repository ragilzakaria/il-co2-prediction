"""Roadmap steps 2.2 and 2.3: read the wB97X-D3/def2-TZVPD single points of HPC job 576022.

Input: calculations/ions_wb97xd3/outputs/<ion>.out (one single point per ion on its r2SCAN-3c
geometry; inputs from scripts/make_ion_singlepoints.py) and the <ion>__pal4 / __pal16 copies of the
parallel test.

Writes
  data/processed/dft_descriptors_wb97xd3.csv   same columns as dft_descriptors_r2scan3c.csv
                                                (Energy, HOMO, LUMO, gap, dipole), plus run time
  calculations/ions_wb97xd3/parallel_timing.csv  wall time on 1, 4 and 16 cores, speed-up, efficiency

Checks per ion: normal termination, SCF converged, closed shell with the expected charge, and for
anions whether the HOMO is bound (negative), which it was not at r2SCAN-3c for [NO3]-.

Run from the repository root:
    python scripts/extract_wb97xd3_descriptors.py
"""
import re
from pathlib import Path

import pandas as pd

from extract_descriptors import dipole_from_components, last_orbital_block, safe_name

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "calculations" / "ions_wb97xd3"
OUT = ROOT / "data" / "processed" / "dft_descriptors_wb97xd3.csv"
TIMING = RUN / "parallel_timing.csv"
HARTREE_EV = 27.211386


def run_time_s(text):
    m = re.search(r"TOTAL RUN TIME:\s*(\d+) days (\d+) hours (\d+) minutes (\d+) seconds (\d+) msec", text)
    if not m:
        return None
    d, h, mi, s, ms = map(int, m.groups())
    return ((d * 24 + h) * 60 + mi) * 60 + s + ms / 1000


def parse(path):
    text = path.read_text(errors="ignore") if path.exists() else ""
    if "ORCA TERMINATED NORMALLY" not in text:
        return None
    method = re.search(r"\|\s*\d+>\s*!\s*(.+)", text).group(1).strip()
    charge = int(re.search(r"Total Charge\s+Charge\s+\.+\s*(-?\d+)", text).group(1))
    mult = int(re.search(r"Multiplicity\s+Mult\s+\.+\s*(\d+)", text).group(1))
    homo, lumo = last_orbital_block(text)
    dip = re.findall(r"Magnitude \(Debye\)\s*:\s*(\d+\.\d+)", text)
    return {
        "method": method,
        "Energy_Hartree": float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", text)[-1]),
        "HOMO_Eh": homo, "LUMO_Eh": lumo, "HOMO_LUMO_Gap_Eh": lumo - homo,
        "Dipole_Moment_Debye": float(dip[-1]) if dip else dipole_from_components(text),
        "charge": charge, "multiplicity": mult,
        "scf_converged": "SCF NOT CONVERGED" not in text and "SCF CONVERGED AFTER" in text,
        "n_atoms": int(re.search(r"Number of atoms\s*\.+\s*(\d+)", text).group(1)),
        "nprocs": int(m.group(1)) if (m := re.search(r"Program running with (\d+) parallel MPI-processes", text)) else 1,
        "run_time_s": run_time_s(text),
    }


def main():
    man = pd.read_csv(RUN / "manifest.csv")
    rows, missing = [], []
    for r in man.itertuples():
        rec = parse(RUN / "outputs" / f"{r.file}.out")
        if rec is None:
            missing.append(r.file)
            continue
        assert rec["charge"] == r.charge and rec["multiplicity"] == 1, (r.file, rec["charge"], rec["multiplicity"])
        rows.append({"ion_abbreviation": r.ion_abbreviation, "file": r.file, "ion_type": r.ion_type, **rec})
    df = pd.DataFrame(rows)

    # parallel test (step 2.2)
    t = df[df.file.str.contains("__pal") | df.ion_abbreviation.isin(df.loc[df.file.str.contains("__pal"),
                                                                           "ion_abbreviation"])]
    timing = t.pivot_table(index=["ion_abbreviation", "n_atoms"], columns="nprocs", values="run_time_s").reset_index()
    for p in (4, 16):
        if p in timing and 1 in timing:
            timing[f"speedup_{p}"] = timing[1] / timing[p]
            timing[f"efficiency_{p}"] = timing[f"speedup_{p}"] / p
    timing.columns = [f"time_{c}cores_s" if isinstance(c, int) else c for c in timing.columns]
    timing.to_csv(TIMING, index=False, float_format="%.2f")
    # the parallel copies give the same numbers; keep the serial run per ion
    for ion, g in t.groupby("ion_abbreviation"):
        assert g.Energy_Hartree.max() - g.Energy_Hartree.min() < 1e-6, ion
    ions = df[~df.file.str.contains("__pal")].drop(columns=["file", "nprocs"])
    ions.to_csv(OUT, index=False)

    print(f"{len(ions)} of {(~man.file.str.contains('__pal')).sum()} ions -> {OUT.relative_to(ROOT)}")
    print(f"SCF not converged: {ions.loc[~ions.scf_converged, 'ion_abbreviation'].tolist()}")
    an = ions[ions.ion_type == "anion"]
    print(f"anions with an unbound HOMO (> 0): {an.loc[an.HOMO_Eh > 0, 'ion_abbreviation'].tolist()}")
    print(f"missing or unfinished: {len(missing)} {missing[:8]}{' ...' if len(missing) > 8 else ''}")
    print("\nParallel test (step 2.2):")
    print(timing.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
