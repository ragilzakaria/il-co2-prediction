"""σ-moments of every ion from its COSMO file (roadmap step 3.4).

Input: calculations/cosmors/cosmo/<ion>.orcacosmo (BP86/def2-TZVPD, CPCM, from cosmors_henry.py cosmo).
Screening charge densities σ are averaged over 0.5 Å as in COSMO-RS (cosmors_henry.averaged_sigma).
σ in e/nm² (COSMO sign convention: σ > 0 on negative parts of the ion, i.e. hydrogen-bond acceptor
regions; an anion's total screening charge is +1).

Per ion (Klamt's σ-moments, the compact descriptors used in COSMO-RS QSPR):
  area_A2       cavity surface area (Å²)
  volume_A3     cavity volume (Å³)
  M2            Σ A σ²          polarity / electrostatic interaction
  M3            Σ A σ³          asymmetry of the σ-profile
  M4            Σ A σ⁴          weight of strongly polar segments
  Macc          Σ A max(σ − σ_hb, 0)    hydrogen-bond acceptor moment
  Mdon          Σ A max(−σ − σ_hb, 0)   hydrogen-bond donor moment
with σ_hb = 0.85 e/nm² (the usual COSMOtherm threshold).

Output: calculations/cosmors/sigma_moments.csv (ion abbreviation, file name, type, moments)

Run from the repository root:
    python scripts/cosmors_sigma_moments.py
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cosmors_henry import BOHR, DEST, MANIFEST, averaged_sigma, read_cosmo  # noqa: E402

SIGMA_HB = 0.85  # e/nm²
MOMENTS = ["area_A2", "volume_A3", "M2", "M3", "M4", "Macc", "Mdon"]


def moments(path):
    m = read_cosmo(path)
    area, sig = averaged_sigma(m)
    sig = sig * 100.0  # e/Å² -> e/nm²
    return dict(area_A2=area.sum(), volume_A3=m["volume"] * BOHR ** 3,
                M2=(area * sig ** 2).sum(), M3=(area * sig ** 3).sum(), M4=(area * sig ** 4).sum(),
                Macc=(area * np.clip(sig - SIGMA_HB, 0, None)).sum(),
                Mdon=(area * np.clip(-sig - SIGMA_HB, 0, None)).sum())


def main():
    man = pd.read_csv(MANIFEST)
    rows = []
    for _, r in man.iterrows():
        path = DEST / "cosmo" / f"{r.file}.orcacosmo"
        if path.exists():
            rows.append(dict(ion_abbreviation=r.ion_abbreviation, file=r.file, ion_type=r.ion_type,
                             **moments(path)))
    out = pd.DataFrame(rows).sort_values(["ion_type", "Macc"])
    out.to_csv(DEST / "sigma_moments.csv", index=False, float_format="%.4f")
    print(f"{len(out)} ions -> {DEST / 'sigma_moments.csv'}")
    with pd.option_context("display.width", 200, "display.max_rows", 200):
        print(out.drop(columns="file").round(2).to_string(index=False))


if __name__ == "__main__":
    main()
