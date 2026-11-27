"""Roadmap step 3.4: CO2 Henry constants in ionic liquids from openCOSMO-RS (ORCA 6.1.1).

ORCA's "! COSMORS(solvent)" runs the BP86/def2-TZVPD gas-phase and CPCM single points that
openCOSMO-RS needs and writes them to .orcacosmo files; the COSMO-RS step itself is the separate
program openCOSMORS. The macOS ORCA build does not ship openCOSMORS, the Mahameru build does
(/mgpfs/apps/nuclear/apps/orca/6.1.1/openCOSMORS; a run takes well under a second, so it is run on
the login node). The workflow:

  cosmo    ORCA COSMO step on the laptop for CO2 and every ion that is needed (current r2SCAN-3c
           geometries from calculations/ions/<ion>.out; CO2 optimised with r2SCAN-3c)
           -> calculations/cosmors/cosmo/<ion>.orcacosmo
  prepare  one openCOSMORS input per ionic liquid -> calculations/cosmors/rs/<IL>.json
  (Mahameru login node: calculations/cosmors/run_rs.sh -> rs/<IL>_out.json)
  results  Henry constants at 25 C -> calculations/cosmors/henry.csv; against Palomar et al. 2011
           (data/literature_reference/raw/palomar2011_cosmors.csv) -> palomar_comparison.csv, and
           against P/x of the points up to 1 MPa at 293-303 K in the curated dataset
           -> dataset_comparison.csv

Ionic liquid as solvent ([C+A] model of Palomar et al.): openCOSMORS gives the solvation free
energy (reference state type 4) only for a pure solvent (one component at mole fraction 1), so the
cation and anion are merged into one pseudo-molecule, the anion moved 100 A away so that the
sigma averaging (radius 0.5 A) never mixes the two surfaces. Mole fractions are then per ion pair,
as in the experimental Henry constants.

Henry constant from the solvation free energy. openCOSMORS returns dGsolv in kcal/mol for the
standard state 1 mol/L gas -> 1 mol/L solution: it subtracts RT ln(V_gas(1 atm) / V_m) from the
mole-fraction quantity, with V_m the solvent molar volume from its own QSPR model (atoms, cavity
area, sigma moments; code/core_functions.hpp of github.com/TUHH-TVT/openCOSMO-RS_cpp). So
    K_H (pressure / mole fraction) = RT / V_m * exp(dGsolv / RT),
with V_m from the same QSPR (reimplemented in qspr_volume; it is additive, so for the merged
pseudo-molecule it is the molar volume of the ion pair).

Run from the repository root:
    python scripts/cosmors_henry.py cosmo [ION ...]   (default: every ion of the systems below)
    python scripts/cosmors_henry.py prepare           (the 110 dataset pairs and 18 Palomar ILs)
    python scripts/cosmors_henry.py results
"""
import argparse
import json
import math
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "calculations" / "cosmors"
COSMO, RS = DEST / "cosmo", DEST / "rs"
IONS = ROOT / "calculations" / "ions"
MANIFEST = ROOT / "calculations" / "conformers_crest" / "manifest.csv"
PALOMAR = ROOT / "data" / "literature_reference" / "raw" / "palomar2011_cosmors.csv"
ORCA = os.path.expanduser("~/orca_6_1_1/orca")
BOHR = 0.529177249
R = 8.314462618
SHIFT_A = 100.0


# ---------------------------------------------------------------- ORCA COSMO step (laptop)

def final_geometry(text):
    block = text.split("CARTESIAN COORDINATES (ANGSTROEM)")[-1].split("\n\n")[0]
    return [r for r in (l.split() for l in block.splitlines()[2:]) if len(r) == 4]


def run_cosmo(name, charge, atoms):
    work = DEST / "_work" / name
    work.mkdir(parents=True, exist_ok=True)
    (work / "x.inp").write_text("! COSMORS(WATER)\n* xyz %d 1\n" % charge
                                + "".join("%s %s %s %s\n" % tuple(a) for a in atoms) + "*\n")
    with open(work / "x.out", "w") as f:   # the openCOSMORS call at the end fails on macOS; expected
        subprocess.run([ORCA, "x.inp"], cwd=work, stdout=f, stderr=subprocess.STDOUT)
    src = work / "x.solute.orcacosmo"
    if not src.exists():
        return "%s FAILED (see %s)" % (name, work / "x.out")
    COSMO.mkdir(parents=True, exist_ok=True)
    text = src.read_text().replace("x.solute_cpcm calculation", "%s calculation" % name, 1)
    (COSMO / ("%s.orcacosmo" % name)).write_text(text)
    e_gas = re.findall(r"FINAL SINGLE POINT ENERGY \(Solute-gas-phase\)\s+(-?\d+\.\d+)",
                       (work / "x.out").read_text())[0]
    (COSMO / ("%s.egas" % name)).write_text(e_gas + "\n")
    subprocess.run(["rm", "-rf", str(work)])
    return "%s done" % name


def cmd_cosmo(ions):
    man = pd.read_csv(MANIFEST).set_index("file")
    jobs = []
    if "CO2" in ions:
        jobs.append(("CO2", 0, [l.split()[:4] for l in (DEST / "CO2_r2scan3c.xyz").read_text().splitlines()[2:]
                               if l.strip()]))
    for ion in ions:
        if ion == "CO2" or (COSMO / ("%s.orcacosmo" % ion)).exists():
            continue
        jobs.append((ion, int(man.loc[ion, "charge"]), final_geometry((IONS / ("%s.out" % ion)).read_text())))
    jobs = [j for j in jobs if not (COSMO / ("%s.orcacosmo" % j[0])).exists()]
    print("%d COSMO runs" % len(jobs))
    with ThreadPoolExecutor(max_workers=os.cpu_count()) as ex:
        for msg in ex.map(lambda j: run_cosmo(*j), jobs):
            print(msg, flush=True)


# ---------------------------------------------------------------- .orcacosmo handling

def read_cosmo(path):
    L = Path(path).read_text().split("\n")
    i_xyz = L.index("#XYZ_FILE")
    n_at = int(L[i_xyz + 1])
    i_c = L.index("#COSMO")

    def hval(tag):
        return next(float(l.split()[0]) for l in L[i_c + 1:] if "# " + tag in l)
    i_atoms = next(k for k, l in enumerate(L) if "CARTESIAN COORDINATES (A.U.)" in l) + 2
    i_surf = next(k for k, l in enumerate(L) if "SURFACE POINTS (A.U.)" in l) + 3
    n_s = int(hval("Number of surface points"))
    i_corr = L.index("#COSMO_corrected")
    return dict(lines=L, n_at=n_at, xyz=L[i_xyz + 3:i_xyz + 3 + n_at], atoms=L[i_atoms:i_atoms + n_at],
                surf=[l.split() for l in L[i_surf:i_surf + n_s]],
                corr_q=[float(x) for x in L[i_corr + 4:i_corr + 4 + n_s]],
                corr_e=float(L[i_corr + 1].split("=")[1]), charge_total=float(L[i_corr + 2].split("=")[1]),
                energy=float(L[L.index("#ENERGY") + 1].split()[-1]), volume=hval("Volume"), area=hval("Area"),
                diel=hval("CPCM dielectric energy"), one_el=hval("One-electron operator energy"))


def averaged_sigma(m, rav=0.5):
    """Segment areas (A^2) and screening charge densities sigma (e/A^2) after the Gaussian averaging of
    COSMO-RS over radius rav (A), as in openCOSMO-RS_cpp averageAndClusterSegments."""
    import numpy as np
    pos = np.array([[float(x) for x in s[:3]] for s in m["surf"]]) * BOHR
    area = np.array([float(s[3]) for s in m["surf"]]) * BOHR ** 2
    sig = np.array(m["corr_q"]) / area
    r2 = area / math.pi
    t = r2 + rav * rav
    avg = np.empty(len(area))
    for i in range(len(area)):
        d2 = ((pos - pos[i]) ** 2).sum(1)
        w = (r2 * rav * rav / t) * np.exp(-d2 / t)
        avg[i] = (sig * w).sum() / w.sum()
    return area, avg


def qspr_volume(m, rav=0.5):
    """Molar volume at 25 C (cm3/mol), as in openCOSMO-RS_cpp averageAndClusterSegments."""
    area, avg = averaged_sigma(m, rav)
    m2 = (avg ** 2 * area).sum() * 1e4
    m4 = (avg ** 4 * area).sum() * 1e8
    n_si = sum(1 for l in m["atoms"] if l.split()[-1] == "14")
    return float(0.9430785419976806 * m["n_at"] + 0.6977322963011842 * m["area"] * BOHR ** 2
                 - 0.3161763939689293 * m2 + 0.032441059832647084 * m4 + 8.113026329415828 * n_si
                 - 0.07066832029215675)


def merge(a, b, out, title):
    """Cation a and anion b as one pseudo-molecule, b moved SHIFT_A along x."""
    sb = SHIFT_A / BOHR
    n_at, n_s = a["n_at"] + b["n_at"], len(a["surf"]) + len(b["surf"])
    xyz_b = ["  %s %.10f %s %s" % (l.split()[0], float(l.split()[1]) + SHIFT_A, l.split()[2], l.split()[3])
             for l in b["xyz"]]
    atoms_b = []
    for l in b["atoms"]:
        p = l.split()
        atoms_b.append("  %16.9f %16.9f %16.9f %16.9f %s" % (float(p[0]) + sb, float(p[1]), float(p[2]),
                                                              float(p[3]), p[4]))
    surf_b = []
    for s in b["surf"]:
        s = list(s)
        s[0] = "%.9f" % (float(s[0]) + sb)
        s[9] = str(int(s[9]) + a["n_at"])
        surf_b.append(s)
    L = a["lines"]
    i_c = L.index("#COSMO")
    i_atoms = next(k for k, l in enumerate(L) if "CARTESIAN COORDINATES (A.U.)" in l) - 1
    i_s = next(k for k, l in enumerate(L) if "SURFACE POINTS (A.U.)" in l)
    repl = {"# Number of atoms": "%d\t\t # Number of atoms " % n_at,
            "# Number of surface points": "%d\t\t # Number of surface points" % n_s,
            "# Volume": "   %.9f\t # Volume" % (a["volume"] + b["volume"]),
            "# Area": "   %.9f\t # Area" % (a["area"] + b["area"]),
            "# CPCM dielectric energy": "    %.9f\t # CPCM dielectric energy" % (a["diel"] + b["diel"]),
            "# One-electron operator energy": "    %.9f\t # One-electron operator energy" % (a["one_el"] + b["one_el"])}
    head = []
    for l in L[i_c + 1:i_atoms]:
        head.append(next((v for k, v in repl.items() if k in l), l))
    o = ["%s calculation : DFT_BP86_def2-TZVPD_SP" % title, "", "#" * 50, "#ENERGY",
         "FINAL SINGLE POINT ENERGY  %.12f" % (a["energy"] + b["energy"]), "", "#" * 50, "#XYZ_FILE",
         str(n_at), "merged " + title] + a["xyz"] + xyz_b + ["", "#" * 50, "#COSMO"] + head
    o += L[i_atoms:i_atoms + 3] + a["atoms"] + atoms_b + [""] + L[i_s - 1:i_s + 3]
    o += ["  ".join(s) for s in a["surf"] + surf_b]
    o += ["", "", "#" * 50, "#COSMO_corrected",
          "Corrected dielectric energy   =  %.9f" % (a["corr_e"] + b["corr_e"]),
          "Total C-PCM charge            =  %.9f" % (a["charge_total"] + b["charge_total"]),
          "C-PCM corrected charges:"] + ["  %.9f" % q for q in a["corr_q"] + b["corr_q"]] + [""]
    Path(out).write_text("\n".join(o))


# ---------------------------------------------------------------- systems, inputs, results

# openCOSMO-RS parameters as written by ORCA 6.1.1 for "! COSMORS(...)" (BP86/def2-TZVPD)
PARAMS = {"Aeff": 5.925, "ln_alpha": 0.202, "ln_CHB": 0.166, "CHBT": 1.5, "SigmaHB": 0.009611, "Rav": 0.5,
          "fCorr": 2.4, "RavCorr": 1.0, "comb_SG_A_std": 41.624, "comb_SG_z_coord": 10.0,
          "dGsolv_eta": -4.448, "dGsolv_omega_ring": 0.263,
          "dGsolv_tau": {"1": 0.02934, "6": 0.02288, "7": 0.0007008, "8": 0.003545, "9": 0.005609,
                         "14": 0.004216, "15": 0.003608, "16": 0.03499, "17": 0.03414, "35": 0.04085,
                         "53": 0.213}}


DATASET = ROOT / "data" / "processed" / "model_dataset_curated.csv"
TEMPS = [298.15]
# Only 25 C: the dGsolv model of openCOSMO-RS has fitted, temperature-independent energy terms
# (parametrised at 25 C), so its temperature dependence for CO2 has the wrong sign (tested on a
# 273-373 K grid: K_H of [BMIM][PF6] fell from 20.8 to 5.0 MPa; measured P/x rises from 4.2 to 12.5).


def systems():
    """All 110 cation-anion pairs of the curated dataset plus the 18 Palomar ILs, at 25 C (and at
    the Palomar temperature, 298.0-303.4 K)."""
    man = pd.read_csv(MANIFEST).set_index("ion_abbreviation").file
    out = {}
    for c, a in pd.read_csv(DATASET)[["cation", "anion"]].drop_duplicates().itertuples(index=False):
        out[(man[c], man[a])] = set(TEMPS)
    for r in pd.read_csv(PALOMAR).dropna(subset=["KH_CO2_exp_MPa"]).itertuples():
        out.setdefault((man[r.cation], man[r.anion]), set(TEMPS)).add(float(r.T_K))
    return {k: sorted(v) for k, v in sorted(out.items())}


def cmd_prepare():
    RS.mkdir(parents=True, exist_ok=True)
    co2 = read_cosmo(COSMO / "CO2.orcacosmo")
    (RS / "CO2.orcacosmo").write_text((COSMO / "CO2.orcacosmo").read_text())
    e_gas = float((COSMO / "CO2.egas").read_text())
    rows = []
    for (cat, an), temps in systems().items():
        name = "%s__%s" % (cat, an)
        mc, ma = read_cosmo(COSMO / ("%s.orcacosmo" % cat)), read_cosmo(COSMO / ("%s.orcacosmo" % an))
        merge(mc, ma, RS / ("%s.orcacosmo" % name), name)
        vm = qspr_volume(read_cosmo(RS / ("%s.orcacosmo" % name)))
        j = dict(PARAMS)
        j["dGsolv_numberOfAtomsInRing"] = [0, 0]   # CO2 has no ring; only the solute's value is used
        j["dGsolv_E_gas"] = [e_gas, 0.0]
        j["componentPaths"] = ["CO2.orcacosmo", "%s.orcacosmo" % name]
        j["calculations"] = [{"concentrations": [[0.0, 1.0]] * len(temps), "temperatures": temps,
                              "reference_state_types": [4] * len(temps), "component_indices": [0, 1]}]
        (RS / ("%s.json" % name)).write_text(json.dumps(j, indent=1))
        rows.append(dict(system=name, cation=cat, anion=an, Vm_cm3mol=round(vm, 2), temperatures=" ".join(
            "%.2f" % t for t in temps)))
    pd.DataFrame(rows).to_csv(DEST / "systems.csv", index=False)
    print("%d inputs in %s; CO2 molar volume (QSPR) %.1f cm3/mol" % (len(rows), RS, qspr_volume(co2)))


def henry_table():
    rows = []
    for s in pd.read_csv(DEST / "systems.csv").itertuples():
        f = RS / ("%s_out.json" % s.system)
        if not f.exists():
            continue
        dg = json.loads(f.read_text())["dGsolv"][0]
        for T, row in zip([float(t) for t in s.temperatures.split()], dg):
            RT = R * T
            rows.append(dict(system=s.system, cation=s.cation, anion=s.anion, T_K=T, Vm_cm3mol=s.Vm_cm3mol,
                             dGsolv_kcalmol=row[0],
                             KH_MPa=RT / (s.Vm_cm3mol * 1e-6) * math.exp(row[0] * 4184 / RT) / 1e6))
    return pd.DataFrame(rows)


def cmd_results():
    import numpy as np
    from scipy.stats import pearsonr, spearmanr
    h = henry_table()
    h.round(5).to_csv(DEST / "henry.csv", index=False)
    print("%d systems, %d temperatures" % (h.system.nunique(), len(h)))
    man = pd.read_csv(MANIFEST).set_index("ion_abbreviation").file

    # 1. Palomar et al. 2011: experimental and COSMOtherm Henry constants
    p = pd.read_csv(PALOMAR).dropna(subset=["KH_CO2_exp_MPa"])
    p["system"] = [man[c] + "__" + man[a] for c, a in zip(p.cation, p.anion)]
    hk = h.assign(T_K=h.T_K.round(2))
    d = p.assign(T_K=p.T_K.round(2)).merge(hk[["system", "T_K", "Vm_cm3mol", "dGsolv_kcalmol", "KH_MPa"]],
                                            on=["system", "T_K"])
    d = d.rename(columns={"KH_MPa": "KH_openCOSMO-RS_[C+A]_MPa"})
    d.round(4).to_csv(DEST / "palomar_comparison.csv", index=False)
    print("\nPalomar et al. 2011, %d ionic liquids with experimental Henry constants" % len(d))
    for col in ("KH_openCOSMO-RS_[C+A]_MPa", "KH_COSMO-RS_[C+A]GAS_MPa", "KH_COSMO-RS_[CA]GAS_MPa"):
        lr = np.log(d[col] / d.KH_CO2_exp_MPa)
        print("  %-28s mean ratio to experiment %.2f, SD of ln ratio %.2f, Spearman %.2f"
              % (col, math.exp(lr.mean()), lr.std(), spearmanr(d[col], d.KH_CO2_exp_MPa)[0]))

    # 2. Curated dataset: apparent Henry constant P/x of the points up to 1 MPa at 293-303 K
    from check_validation import TARGET
    x = pd.read_csv(DATASET)
    x = x[(x.Pressure_kPa > 0) & (x.Pressure_kPa <= 1000) & (x[TARGET] > 0)
          & ((x.Temperature_K - 298.15).abs() <= 5)].copy()
    x["system"] = [man[c] + "__" + man[a] for c, a in zip(x.cation, x.anion)]
    x["KH_app_MPa"] = x.Pressure_kPa / 1000 / x[TARGET]
    h298 = h[h.T_K.round(2) == 298.15].set_index("system").KH_MPa
    per = x.groupby("system").agg(n_points=("KH_app_MPa", "size"), KH_app_median=("KH_app_MPa", "median"))
    per["KH_openCOSMO-RS_298K"] = h298.reindex(per.index)
    per["ratio"] = per["KH_openCOSMO-RS_298K"] / per.KH_app_median
    per.round(4).to_csv(DEST / "dataset_comparison.csv")
    lr = np.log(per.ratio)
    print("\nCurated dataset, points up to 1 MPa at 293-303 K: %d points, %d ionic liquids" % (len(x), len(per)))
    print("  per IL (median P/x): mean ratio %.2f, SD of ln ratio %.2f, Spearman %.2f, Pearson(ln) %.2f"
          % (math.exp(lr.mean()), lr.std(), spearmanr(per.KH_app_median, per["KH_openCOSMO-RS_298K"])[0],
             pearsonr(np.log(per.KH_app_median), np.log(per["KH_openCOSMO-RS_298K"]))[0]))
    print("  middle 80%% of the per-IL ratios: %.2f to %.2f" % tuple(np.exp(lr.quantile([0.1, 0.9]))))
    print("  largest deviations (ours / median P over x):")
    for n, r in per.reindex(lr.abs().sort_values(ascending=False).index).head(8).iterrows():
        print("    %-24s %.2f  (%d points)" % (n, r.ratio, r.n_points))

if __name__ == "__main__":
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["cosmo", "prepare", "results"])
    ap.add_argument("ions", nargs="*")
    a = ap.parse_args()
    if a.step == "cosmo":
        cmd_cosmo(a.ions or sorted({"CO2"} | {x for k in systems() for x in k}))
    elif a.step == "prepare":
        cmd_prepare()
    else:
        cmd_results()
