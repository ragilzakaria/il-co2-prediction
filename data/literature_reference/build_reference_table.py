"""Combine the literature reference values in raw/ into one long table, and report coverage.

Inputs (raw/, one file per source, extracted from the papers listed in README.md):
  karu2018_pairs.csv, il2013_identified.csv, xtbdb_pairs.csv, xtbdb_ions.csv,
  hunt2025_ion_IE.csv, palomar2011_cosmors.csv, doherty2017_md.csv, kerle2017_md_henry.csv

Outputs:
  reference_values.csv  one row per (system, quantity, method, source, conditions)
  coverage_by_pair.csv  for every cation-anion pair in the training data: which quantities have references
  coverage_by_ion.csv   the same for single ions

Usage:  python data/literature_reference/build_reference_table.py
"""
import collections
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw"
DATA = HERE.parent

SRC = {
    "karu2018": ("Karu et al. 2018, Int. J. Quantum Chem. 118, e25582", "10.1002/qua.25582",
                 "data: github.com/doublelayer/SCAN (Zenodo 10.5281/zenodo.292628)"),
    "il2013": ("Zahn, MacFarlane, Izgorodina 2013, Phys. Chem. Chem. Phys. 15, 13664 (IL-2013 set)",
               "10.1039/C3CP51682B", "ESI file IL-2013-set.xyz; ion identities assigned by graph matching"),
    "xtbdb": ("Stavarache, Calero, Vicent-Luna 2026, arXiv:2610.15190 (ionic-liquids-database v1.0)",
              "10.5281/zenodo.20539350", "CC0; lowest-energy of 10 DFT conformers per pair"),
    "dos2025": ("Phys. Chem. Chem. Phys. 2025, 27, 9068 (IL density of states from lone-ion calculations)",
                "10.1039/D5CP00214A", "ESI Tables S12-S13"),
    "palomar2011": ("Palomar et al. 2011, Ind. Eng. Chem. Res. 50, 3452", "10.1021/ie101572m",
                    "SI Tables S3-S4"),
    "doherty2017": ("Doherty et al. 2017, J. Chem. Theory Comput. 13, 6131", "10.1021/acs.jctc.7b00520",
                    "Tables 2 and 6 (open access)"),
    "kerle2017": ("Kerlé et al. 2017, Phys. Chem. Chem. Phys. 19, 1770", "10.1039/C6CP06792A",
                  "Table 2 (open access)"),
}
COLS = ["cation", "anion", "ion", "category", "quantity", "value", "uncertainty", "units", "method",
        "conditions", "experiment_same_conditions", "source", "doi", "location"]


def rows():
    k, (s, d, loc) = "karu2018", SRC["karu2018"]
    for r in csv.DictReader(open(RAW / "karu2018_pairs.csv")):
        if not r["cation"].startswith("["):
            continue  # triethylpropylammonium is not in our set
        base = dict(cation=r["cation"], anion=r["anion"], category="QM ion pair", source=s, doi=d, location=loc)
        yield {**base, "quantity": "interaction energy", "value": r["E_int_DLPNO-CCSD(T)/CBS_kJmol"], "units": "kJ/mol",
               "method": "DLPNO-CCSD(T)/CBS (TightPNO) on B2PLYP-D3/def2-TZVPP geometry", "conditions": "gas phase, one minimum"}
        yield {**base, "quantity": "interaction energy", "value": r["E_int_B2PLYP-D3+CP_kJmol"], "units": "kJ/mol",
               "method": "B2PLYP-D3(BJ)/def2-TZVPP + counterpoise", "conditions": "gas phase, one minimum"}
        yield {**base, "quantity": "ion-pair dipole moment", "value": r["dipole_REF_D"], "units": "D",
               "method": "reference level (see paper)", "conditions": "gas phase, one minimum"}

    s, d, loc = SRC["il2013"]
    by = collections.defaultdict(list)
    for r in csv.DictReader(open(RAW / "il2013_identified.csv")):
        ions = r["ions"].split(" | ")
        if any(i.startswith("?") for i in ions):
            continue
        by[tuple(ions)].append(float(r["E_int_CCSDT_CBS_kJmol"]))
    for ions, v in by.items():
        cat = next(i for i in ions if i in CATIONS)
        an = next(i for i in ions if i != cat)
        yield dict(cation=cat, anion=an, category="QM ion pair", quantity="interaction energy",
                   value=f"{min(v):.1f}", units="kJ/mol",
                   method="CCSD(T)/CBS, counterpoise corrected (MP2/CBS + CCSD(T) small-basis correction)",
                   conditions=f"gas phase; most stable of {len(v)} structures (range {min(v):.1f} to {max(v):.1f})",
                   source=s, doi=d, location=loc)

    s, d, loc = SRC["xtbdb"]
    for r in csv.DictReader(open(RAW / "xtbdb_pairs.csv")):
        base = dict(cation=r["cation"], anion=r["anion"], category="QM ion pair", source=s, doi=d, location=loc,
                    conditions=f"gas phase; lowest of {r['n_dft_conformers']} conformers")
        yield {**base, "quantity": "binding energy (vs relaxed ions)", "value": r["E_bind_min_PBE-D3BJ_DZP_noCP_kJmol"],
               "units": "kJ/mol", "method": "PBE-D3(BJ)/DZP (ADF, frozen core), no counterpoise"}
        yield {**base, "quantity": "HOMO-LUMO gap", "value": r["HOMO_LUMO_gap_lowestE_eV"], "units": "eV",
               "method": "PBE-D3(BJ)/DZP (ADF)"}
    for r in csv.DictReader(open(RAW / "xtbdb_ions.csv")):
        yield dict(ion=r["ion"], category="QM single ion", quantity="HOMO-LUMO gap",
                   value=r["HOMO_LUMO_gap_PBE-D3BJ_DZP_eV"], units="eV", method="PBE-D3(BJ)/DZP (ADF)",
                   conditions="gas phase, one conformer", source=s, doi=d, location=loc)

    s, d, loc = SRC["dos2025"]
    for r in csv.DictReader(open(RAW / "hunt2025_ion_IE.csv")):
        yield dict(ion=r["ion"], category="QM single ion", quantity="ionisation energy in IL (≈ -E_HOMO)",
                   value=r["Ei_IL_lone-ion-SMD_eV"], units="eV", method="B3LYP-D3(BJ)/6-311+G(d,p), SMD continuum",
                   conditions="lone ion in IL-like continuum, aligned to valence XPS",
                   source=s, doi=d, location=r["location"])

    s, d, _ = SRC["palomar2011"]
    for r in csv.DictReader(open(RAW / "palomar2011_cosmors.csv")):
        for col, model in (("KH_COSMO-RS_[CA]GAS_MPa", "ion-pair [CA], gas-phase geometry (model used in SI Tables S4-S5)"),
                           ("KH_COSMO-RS_[C+A]GAS_MPa", "separate ions [C+A], gas-phase geometry")):
            if r[col]:
                yield dict(cation=r["cation"], anion=r["anion"], category="COSMO-RS",
                           quantity="CO2 Henry constant", value=r[col], units="MPa",
                           method=f"COSMO-RS, {model}", conditions=f"T = {r['T_K']} K",
                           experiment_same_conditions=r["KH_CO2_exp_MPa"], source=s, doi=d, location=r["location"])

    s, d, _ = SRC["doherty2017"]
    for r in csv.DictReader(open(RAW / "doherty2017_md.csv")):
        for col, ff in (("sim_OPLS-2009IL", "OPLS-2009IL (full charges)"), ("sim_0.8xOPLS-2009IL", "OPLS-2009IL, charges x0.8")):
            yield dict(cation=r["cation"], anion=r["anion"], category="MD bulk liquid", quantity=r["property"],
                       value=r[col], units=r["units"], method=f"MD, {ff}", conditions="298 K",
                       experiment_same_conditions=r["experiment"], source=s, doi=d, location=r["location"])

    s, d, loc = SRC["kerle2017"]
    for r in csv.DictReader(open(RAW / "kerle2017_md_henry.csv")):
        yield dict(cation=r["cation"], anion=r["anion"], category="MD bulk liquid", quantity="CO2 Henry constant",
                   value=r["H_CO2_sim_bar"], uncertainty=r["uncertainty_bar"], units="bar",
                   method="MD, non-polarisable force field; solvation free energy (Widom/BAR)",
                   conditions=f"T = {r['T_K']} K, infinite dilution", source=s, doi=d, location=loc)


IONS = list(csv.DictReader(open(DATA / "ions" / "ions_with_smiles.csv")))
CATIONS = {r["ion_abbreviation"] for r in IONS if r["ion_type"] == "cation"}


def main():
    out = [{c: r.get(c, "") for c in COLS} for r in rows()]
    with open(HERE / "reference_values.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        w.writerows(out)

    pairs = collections.Counter()
    with open(DATA / "processed" / "final_referenced_model_dataset.csv") as f:
        for r in csv.DictReader(f):
            pairs[(r["cation"], r["anion"])] += 1
    have = collections.defaultdict(set)
    ion_have = collections.defaultdict(set)
    for r in out:
        key = f"{r['category']}: {r['quantity']}"
        if r["cation"]:
            have[(r["cation"], r["anion"])].add(key)
            ion_have[r["cation"]].add("in a referenced pair")
            ion_have[r["anion"]].add("in a referenced pair")
        else:
            ion_have[r["ion"]].add(key)

    with open(HERE / "coverage_by_pair.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cation", "anion", "data_points", "n_reference_types", "reference_types"])
        for (c, a), n in pairs.most_common():
            k = sorted(have.get((c, a), ()))
            w.writerow([c, a, n, len(k), "; ".join(k)])
    with open(HERE / "coverage_by_ion.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["ion", "type", "in_training_data", "reference_types"])
        used = {c for c, _ in pairs} | {a for _, a in pairs}
        for r in IONS:
            i = r["ion_abbreviation"]
            w.writerow([i, r["ion_type"], "yes" if i in used else "no", "; ".join(sorted(ion_have.get(i, ())))])

    covered = [p for p in pairs if have.get(p)]
    pts = sum(pairs[p] for p in covered)
    print(f"{len(out)} reference values from {len(SRC)} sources")
    print(f"training pairs with at least one reference: {len(covered)}/{len(pairs)} "
          f"({pts}/{sum(pairs.values())} data points)")
    for cat in ("QM ion pair", "COSMO-RS", "MD bulk liquid"):
        n = sum(1 for p in pairs if any(k.startswith(cat) for k in have.get(p, ())))
        print(f"  {cat}: {n} pairs")
    single = {r["ion"] for r in out if r["ion"]}
    print(f"single ions with orbital/IE references: {len(single & used)}/{len(used)} training ions")


if __name__ == "__main__":
    main()
