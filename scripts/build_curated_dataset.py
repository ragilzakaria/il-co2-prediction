"""Build the curated CO2-solubility dataset directly from the raw NIST ILThermo records.

Replaces the loose parser of scripts/build_referenced_dataset.py. Input: the ILThermo data sets
saved as received in data/raw/ilthermo/<set>.json (one per row of data/raw/matched_ilthermo_targets.csv).
Every data set and every point is either kept or excluded with a reason.

Rules
  Set level
    - binary system CO2 + the IL named in the target table (all 462 sets are binary)
    - the IL is the one we computed (acetylcholine [Tf2N] was matched to choline [Tf2N]: excluded)
    - temperature, pressure (total or CO2 partial) and a CO2 composition of the IL-rich liquid
      phase must all be columns of the set
    - composition given as mole fraction, molality, mole ratio, mass ratio or weight fraction;
      molarity (needs the solution density) is excluded
    - a composition of the IL in the CO2 phase, or of the CO2-rich phase, is excluded
    ILThermo's phase labels are not always right: a CO2 composition labelled "Gas" whose values
    are well below 1 can only be the liquid (the vapour over an IL is almost pure CO2); it is kept
    and noted. "Liquid mixture 1" is the IL-rich liquid of a two-liquid system and is kept.
  Point level
    - x_CO2 >= 0.95: CO2-rich phase, not CO2 dissolved in the IL
    - repeated row within a set: dropped
    - the same (IL, T, P, x) reported by two publications: kept once, from the earlier one
  Hand-checked
    - data sets listed in data/curation/excluded_sets.csv (found by scripts/check_data_consistency.py
      to disagree with the consensus of other publications; the evidence is in the file)
Molar masses of the ions (for molality and mass-based units) come from the optimised ORCA
structures in calculations/ions/ (the corrected species).

Output
  data/processed/curated_co2_solubility_dataset.csv   kept points, with set id, citation, verified DOI,
                                                       method, original unit and uncertainty
  data/processed/curated_exclusions.csv               excluded points and sets, with the reason

Run from the repository root:
    python scripts/build_curated_dataset.py
"""
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "ilthermo"
TARGETS = ROOT / "data" / "raw" / "matched_ilthermo_targets.csv"
IONS = ROOT / "data" / "ions" / "ions_with_smiles.csv"
ION_OUT = ROOT / "calculations" / "ions"
DOIS = ROOT / "data" / "references" / "doi_verification.csv"
OUT = ROOT / "data" / "processed" / "curated_co2_solubility_dataset.csv"
EXCL = ROOT / "data" / "processed" / "curated_exclusions.csv"
# Data sets removed by hand after scripts/check_data_consistency.py, with the evidence
MANUAL = ROOT / "data" / "curation" / "excluded_sets.csv"

M_CO2 = 44.0095
X_MAX = 0.95
WRONG_IL = {"acetylcholine bis(trifluoromethylsulfonyl)imide":
            "acetylcholine [Tf2N], matched by mistake to choline [Tf2N]"}
ATOMIC_MASS = {"H": 1.008, "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998,
               "P": 30.974, "S": 32.06, "Cl": 35.45}


def safe_name(abbr):
    return (abbr.strip("[]").replace("+", "_p").replace("-", "_m")
            .replace("(", "_").replace(")", "_"))


def ion_molar_masses():
    masses = {"[Cl]": 35.45}
    for abbr in pd.read_csv(IONS).ion_abbreviation.unique():
        f = ION_OUT / f"{safe_name(abbr)}.out"
        if not f.exists():
            continue
        block = f.read_text(errors="ignore").split("CARTESIAN COORDINATES (ANGSTROEM)")[-1].split("\n\n")[0]
        atoms = [l.split()[0] for l in block.splitlines()[2:] if len(l.split()) == 4]
        masses[abbr] = sum(ATOMIC_MASS[a] for a in atoms)
    return masses


def cell(c):
    v = float(c[0])
    u = float(c[1]) if len(c) > 1 and c[1] not in (None, "") else None
    return v, u


def classify(dhead):
    """Return (t_idx, p_idx, p_kind, x_idx, x_kind, phase) or an exclusion reason (str)."""
    t_idx = p_idx = None
    p_kind = None
    comps = []
    for i, h in enumerate(dhead):
        name, phase = h[0], (h[1] if len(h) > 1 else None)
        low = name.lower()
        if low.startswith("temperature, k"):
            t_idx = i
        elif "pressure" in low and "kpa" in low:
            p_idx, p_kind = i, ("CO2 partial pressure" if "partial" in low else "pressure")
        elif "carbon dioxide" in low:
            comps.append((i, low, phase or ""))
        elif any(k in low for k in ("mole fraction of", "fraction of")):
            return "composition of the IL in the CO2 phase, not CO2 in the IL"
    if t_idx is None or p_idx is None:
        return "no temperature or no pressure column"
    if not comps:
        return "no CO2 composition column"
    # prefer the liquid phase, then an unlabelled/'Gas' column (checked on values later)
    comps.sort(key=lambda c: (0 if c[2].startswith("Liquid") and "mixture 2" not in c[2] else
                              1 if "mixture 2" not in c[2] else 2))
    i, low, phase = comps[0]
    if "mixture 2" in phase:
        return "CO2-rich liquid phase (Liquid mixture 2)"
    if "molarity" in low:
        return "molarity (mol/dm3); conversion needs the solution density"
    kind = ("mole fraction" if "mole fraction" in low else "molality" if "molality" in low else
            "mole ratio" if "mole ratio" in low else "mass ratio" if "mass ratio" in low else
            "weight fraction" if "weight fraction" in low else None)
    if kind is None:
        return f"unknown composition unit: {low}"
    return t_idx, p_idx, p_kind, i, kind, phase


def to_mole_fraction(v, kind, m_il):
    if kind == "mole fraction":
        return v
    if kind == "molality":            # mol CO2 per kg IL
        return v / (v + 1000.0 / m_il)
    if kind == "mole ratio":          # mol CO2 per mol IL
        return v / (1.0 + v)
    if kind == "mass ratio":          # g CO2 per g IL
        n = v / M_CO2
        return n / (n + 1.0 / m_il)
    if kind == "weight fraction":     # g CO2 per g solution
        n = v / M_CO2
        return n / (n + (1.0 - v) / m_il)
    raise ValueError(kind)


def main():
    targets = pd.read_csv(TARGETS)
    masses = ion_molar_masses()
    ver = pd.read_csv(DOIS)
    doi_of = dict(zip(ver.ref_citation, ver.doi_verified))
    manual = pd.read_csv(MANUAL).set_index("set_id").reason.to_dict() if MANUAL.exists() else {}
    kept, excluded = [], []
    for t in targets.itertuples():
        d = json.load(open(RAW / f"{t.id}.json"))
        comps = [c if isinstance(c, str) else c.get("name") for c in d["components"]]
        il_name = next((c for c in comps if c != "carbon dioxide"), "")
        cit, title = d["ref"]["full"], d["ref"].get("title", "")
        base = {"set_id": t.id, "cation": t.cation, "anion": t.anion, "ref_citation": cit}

        def drop_set(reason):
            excluded.append({**base, "level": "set", "n_points": len(d["data"]), "reason": reason})

        if len(comps) != 2 or "carbon dioxide" not in comps:
            drop_set("not a binary CO2 + IL system"); continue
        if il_name in WRONG_IL:
            drop_set(WRONG_IL[il_name]); continue
        if t.id in manual:
            drop_set("inconsistent with other publications: " + manual[t.id]); continue
        c = classify(d["dhead"])
        if isinstance(c, str):
            drop_set(c); continue
        t_idx, p_idx, p_kind, x_idx, kind, phase = c
        m_il = masses.get(t.cation, 0) + masses.get(t.anion, 0) if t.cation in masses and t.anion in masses else None
        if kind != "mole fraction" and m_il is None:
            drop_set("no molar mass for unit conversion"); continue
        raw_x = [float(row[x_idx][0]) for row in d["data"]]
        note = ""
        if phase == "Gas":
            if pd.Series(raw_x).median() > X_MAX:
                drop_set("CO2 composition of the gas phase"); continue
            note = "labelled Gas in ILThermo; values are liquid-phase solubility"
        for row in d["data"]:
            T, _ = cell(row[t_idx])
            P, _ = cell(row[p_idx])
            v, u = cell(row[x_idx])
            x = to_mole_fraction(v, kind, m_il)
            ux = None if u is None else (u if kind == "mole fraction" else u * x / v if v else None)
            rec = {**base, "Temperature_K": round(T, 2), "Pressure_kPa": round(P, 3),
                   "CO2_solubility": round(x, 6), "solubility_uncertainty": ux,
                   "pressure_kind": p_kind, "original_unit": kind, "phase_label": phase,
                   "note": note, "exp_method": d.get("expmeth") or "unspecified", "ref_title": title}
            if not (0 < x < X_MAX):
                excluded.append({**base, "level": "point", "n_points": 1, "Temperature_K": T,
                                 "Pressure_kPa": P, "CO2_solubility": x,
                                 "reason": "x_CO2 >= 0.95: CO2-rich phase" if x >= X_MAX else "x_CO2 <= 0"})
                continue
            kept.append(rec)
    df = pd.DataFrame(kept)
    key = ["cation", "anion", "Temperature_K", "Pressure_kPa", "CO2_solubility"]
    # repeated rows within a set
    dup_in = df.duplicated(key + ["set_id"])
    for r in df[dup_in].itertuples():
        excluded.append({"set_id": r.set_id, "cation": r.cation, "anion": r.anion, "ref_citation": r.ref_citation,
                         "level": "point", "n_points": 1, "Temperature_K": r.Temperature_K,
                         "Pressure_kPa": r.Pressure_kPa, "CO2_solubility": r.CO2_solubility,
                         "reason": "repeated row within the data set"})
    df = df[~dup_in]
    # same values in two publications: keep the earlier one
    df["ref_year"] = df.ref_citation.str.extract(r"\((\d{4})\)")[0].astype(int)
    df = df.sort_values(["ref_year", "set_id"])
    dup_x = df.duplicated(key)
    first = df.drop_duplicates(key).set_index(key).ref_citation
    for r in df[dup_x].itertuples():
        excluded.append({"set_id": r.set_id, "cation": r.cation, "anion": r.anion, "ref_citation": r.ref_citation,
                         "level": "point", "n_points": 1, "Temperature_K": r.Temperature_K,
                         "Pressure_kPa": r.Pressure_kPa, "CO2_solubility": r.CO2_solubility,
                         "reason": "same values already reported in " + first.loc[tuple(getattr(r, k) for k in key)][:60]})
    df = df[~dup_x]
    df["ref_doi"] = df.ref_citation.map(doi_of)
    missing_doi = df.ref_doi.isna().sum()
    df = df.sort_values(["cation", "anion", "ref_year", "set_id", "Temperature_K", "Pressure_kPa"])
    cols = ["set_id", "cation", "anion", "Temperature_K", "Pressure_kPa", "CO2_solubility",
            "solubility_uncertainty", "pressure_kind", "original_unit", "phase_label", "note",
            "exp_method", "ref_year", "ref_citation", "ref_title", "ref_doi"]
    df[cols].to_csv(OUT, index=False)
    ex = pd.DataFrame(excluded)
    ex.to_csv(EXCL, index=False)
    print(f"kept {len(df)} points, {df.groupby(['cation', 'anion']).ngroups} ILs, "
          f"{df.ref_citation.nunique()} publications, {df.set_id.nunique()} data sets; "
          f"points without verified DOI: {missing_doi}")
    print("excluded points by reason:")
    print(ex.groupby("reason").n_points.sum().sort_values(ascending=False).to_string())
    print("original units of kept points:", df.original_unit.value_counts().to_dict())


if __name__ == "__main__":
    main()
