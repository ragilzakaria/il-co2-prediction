"""Apply the verified DOIs to the referenced dataset and rebuild its bibliography.

Reads data/references/doi_verification.csv (from scripts/verify_dataset_references.py) and
  1. replaces ref_doi in data/processed/referenced_co2_solubility_dataset.csv by the verified DOI
     (only that column changes);
  2. rewrites data/references/dataset_bibliography.csv (one row per publication, with the DOI
     check status) and data/references/dataset_references.bib (BibTeX, authors joined by "and",
     journal, volume and pages from the Crossref record where available).

Afterwards run scripts/merge_referenced_model_dataset.py (with and without --thesis) so the model
tables carry the same DOIs.

Run from the repository root:
    python scripts/apply_verified_references.py
"""
import json
import re
import unicodedata
import urllib.parse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "processed" / "referenced_co2_solubility_dataset.csv"
REFS = ROOT / "data" / "references"
VERIFY = REFS / "doi_verification.csv"
CACHE = REFS / ".crossref_cache.json"


def crossref_record(doi, cache):
    return cache.get("https://api.crossref.org/works/" + urllib.parse.quote(doi)) if doi else None


def bib_escape(s):
    return str(s).replace("&", r"\&").replace("%", r"\%").replace("_", r"\_")


def strip_tags(s):
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", str(s))).strip()


def key_for(citation, year, used):
    first = unicodedata.normalize("NFKD", citation.split(",")[0]).encode("ascii", "ignore").decode()
    base = re.sub(r"[^a-z]", "", first.lower()) + str(year)
    key, n = base, 1
    while key in used:
        n += 1
        key = f"{base}{chr(ord('a') + n - 1)}"
    used.add(key)
    return key


def main():
    v = pd.read_csv(VERIFY)
    missing = v[v.doi_verified.isna() | (v.doi_verified.astype(str).str.strip() == "")]
    if len(missing):
        raise SystemExit(f"{len(missing)} publications have no verified DOI; resolve them first:\n"
                         + "\n".join(missing.ref_citation))
    doi_of = dict(zip(v.ref_citation, v.doi_verified.str.lower()))

    df = pd.read_csv(DATA)
    before = df.ref_doi.copy()
    df["ref_doi"] = df.ref_citation.map(doi_of)
    assert df.ref_doi.notna().all()
    df.to_csv(DATA, index=False)
    changed = (before.fillna("") != df.ref_doi).sum()
    print(f"ref_doi changed for {changed} of {len(df)} points; {df.ref_doi.nunique()} distinct DOIs")

    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    pubs = (df.assign(il=df.cation + df.anion)
              .groupby(["ref_citation", "ref_title", "ref_year", "ref_doi"])
              .agg(data_points=("entry_id", "size"),
                   ionic_liquids=("il", lambda s: ", ".join(sorted(set(s)))),
                   num_unique_ils=("il", "nunique"),
                   exp_methods=("exp_method", lambda s: "; ".join(sorted(set(s.astype(str))))))
              .reset_index()
              .merge(v[["ref_citation", "status", "note"]].rename(columns={"status": "doi_check"}),
                     on="ref_citation", how="left")
              .sort_values(["data_points"], ascending=False))
    pubs.to_csv(REFS / "dataset_bibliography.csv", index=False)

    used, entries = set(), []
    for p in pubs.itertuples():
        rec = crossref_record(p.ref_doi, cache) or {}
        authors = [f"{a.get('family', '')}, {a.get('given', '')}".strip(", ") for a in rec.get("author", [])]
        if thesis_author := {"Sanchez, L.M. Galan": "Gal{\\'a}n S{\\'a}nchez, L. M."}.get(p.ref_citation.split(" (")[0]):
            authors = [thesis_author]  # ILThermo inverts the Spanish double surname
        if not authors:  # fall back to the ILThermo author list "Surname, I.; Surname, I. (year) ..."
            authors = [a.strip() for a in p.ref_citation.split("(")[0].split(";") if a.strip()]
        thesis = "Thesis" in p.ref_citation
        fields = {
            "author": " and ".join(authors),
            "title": "{" + strip_tags((rec.get("title") or [p.ref_title])[0]) + "}",
            "year": str(p.ref_year),
            "doi": p.ref_doi,
        }
        if thesis:
            fields["school"] = "Eindhoven University of Technology"
        else:
            fields["journal"] = (rec.get("container-title") or [re.sub(r"^.*?\)\s*", "", p.ref_citation).split(" ")[0]])[0]
            for k_cr, k_bib in (("volume", "volume"), ("issue", "number"), ("page", "pages")):
                if rec.get(k_cr):
                    fields[k_bib] = str(rec[k_cr]).replace("-", "--")
            if not rec.get("page") and rec.get("article-number"):
                fields["pages"] = rec["article-number"]
        fields["note"] = (f"CO2 solubility: {p.data_points} points for {p.ionic_liquids} "
                          f"(source record: NIST ILThermo; DOI check: {p.doi_check})")
        body = ",\n".join(f"  {k:<7} = {{{bib_escape(val) if k not in ('doi', 'title') else val}}}"
                          for k, val in fields.items())
        entries.append(f"@{'phdthesis' if thesis else 'article'}{{{key_for(p.ref_citation, p.ref_year, used)},\n{body}\n}}")
    (REFS / "dataset_references.bib").write_text(
        "% Sources of the referenced CO2-solubility dataset (data/processed/referenced_co2_solubility_dataset.csv).\n"
        "% Generated by scripts/apply_verified_references.py from data/references/doi_verification.csv.\n\n"
        + "\n\n".join(entries) + "\n")
    print(f"wrote {len(pubs)} publications to dataset_bibliography.csv and dataset_references.bib")


if __name__ == "__main__":
    main()
