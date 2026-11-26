"""Verify the DOI of every publication in the referenced CO2-solubility dataset.

The citations and titles in data/processed/referenced_co2_solubility_dataset.csv come from
NIST ILThermo, but the DOIs were added by scripts/build_referenced_dataset.py from the first hit of
a Crossref search, without a check. This script checks each DOI against the Crossref record of the
article (title, year, first-author surname) and, when it does not match or is missing, searches
Crossref for the article and accepts a hit only if title, year and first author agree.

Output: data/references/doi_verification.csv, one row per publication:
  status = ok            the stored DOI is the article
           si_suffix     the stored DOI points to a Supporting Information file of the article
           corrected     the stored DOI was wrong or missing; doi_verified is the article found
           manual        set by hand in data/references/doi_overrides.csv (reason in the note column)
           unresolved    no Crossref record agrees with the ILThermo citation (doi_verified empty)
A match needs title similarity >= 0.85, the same year (+-1), the same volume and the same first
author (or a near-identical title). vol_page records whether volume and first page also agree.

Run from the repository root (network access to api.crossref.org needed):
    python scripts/verify_dataset_references.py
"""
import difflib
import json
import re
import ssl
import time
import unicodedata
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

try:  # python.org builds on macOS ship without CA certificates
    import certifi
    SSL_CTX = ssl.create_default_context(cafile=certifi.where())
except ImportError:
    SSL_CTX = ssl.create_default_context()

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "processed" / "referenced_co2_solubility_dataset.csv"
CURATED = ROOT / "data" / "processed" / "curated_co2_solubility_dataset.csv"
OUT = ROOT / "data" / "references" / "doi_verification.csv"
CACHE = ROOT / "data" / "references" / ".crossref_cache.json"
# Checked by hand where the ILThermo record itself is wrong or the work has no Crossref DOI.
OVERRIDES = ROOT / "data" / "references" / "doi_overrides.csv"
# DOIs as first written by scripts/build_referenced_dataset.py (unchecked Crossref search), kept so
# that the check still compares against them after apply_verified_references.py has fixed the data.
FROM_SEARCH = ROOT / "data" / "references" / "doi_from_search.csv"
UA = "il-co2-prediction reference check (python urllib)"


def norm(s):
    s = unicodedata.normalize("NFKD", str(s)).encode("ascii", "ignore").decode()
    s = re.sub(r"<[^>]+>", "", s)                      # Crossref titles may carry <sub> tags
    return re.sub(r"[^a-z0-9]+", " ", s.lower()).strip()


def similar(a, b):
    return difflib.SequenceMatcher(None, norm(a), norm(b)).ratio()


def first_surname(citation):
    # ILThermo style: "Jung, Y.-H.; Jung, J.-Y.; ... (2012) J. Chem. Eng. Data ..."
    return norm(citation.split(",")[0]).split(" ")[-1] if citation else ""


class Crossref:
    def __init__(self):
        self.cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}

    def get(self, url):
        if url not in self.cache:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            for attempt in range(3):
                try:
                    with urllib.request.urlopen(req, timeout=30, context=SSL_CTX) as r:
                        self.cache[url] = json.load(r)["message"]
                    break
                except urllib.error.HTTPError as e:
                    if e.code == 404:
                        self.cache[url] = None
                        break
                    time.sleep(3 * (attempt + 1))
                except Exception:
                    time.sleep(3 * (attempt + 1))
            else:
                return None
            time.sleep(0.3)
        return self.cache[url]

    def work(self, doi):
        return self.get("https://api.crossref.org/works/" + urllib.parse.quote(doi))

    def search(self, title, citation):
        q = urllib.parse.quote(f"{title} {citation}")
        msg = self.get(f"https://api.crossref.org/works?query.bibliographic={q}&rows=8")
        return (msg or {}).get("items", [])

    def save(self):
        CACHE.write_text(json.dumps(self.cache))


def describe(item):
    if not item:
        return {}
    year = (item.get("issued", {}).get("date-parts") or [[None]])[0][0]
    authors = item.get("author") or [{}]
    return {"cr_title": (item.get("title") or [""])[0], "cr_year": year,
            "cr_first_author": authors[0].get("family", ""),
            "cr_journal": (item.get("container-title") or [""])[0], "cr_type": item.get("type", "")}


def volume_page(citation):
    # "... (2012) J. Chem. Eng. Data 57(12), 3321-3329." -> ("57", "3321")
    m = re.search(r"\)\s.*?\b(\d+)(?:\([^)]*\))?,\s*[A-Za-z]?(\d+)", str(citation))
    return (m.group(1), m.group(2)) if m else (None, None)


def digits(x):
    return re.sub(r"\D", "", str(x or "")).lstrip("0")


def vol_page_check(item, citation):
    """'ok', 'page differs' or 'volume differs' between the ILThermo citation and Crossref."""
    vol, page = volume_page(citation)
    cr_vol, cr_page = item.get("volume"), (item.get("page") or item.get("article-number") or "").split("-")[0]
    if vol and cr_vol and digits(vol) not in [digits(x) for x in re.findall(r"\d+", str(cr_vol))]:
        return "volume differs"
    if page and cr_page and digits(page) != digits(cr_page):
        return "page differs"
    return "ok"


def agrees(item, title, year, surname, citation):
    d = describe(item)
    if not d:
        return False, 0.0
    sim = similar(d["cr_title"], title)
    year_ok = d["cr_year"] is not None and abs(int(d["cr_year"]) - int(year)) <= 1
    author_ok = bool(surname) and norm(d["cr_first_author"]).split(" ")[-1:] == [surname]
    vol_ok = vol_page_check(item, citation) != "volume differs"
    if d["cr_type"] == "posted-content":  # a preprint (e.g. SSRN), not the published article
        return False, sim
    return (sim >= 0.85 and year_ok and vol_ok and (author_ok or sim >= 0.97)), sim


def main():
    df = pd.read_csv(DATA)
    if CURATED.exists():  # publications that only the curated (strict ILThermo) build contains
        cur = pd.read_csv(CURATED).rename(columns={"set_id": "entry_id"})
        extra = cur[~cur.ref_citation.isin(df.ref_citation)].assign(ref_doi=None)
        df = pd.concat([df, extra[df.columns.intersection(extra.columns)]], ignore_index=True)
    pubs = (df.groupby(["ref_citation", "ref_title", "ref_year"], dropna=False)
              .agg(doi_stored=("ref_doi", "first"), n_points=("entry_id", "size"),
                   ils=("cation", lambda s: len(set(zip(s, df.loc[s.index, "anion"])))))
              .reset_index())
    overrides = {r.ref_citation: r for r in pd.read_csv(OVERRIDES).itertuples()} if OVERRIDES.exists() else {}
    if FROM_SEARCH.exists():
        searched = pd.read_csv(FROM_SEARCH).fillna("")
        pubs["doi_stored"] = pubs.ref_citation.map(dict(zip(searched.ref_citation, searched.doi_from_search)))
    cr = Crossref()
    rows = []
    for i, p in pubs.iterrows():
        title, cit, year = p.ref_title, p.ref_citation, p.ref_year
        surname = first_surname(cit)
        stored = p.doi_stored if isinstance(p.doi_stored, str) and p.doi_stored.strip() else ""
        base = re.sub(r"\.s\d{3}$", "", stored)
        status, verified, item, sim = "unresolved", "", None, 0.0
        if cit in overrides:
            o = overrides[cit]
            item = cr.work(o.doi) if isinstance(o.doi, str) and o.doi else None
            rows.append({"ref_citation": cit, "ref_title": title, "ref_year": year, "n_points": p.n_points,
                         "doi_stored": stored, "doi_verified": o.doi if isinstance(o.doi, str) else "",
                         "status": "manual", "title_similarity": round(similar(describe(item).get("cr_title", ""), title), 3) if item else None,
                         "vol_page": vol_page_check(item, cit) if item else "", "note": o.reason,
                         **describe(item)})
            continue
        if base:
            item = cr.work(base)
            ok, sim = agrees(item, title, year, surname, cit)
            if ok:
                status, verified = ("si_suffix" if base != stored else "ok"), item["DOI"].lower()
        if not verified:
            best = None
            for cand in cr.search(title, cit):
                if cand.get("type") not in ("journal-article", "proceedings-article", "book-chapter"):
                    continue
                ok, s = agrees(cand, title, year, surname, cit)
                if ok and (best is None or s > best[1]):
                    best = (cand, s)
            if best:
                item, sim = best
                verified = item["DOI"].lower()
                status = "ok" if verified == stored.lower() else ("si_suffix" if verified == base.lower() else "corrected")
            else:
                item = cr.work(base) if base else None
        rows.append({"ref_citation": cit, "ref_title": title, "ref_year": year, "n_points": p.n_points,
                     "doi_stored": stored, "doi_verified": verified, "status": status,
                     "title_similarity": round(sim, 3), "vol_page": vol_page_check(item, cit) if item else "",
                     "note": "", **describe(item)})
        if i % 20 == 0:
            cr.save()
            print(f"{i}/{len(pubs)}")
    cr.save()
    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(out.status.value_counts().to_string())
    print("points by status:", out.groupby("status").n_points.sum().to_dict())
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
