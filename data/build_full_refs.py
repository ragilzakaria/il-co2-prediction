#!/usr/bin/env python3
"""
build_full_refs.py — Merge thesis references, dataset references, and foundational citations
into workspace/refs.bib and workspace/citation_pool.json.
"""
import re
import json
import os

# 1. Read dataset references
with open("il-co2-prediction/data/dataset_references.bib") as f:
    dataset_bib = f.read()

# 2. Read thesis references
with open("il-co2-prediction/data/thesis_references_extracted.txt") as f:
    thesis_lines = f.readlines()

cleaned = []
for line in thesis_lines:
    stripped = line.strip()
    if re.match(r'^\d+$', stripped) or stripped == 'REFERENCES':
        continue
    if stripped.startswith('Attachment 1') or stripped.startswith('Table 1.') or stripped.startswith('No Cation'):
        break
    cleaned.append(line)

text = ''.join(cleaned)
raw_entries = re.split(r'\n(?=[A-Z][a-zA-Z\s\.,\-]+,\s+[A-Z]\.)', text)

thesis_bib_entries = []
seen_keys = set(re.findall(r"^@\w+\{\s*([^,\s]+)", dataset_bib, re.M))

for entry in raw_entries:
    e = ' '.join(entry.split())
    if not e:
        continue
    ym = re.search(r'\((\d{4})\)', e)
    if not ym:
        continue
    year = ym.group(1)
    author_part = e[:ym.start()].strip()
    rest = e[ym.end():].strip().lstrip('.').strip()
    
    dm = re.search(r'https?://doi\.org/([^\s]+)', rest)
    doi = dm.group(1) if dm else ''
    if dm:
        rest = rest[:dm.start()].strip()
        
    parts = rest.split('. ')
    title = parts[0].strip().rstrip('.')
    journal = '. '.join(parts[1:]).strip() if len(parts) > 1 else 'Thesis Literature'
    
    first_author = author_part.split(',')[0].strip().lower()
    first_author = re.sub(r'[^a-z]', '', first_author)
    first_word = re.sub(r'[^a-z]', '', title.split()[0].lower()) if title.split() else 'ref'
    base_key = f"{first_author}_{year}_{first_word}"
    key = base_key
    counter = 1
    while key in seen_keys:
        counter += 1
        key = f"{base_key}_{counter}"
    seen_keys.add(key)
    
    thesis_bib_entries.append({
        'key': key,
        'author': author_part,
        'year': year,
        'title': title,
        'journal': journal,
        'doi': doi,
        'raw': e
    })

# Format thesis entries as BibTeX
extra_bib = []
for b in thesis_bib_entries:
    doi_field = f"  doi = {{{b['doi']}}},\n" if b['doi'] else ""
    entry_tex = f"""@article{{{b['key']},
  author = {{{b['author']}}},
  title = {{{{{b['title']}}}}},
  journal = {{{b['journal']}}},
  year = {{{b['year']}}},
{doi_field}  note = {{Referenced from undergraduate thesis literature}}
}}
"""
    extra_bib.append(entry_tex)

# Foundational references
foundational = """
@article{neese_2012_orca,
  author = {Neese, Frank},
  title = {{The ORCA Program System}},
  journal = {WIREs Comput. Mol. Sci.},
  volume = {2},
  number = {1},
  pages = {73--78},
  year = {2012},
  doi = {10.1002/wcms.81}
}

@article{becke_1993_density,
  author = {Becke, Axel D.},
  title = {{Density-functional thermochemistry. III. The role of exact exchange}},
  journal = {J. Chem. Phys.},
  volume = {98},
  number = {7},
  pages = {5648--5652},
  year = {1993},
  doi = {10.1063/1.464913}
}

@article{weigend_2005_balanced,
  author = {Weigend, Florian and Ahlrichs, Reinhart},
  title = {{Balanced basis sets of split valence, triple zeta valence and quadruple zeta valence quality for H to Rn: Design and assessment of accuracy}},
  journal = {Phys. Chem. Chem. Phys.},
  volume = {7},
  number = {18},
  pages = {3297--3305},
  year = {2005},
  doi = {10.1039/B508541A}
}

@article{friedman_2001_greedy,
  author = {Friedman, Jerome H.},
  title = {{Greedy function approximation: A gradient boosting machine}},
  journal = {Ann. Statist.},
  volume = {29},
  number = {5},
  pages = {1189--1232},
  year = {2001},
  doi = {10.1214/aos/1013203451}
}

@article{hoerl_1970_ridge,
  author = {Hoerl, Arthur E. and Kennard, Robert W.},
  title = {{Ridge Regression: Biased Estimation for Nonorthogonal Problems}},
  journal = {Technometrics},
  volume = {12},
  number = {1},
  pages = {55--67},
  year = {1970},
  doi = {10.1080/00401706.1970.10488634}
}

@article{kingma_2014_adam,
  author = {Kingma, Diederik P. and Ba, Jimmy},
  title = {{Adam: A Method for Stochastic Optimization}},
  journal = {arXiv preprint arXiv:1412.6980},
  year = {2014}
}

@article{landrum_2013_rdkit,
  author = {Landrum, Greg},
  title = {{RDKit: Open-Source Cheminformatics Software}},
  journal = {Online},
  year = {2013},
  url = {https://www.rdkit.org}
}

@article{frenkel_2005_ilthermo,
  author = {Frenkel, Michael and Chirico, Robert D. and Diky, Vladimir and Dong, Qian and Marsh, Kenneth N. and Dymond, John H. and Wakeham, William A. and Stein, Stephen E. and Sievers, Edgar},
  title = {{ILThermo: A Web-Based Database for Thermodynamics Properties of Ionic Liquids}},
  journal = {J. Chem. Eng. Data},
  volume = {50},
  number = {5},
  pages = {1497--1500},
  year = {2005},
  doi = {10.1021/je0501844}
}
"""

all_bib = dataset_bib.strip() + "\n\n" + "".join(extra_bib).strip() + "\n\n" + foundational.strip() + "\n"

with open("workspace/refs.bib", "w") as f:
    f.write(all_bib)

total_keys = len(re.findall(r"^@\w+\{\s*([^,\s]+)", all_bib, re.M))
print(f"Written workspace/refs.bib with {total_keys} total entries.")

# Build citation_pool.json
pool = {}
# Parse entries from all_bib
pattern = re.compile(r"@(\w+)\{\s*([^,\s]+),([^@]*?)\n\}", re.S)
for match in pattern.finditer(all_bib):
    entry_type, key, body = match.groups()
    fields = {}
    for line in body.split("\n"):
        fm = re.match(r"^\s*(\w+)\s*=\s*[\{|\"](.*)[\}|\"],?\s*$", line)
        if fm:
            fields[fm.group(1).lower()] = fm.group(2).strip().strip("{}")
    
    title = fields.get("title", key)
    authors = fields.get("author", "Unknown")
    year = fields.get("year", "2024")
    doi = fields.get("doi", "")
    note = fields.get("note", "")
    abstract = f"{title}. By {authors} ({year}). {note}".strip()
    
    pool[key] = {
        "key": key,
        "title": title,
        "authors": authors,
        "year": year,
        "doi": doi,
        "abstract": abstract
    }

with open("workspace/citation_pool.json", "w") as f:
    json.dump(pool, f, indent=2)

print(f"Written workspace/citation_pool.json with {len(pool)} entries.")
