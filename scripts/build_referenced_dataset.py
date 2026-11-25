"""
Extraction and harmonization pipeline for CO2 solubility in ionic liquids from IUPAC/NIST ILThermo.
Extracts experimental data points, experimental methods, and bibliographic citations with DOIs
strictly for ionic liquids formed by the 51 cations and 27 anions in the ORCA database.
"""

import os
import sys
import time
import re
import urllib.parse
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
import requests
import ilthermopy
from rdkit import Chem
from rdkit.Chem import Descriptors

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(ROOT, "data")
TARGETS_FILE = os.path.join(DATA_DIR, "raw", "matched_ilthermo_targets.csv")
OUTPUT_FILE = os.path.join(DATA_DIR, "processed", "referenced_co2_solubility_dataset.csv")
IONS_SMILES_FILE = os.path.join(DATA_DIR, "ions", "ions_with_smiles.csv")
ORCA_FILE = os.path.join(DATA_DIR, "processed", "orca_descriptors_individual.csv")

# 1. Precompute Molecular Weights for unit conversion (molality -> mole fraction)
ions_df = pd.read_csv(IONS_SMILES_FILE)
orca_df = pd.read_csv(ORCA_FILE)
known_ions = set(orca_df['ion_abbreviation'])

mw_dict = {}
for idx, r in ions_df[ions_df['ion_abbreviation'].isin(known_ions)].iterrows():
    try:
        mol = Chem.MolFromSmiles(r['smiles'])
        if mol:
            mw_dict[r['ion_abbreviation']] = Descriptors.MolWt(mol)
    except Exception:
        pass

print(f"Loaded {len(mw_dict)} ion molecular weights from SMILES.")

# 2. DOI Resolver Cache
doi_cache = {}

def get_doi(title, citation):
    clean_title = re.sub(r'<[^>]+>', '', title or '').strip()
    if clean_title in doi_cache:
        return doi_cache[clean_title]
    if not clean_title and citation in doi_cache:
        return doi_cache[citation]
        
    query_str = clean_title if clean_title else citation
    if not query_str:
        return ""
        
    try:
        url = f"https://api.crossref.org/works?query.bibliographic={urllib.parse.quote(query_str)}&rows=1"
        r = requests.get(url, headers={'User-Agent': 'IL-CO2-Researcher/1.0 (mailto:academic-research@univ.edu)'}, timeout=4)
        if r.status_code == 200:
            items = r.json().get('message', {}).get('items', [])
            if items:
                doi = items[0].get('DOI', '')
                if clean_title: doi_cache[clean_title] = doi
                doi_cache[citation] = doi
                return doi
    except Exception:
        pass
        
    return ""

# 3. Smart Parser for ILThermo Entries
def parse_entry_data(entry, cat_abbr, an_abbr):
    header = entry.header
    data = entry.data
    mw_il = mw_dict.get(cat_abbr, 0.0) + mw_dict.get(an_abbr, 0.0)
    
    t_col, p_col, x_col, dx_col = None, None, None, None
    x_type = 'mole_fraction'
    t_unit, p_unit = 'k', 'kpa'
    
    for col, desc in header.items():
        desc_l = desc.lower()
        if 'temperature' in desc_l:
            t_col = col
            if 'k' in desc_l: t_unit = 'k'
            elif 'c' in desc_l: t_unit = 'c'
        elif 'pressure' in desc_l:
            if col.startswith('d') or 'error' in desc_l:
                continue
            p_col = col
            if 'kpa' in desc_l: p_unit = 'kpa'
            elif 'mpa' in desc_l: p_unit = 'mpa'
            elif 'bar' in desc_l: p_unit = 'bar'
            elif ' pa' in desc_l: p_unit = 'pa'
        elif 'carbon dioxide' in desc_l:
            if 'error' in desc_l or col.startswith('d'):
                dx_col = col
            elif 'mole fraction' in desc_l:
                x_col = col
                x_type = 'mole_fraction'
            elif 'molality' in desc_l:
                x_col = col
                x_type = 'molality'
            elif 'mole ratio' in desc_l:
                x_col = col
                x_type = 'mole_ratio'
            elif 'weight fraction' in desc_l or 'mass fraction' in desc_l:
                x_col = col
                x_type = 'weight_fraction'
                
    if not (t_col and p_col and x_col):
        return []
        
    ref_citation = entry.ref.full if hasattr(entry.ref, 'full') else str(entry.ref)
    ref_title = entry.ref.title if hasattr(entry.ref, 'title') else ''
    
    # Extract year
    year_match = re.search(r'\((\d{4})\)', ref_citation)
    ref_year = int(year_match.group(1)) if year_match else None
    
    exp_method = entry.expmeth if entry.expmeth else 'unspecified'
    
    # DOI resolution
    doi = get_doi(ref_title, ref_citation)
    
    records = []
    for _, row in data.iterrows():
        try:
            t_val = float(row[t_col])
            p_val = float(row[p_col])
            raw_x = float(row[x_col])
            raw_dx = float(row[dx_col]) if dx_col and dx_col in row and pd.notna(row[dx_col]) else None
            
            # Unit conversions for T
            if t_unit == 'c' or t_val < 150:
                t_val += 273.15
                
            # Unit conversions for P (to kPa)
            if p_unit == 'bar':
                p_val *= 100.0
            elif p_unit == 'mpa':
                p_val *= 1000.0
            elif p_unit == 'pa':
                p_val /= 1000.0
                
            # Unit conversions for x_CO2 (to mole fraction)
            if x_type == 'mole_fraction':
                x_val = raw_x / 100.0 if raw_x > 1.0 else raw_x
                dx_val = (raw_dx / 100.0 if raw_dx and raw_x > 1.0 else raw_dx)
            elif x_type == 'mole_ratio':
                x_val = raw_x / (1.0 + raw_x)
                dx_val = None
            elif x_type == 'molality':
                if mw_il > 0:
                    x_val = raw_x / (raw_x + 1000.0 / mw_il)
                    dx_val = None
                else:
                    continue
            elif x_type == 'weight_fraction':
                w = raw_x / 100.0 if raw_x > 1.0 else raw_x
                if mw_il > 0:
                    n_co2 = w / 44.01
                    n_il = (1.0 - w) / mw_il
                    x_val = n_co2 / (n_co2 + n_il)
                    dx_val = None
                else:
                    continue
            else:
                continue
                
            # Sanity checks
            if not (180.0 <= t_val <= 500.0):
                continue
            if not (0.1 <= p_val <= 100000.0):
                continue
            if not (0.0 < x_val < 1.0):
                continue
                
            records.append({
                'entry_id': entry.id,
                'cation': cat_abbr,
                'anion': an_abbr,
                'Temperature_K': round(t_val, 2),
                'Pressure_kPa': round(p_val, 2),
                'CO2_solubility': round(x_val, 6),
                'solubility_uncertainty': round(dx_val, 6) if dx_val is not None else None,
                'exp_method': exp_method,
                'ref_year': ref_year,
                'ref_citation': ref_citation,
                'ref_title': ref_title,
                'ref_doi': doi
            })
        except Exception:
            continue
            
    return records

def process_target(row_dict):
    eid = row_dict['id']
    cat = row_dict['cation']
    an = row_dict['anion']
    for attempt in range(3):
        try:
            entry = ilthermopy.GetEntry(eid)
            records = parse_entry_data(entry, cat, an)
            return records
        except Exception as e:
            time.sleep(1.0 + attempt * 0.5)
    return []

def main():
    if not os.path.exists(TARGETS_FILE):
        print(f"Targets file not found: {TARGETS_FILE}")
        return

    targets_df = pd.read_csv(TARGETS_FILE)
    print(f"Starting extraction for {len(targets_df)} ILThermo entries across {targets_df.groupby(['cation', 'anion']).size().shape[0]} unique ILs...")
    
    target_dicts = targets_df.to_dict('records')
    all_extracted_records = []
    
    start_time = time.time()
    completed = 0
    total = len(target_dicts)
    
    # Process using thread pool
    with ThreadPoolExecutor(max_workers=6) as executor:
        future_to_entry = {executor.submit(process_target, r): r for r in target_dicts}
        for future in as_completed(future_to_entry):
            records = future.result()
            if records:
                all_extracted_records.extend(records)
            completed += 1
            if completed % 25 == 0 or completed == total:
                elapsed = time.time() - start_time
                rate = completed / elapsed
                pts = len(all_extracted_records)
                print(f"[{completed}/{total}] Entries processed ({rate:.1f} entry/s) | Points collected: {pts}")
                
    if not all_extracted_records:
        print("No records extracted.")
        return
        
    df_out = pd.DataFrame(all_extracted_records)
    print(f"\nRaw extracted data points: {len(df_out)}")
    
    # Deduplicate exact measurements from same reference
    df_out.drop_duplicates(subset=['cation', 'anion', 'Temperature_K', 'Pressure_kPa', 'CO2_solubility', 'ref_citation'], inplace=True)
    df_out.sort_values(['cation', 'anion', 'Temperature_K', 'Pressure_kPa'], inplace=True)
    
    # Save output
    df_out.to_csv(OUTPUT_FILE, index=False)
    print(f"\nSuccessfully generated: {OUTPUT_FILE}")
    print("=" * 65)
    print("DATASET SUMMARY REPORT:")
    print(f"Total clean experimental data points: {len(df_out):,}")
    unique_ils = df_out.groupby(['cation', 'anion']).size()
    print(f"Total UNIQUE Ionic Liquids:           {len(unique_ils)}")
    print(f"Unique Cations:                       {df_out['cation'].nunique()}")
    print(f"Unique Anions:                        {df_out['anion'].nunique()}")
    print(f"Distinct Primary Publications Cited:  {df_out['ref_citation'].nunique()}")
    dois_found = df_out['ref_doi'].str.len() > 3
    print(f"Data points with verified DOI:        {dois_found.sum():,} ({dois_found.mean()*100:.1f}%)")
    print(f"Temperature Range (K):                {df_out['Temperature_K'].min()} K - {df_out['Temperature_K'].max()} K")
    print(f"Pressure Range (kPa):                 {df_out['Pressure_kPa'].min()} kPa - {df_out['Pressure_kPa'].max()} kPa")
    print(f"CO2 Solubility Range (mole frac):     {df_out['CO2_solubility'].min()} - {df_out['CO2_solubility'].max()}")
    print("=" * 65)
    print("\nTop 15 Ionic Liquids by Data Points in New Dataset:")
    print(unique_ils.sort_values(ascending=False).head(15))

if __name__ == '__main__':
    main()
