"""
Merges the referenced NIST experimental dataset with per-ion DFT descriptors and elemental
atom counts, producing the model-ready training table with full provenance.

Two descriptor sets:
  r2scan3c (default)  data/processed/dft_descriptors_r2scan3c.csv, corrected parser and the
                      reruns of 25 Nov 2026  ->  data/processed/final_referenced_model_dataset_r2scan3c.csv
  thesis              data/processed/orca_descriptors_individual.csv, the values the thesis model used
                      ->  data/processed/final_referenced_model_dataset.csv (reproduces the existing file)

Usage, from the repository root:
    python scripts/merge_referenced_model_dataset.py               # corrected descriptors
    python scripts/merge_referenced_model_dataset.py --thesis      # thesis descriptors
"""

import argparse
import os
import pandas as pd
import numpy as np

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
EXP_FILE = os.path.join(DATA_DIR, "processed", "referenced_co2_solubility_dataset.csv")
DESC_FILE_THESIS = os.path.join(DATA_DIR, "processed", "orca_descriptors_individual.csv")
DESC_FILE = os.path.join(DATA_DIR, "processed", "dft_descriptors_r2scan3c.csv")
ATOM_FILE = os.path.join(DATA_DIR, "ions", "ions_with_atomic_counts.csv")
OUT_FILE_THESIS = os.path.join(DATA_DIR, "processed", "final_referenced_model_dataset.csv")
OUT_FILE = os.path.join(DATA_DIR, "processed", "final_referenced_model_dataset_r2scan3c.csv")
DESCRIPTOR_COLUMNS = ["ion_abbreviation", "Energy_Hartree", "HOMO_Eh", "LUMO_Eh", "HOMO_LUMO_Gap_Eh",
                      "Dipole_Moment_Debye"]

def standardize_name(series):
    return series.str.strip().str.replace(r'[+\[\]()-]', '', regex=True)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--thesis", action="store_true", help="use the thesis descriptor file")
    args = ap.parse_args()
    desc_file = DESC_FILE_THESIS if args.thesis else DESC_FILE
    out_file = OUT_FILE_THESIS if args.thesis else OUT_FILE

    if not os.path.exists(EXP_FILE):
        print(f"Experimental dataset not found yet: {EXP_FILE}")
        return

    print("Loading datasets...")
    df_exp = pd.read_csv(EXP_FILE)
    df_desc = pd.read_csv(desc_file)
    if not args.thesis:
        df_desc = df_desc[DESCRIPTOR_COLUMNS]
    df_atoms = pd.read_csv(ATOM_FILE)

    # 1. Ensure HOMO-LUMO gaps are calculated for individual ions
    if 'HOMO_LUMO_Gap_Eh' not in df_desc.columns:
        df_desc['HOMO_LUMO_Gap_Eh'] = (df_desc['LUMO_Eh'] - df_desc['HOMO_Eh']).abs()

    # 2. Combine ORCA descriptors with atom counts
    atom_cols_drop = [c for c in ['ion_name', 'ion_type', 'smiles'] if c in df_atoms.columns]
    df_atoms_only = df_atoms.drop(columns=atom_cols_drop)
    
    df_all_features = pd.merge(df_desc, df_atoms_only, on='ion_abbreviation', how='left')
    df_all_features = pd.merge(df_all_features, df_atoms[['ion_abbreviation', 'ion_name', 'ion_type']], on='ion_abbreviation', how='left')

    # 3. Standardize merge keys
    df_all_features['merge_key'] = standardize_name(df_all_features['ion_abbreviation'])
    df_exp['cation_merge_key'] = standardize_name(df_exp['cation'])
    df_exp['anion_merge_key'] = standardize_name(df_exp['anion'])

    # Separate into cation and anion feature tables
    df_cat = df_all_features[df_all_features['ion_type'] == 'cation'].copy()
    df_an = df_all_features[df_all_features['ion_type'] == 'anion'].copy()

    df_cat.columns = ['cation_' + col if col not in ['ion_abbreviation', 'merge_key'] else col for col in df_cat.columns]
    df_an.columns = ['anion_' + col if col not in ['ion_abbreviation', 'merge_key'] else col for col in df_an.columns]

    # 4. Perform inner merge with experimental data
    df_merged = pd.merge(df_exp, df_cat, left_on='cation_merge_key', right_on='merge_key', how='inner')
    df_final = pd.merge(df_merged, df_an, left_on='anion_merge_key', right_on='merge_key', how='inner')

    # Drop temporary merge keys
    drop_temp = [c for c in ['cation_merge_key', 'anion_merge_key', 'merge_key_x', 'merge_key_y', 'cation_ion_abbreviation', 'anion_ion_abbreviation'] if c in df_final.columns]
    df_final.drop(columns=drop_temp, inplace=True)
    if args.thesis:
        df_final.fillna(0, inplace=True)  # thesis behaviour, kept so the old file is reproduced
    else:
        desc = [c for c in df_final.columns if c.endswith(("_Eh", "_Debye", "_Hartree"))]
        missing = df_final[desc].isna().any(axis=1)
        if missing.any():
            bad = df_final.loc[missing, ["cation", "anion"]].drop_duplicates().values.tolist()
            raise SystemExit(f"Missing descriptors for {missing.sum()} rows, pairs: {bad}")
        df_final.fillna(0, inplace=True)  # atom counts only

    # 5. Organize column order
    id_cols = ['cation', 'cation_ion_name', 'anion', 'anion_ion_name']
    ref_cols = [c for c in ['entry_id', 'ref_year', 'ref_citation', 'ref_title', 'ref_doi', 'exp_method', 'solubility_uncertainty'] if c in df_final.columns]
    atom_count_cols = sorted([col for col in df_final.columns if 'num_' in col])
    condition_cols = ['Temperature_K', 'Pressure_kPa']
    descriptor_cols = sorted([col for col in df_final.columns if col.endswith(('_Eh', '_Debye', '_au', '_A3', '_Hartree', '_Gap_Eh'))])
    target_col = ['CO2_solubility']

    final_order = id_cols + ref_cols + condition_cols + descriptor_cols + atom_count_cols + target_col
    final_order_existing = [col for col in final_order if col in df_final.columns]
    df_final = df_final[final_order_existing]

    df_final.to_csv(out_file, index=False)
    print(f"\nFinal model-ready dataset saved to: {out_file}")
    print(f"Shape: {df_final.shape[0]} rows, {df_final.shape[1]} columns.")
    print(f"Unique IL pairs in final model dataset: {df_final.groupby(['cation', 'anion']).size().shape[0]}")

if __name__ == '__main__':
    main()
