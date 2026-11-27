"""Why does [MDEGSO4] get much worse with the wB97X-D3 descriptors, and is the anion LUMO worth keeping?

scripts/check_descriptor_noise.py --add-wb97xd3 showed that switching the 10 ion descriptors from
r2SCAN-3c to wB97X-D3/def2-TZVPD raises the leave-one-anion-out error of [MDEGSO4] from 0.041 to
0.140 in x(CO2), far outside the shuffled-descriptor noise band. Two tests, gradient boosting with the
settings of scripts/check_descriptor_noise.py:

1. Swap test, [MDEGSO4] fold only: start from one level and take a single anion property from the
   other level. The property whose swap moves the whole error is the cause.
2. Without the anion LUMO and HOMO-LUMO gap: wB97X-D3 descriptors minus these two columns, all 23
   anion folds, placed against the noise band of paper/data/descriptor_noise_summary.csv.
Also printed: the anions nearest to [MDEGSO4] in standardised anion-descriptor space, and how the
anion LUMO follows the number of atoms at each level.

Output: paper/data/anion_lumo_swap.csv        [MDEGSO4] fold: MAE and mean prediction per variant
        paper/data/anion_lumo_dropped.csv     per anion: effect of wB97X-D3 without anion LUMO/gap

Run from the repository root (about 15 min; one thread, because the OpenMP threads of gradient
boosting stall when ORCA runs on the same machine):
    python scripts/check_anion_lumo.py
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_validation import DFT, SEED, STATE, TARGET, load  # noqa: E402
from check_descriptor_noise import per_anion_mae  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
DATA = os.path.join(ROOT, "data", "processed", "model_dataset_curated.csv")
WB97X = os.path.join(ROOT, "data", "processed", "dft_descriptors_wb97xd3.csv")
MANIFEST = os.path.join(ROOT, "calculations", "conformers_crest", "manifest.csv")
OUT = os.path.join(ROOT, "paper", "data")
PROPS = ["Energy_Hartree", "HOMO_Eh", "LUMO_Eh", "HOMO_LUMO_Gap_Eh", "Dipole_Moment_Debye"]
HARTREE_EV = 27.211386
ANION = "[MDEGSO4]"


def main():
    df, _, counts = load(DATA)
    wb = pd.read_csv(WB97X).set_index("ion_abbreviation")
    for side in ("cation", "anion"):
        for p in PROPS:
            df[f"{side}_{p}__wb"] = df[side].map(wb[p])
    y = df[TARGET].to_numpy(float)
    wbcols = [f"{c}__wb" for c in DFT]

    # 1. swap test on the [MDEGSO4] fold
    te = (df.anion == ANION).to_numpy()
    rows = []

    def fold(cols, label):
        X = df[cols].to_numpy(float)
        m = HistGradientBoostingRegressor(max_iter=500, learning_rate=0.05, random_state=SEED)
        p = m.fit(X[~te], y[~te]).predict(X[te])
        rows.append(dict(variant=label, MAE=np.abs(p - y[te]).mean(), mean_pred=p.mean(), mean_true=y[te].mean()))
        print(f"{label:48s} MAE {rows[-1]['MAE']:.3f}  mean prediction {p.mean():.3f} (true {y[te].mean():.3f})",
              flush=True)

    print(f"{ANION} fold, gradient boosting:")
    fold(STATE + counts, "counts")
    fold(STATE + DFT + counts, "r2SCAN-3c")
    fold(STATE + wbcols + counts, "wB97X-D3")
    for p in PROPS:
        fold(STATE + [f"anion_{p}" if c == f"anion_{p}__wb" else c for c in wbcols] + counts,
             f"wB97X-D3, anion {p} from r2SCAN-3c")
    for p in PROPS:
        fold(STATE + [f"anion_{p}__wb" if c == f"anion_{p}" else c for c in DFT] + counts,
             f"r2SCAN-3c, anion {p} from wB97X-D3")
    os.makedirs(OUT, exist_ok=True)
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "anion_lumo_swap.csv"), index=False, float_format="%.4f")

    # where [MDEGSO4] sits: nearest anions, and the LUMO against size
    mean_x = df.groupby("anion")[TARGET].mean()
    n_atoms = pd.read_csv(MANIFEST).set_index("ion_abbreviation").n_atoms
    print(f"\nNearest anions (standardised HOMO, LUMO, gap, dipole); mean x(CO2) of {ANION}: {mean_x[ANION]:.3f}")
    for lab, cols in (("r2SCAN-3c", [f"anion_{p}" for p in PROPS[1:]]),
                      ("wB97X-D3", [f"anion_{p}__wb" for p in PROPS[1:]])):
        t = df.groupby("anion")[cols].first()
        z = (t - t.mean()) / t.std()
        d = np.sqrt(((z - z.loc[ANION]) ** 2).sum(1)).drop(ANION).sort_values()
        print(f"  {lab:10s} " + ", ".join(f"{a} (x {mean_x[a]:.3f})" for a in d.index[:4]))
        lumo = t[cols[1]] * HARTREE_EV
        print(f"  {'':10s} anion LUMO {lumo.min():.2f} to {lumo.max():.2f} eV, lowest {lumo.idxmin()}; "
              f"rank correlation with number of atoms {lumo.rank().corr(n_atoms[lumo.index].rank()):+.2f}")

    # 2. wB97X-D3 without the anion LUMO and gap, all folds
    cols = [c for c in wbcols if c not in ("anion_LUMO_Eh__wb", "anion_HOMO_LUMO_Gap_Eh__wb")]
    mae, pooled = per_anion_mae(df, STATE + cols + counts, y)
    s = pd.read_csv(os.path.join(OUT, "descriptor_noise_summary.csv")).set_index("anion")
    s["effect_wb97xd3_no_anion_lumo"] = mae - s.MAE_counts
    s["outside_null_no_anion_lumo"] = ~s.effect_wb97xd3_no_anion_lumo.between(s.null_p05, s.null_p95)
    keep = ["n_points", "MAE_counts", "null_p05", "null_p95", "effect_r2scan3c", "effect_wb97xd3",
            "effect_wb97xd3_no_anion_lumo", "outside_null_no_anion_lumo"]
    s = s[keep].sort_values("effect_wb97xd3_no_anion_lumo")
    s.reset_index().to_csv(os.path.join(OUT, "anion_lumo_dropped.csv"), index=False, float_format="%.4f")
    print(f"\nAll anion folds, pooled MAE: wB97X-D3 without anion LUMO/gap {pooled:.4f}")
    print(s.round(3).to_string())
    for lev in ("r2scan3c", "wb97xd3", "wb97xd3_no_anion_lumo"):
        e = s[f"effect_{lev}"]
        print(f"{lev:22s} better than the noise band: {list(s.index[e < s.null_p05])}; "
              f"worse: {list(s.index[e > s.null_p95])}")


if __name__ == "__main__":
    main()
