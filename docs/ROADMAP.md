# Research roadmap: CO₂ solubility in ionic liquids

Last updated: 25 November 2026. Time estimates are rough and assume HPC queue times like today's.

## The goal in one sentence

Predict CO₂ solubility in ionic liquids, **including liquids the model has never seen**, from computed properties of the ions and the liquid. Validate every computed number against the literature, and turn the workflow into a reusable QM → MD → ML pipeline for BRIN's computational chemistry team.

## The big picture

```
Phase 0  Thesis + clean-up             DONE
Phase 1  Fix the thesis-level data     NOW        → honest baseline (paper v1 numbers)
Phase 2  Choose and validate methods   NOW        → documented QM protocol
Phase 3  New QM descriptors            next       → physical descriptors for all ions and pairs
Phase 4  Molecular dynamics            parallel   → bulk-liquid properties and features
Phase 5  ML v2                         after 3–4  → final model with uncertainty
Phase 6  Paper + pipeline              end        → publication and pipeline handed to BRIN
```

Phases 1 and 2 run at the same time now. Phase 4 can start during Phase 3.

---

## Phase 0: done (thesis and November clean-up)

- Thesis Random Forest: R² = 0.986 on a random split, but about 0.94 when whole ionic liquids are held out, and about 0.83 when whole anions are held out.
- Problems found and documented:
  - the descriptor parser read the wrong lines
  - six ions were built from the wrong SMILES
  - [PF6] broke apart during its optimisation
  - the random train/test split was optimistic
- Referenced dataset: 12,297 points from NIST ILThermo, with a citation for every point.
- Literature reference set: 518 published computed values (`data/literature_reference/`).
- Repo restructured, ORCA working on the laptop and on Mahameru.

## Phase 1: fix the thesis-level data (now, about 1 week)

Goal: an **honest baseline**, using the same model type as the thesis but with correct inputs and honest validation.

| Step | What | Who | Status |
|---|---|---|---|
| 1.1 | Rerun the 9 broken ions | ORCA (laptop + HPC) | **done** 25 Nov |
| 1.2 | Re-extract descriptors; literature check | Claude | **done**: anion HOMO rank ρ = 0.82 (thesis −0.25) |
| 1.3 | Rebuild the training table with corrected descriptors | Claude | **done**: `final_referenced_model_dataset_r2scan3c.csv` |
| 1.4 | Retrain with leave-IL-out and leave-anion-out validation; compare thesis descriptors, corrected descriptors, and T + P only | Claude (laptop) | **done** 25 Nov: `scripts/compare_descriptor_sets.py`, `paper/data/descriptor_comparison.csv` |
| 1.5 | Update manuscript sections 5–7 with the honest numbers | Claude + you | after 1.4 |

**Decision point A:** do the DFT descriptors beat T + P alone for unseen anions? If they barely do, that motivates Phase 3: better, more physical descriptors.

**Answer (25 Nov 2026), leave-anion-out R², mean of 5 folds (Random Forest / gradient boosting):**

| Features | R² |
|---|---|
| T + P only | 0.78 / 0.80 |
| + thesis DFT descriptors | 0.81 / 0.86 (unstable across folds, SD ≈ 0.1) |
| + corrected DFT descriptors | **0.89 / 0.91** (stable, SD ≈ 0.02–0.03) |
| + atom counts only (no QM) | 0.88 / 0.90 |
| + corrected DFT + atom counts | 0.89 / 0.91 |

- The corrections matter: the corrected descriptors beat the thesis ones by 0.05–0.08 and are far more stable across folds.
- The isolated-ion DFT descriptors are **no better than simple atom counts.** At this level they mostly encode ion size and identity.
- **Conclusion: go ahead with Phase 3 and Phase 4.** Getting beyond atom counts needs liquid-level physics: COSMO-RS, ion-pair and anion–CO₂ binding, and MD.

## Phase 2: choose and validate the QM method (now, about 1–2 weeks)

Goal: a QM protocol that is **proven** against high-level references before it is used on everything.

| Step | What | Check against |
|---|---|---|
| 2.1 | Ion-pair benchmark: 66 published structures, r²SCAN-3c and ωB97X-D3/def2-TZVPD (± counterpoise) | CCSD(T) values from Karu 2018 and IL-2013. Target: mean absolute deviation < 5 kJ/mol |
| 2.2 | Test parallel ORCA on Mahameru (`openmpi4/4.1.4`, `%pal`) | speed-up versus one core |
| 2.3 | Anion orbitals with diffuse functions (ωB97X-D3/def2-TZVPD single points) | published ionisation energies (rank order) |
| 2.4 | Write the protocol down (`docs/QM_PROTOCOL.md`) | — |

**Decision point B:** which method to use for energies and for orbitals.

## Phase 3: new QM descriptors (about 3–5 weeks, HPC)

Goal: replace the 236 sparse atom-count features with a **small set of physical descriptors**.

| Step | What | Gives |
|---|---|---|
| 3.1 | CREST/GFN2-xTB conformer search for all 84 ions (and new candidate ions) | realistic lowest-energy structures |
| 3.2 | r²SCAN-3c optimisation and frequencies of the lowest conformers | geometries, confirmed minima |
| 3.3 | Single points with the method chosen in Phase 2 | HOMO/LUMO, ESP extremes, molecular volume, dipole, polarisability |
| 3.4 | COSMO-RS (openCOSMO-RS in ORCA) | σ-profiles → σ-moments (compact features), plus COSMO-RS CO₂ Henry constants as a physics baseline |
| 3.5 | Ion pairs and anion–CO₂ complexes for the 109 training pairs | cation–anion and anion–CO₂ binding energies |
| 3.6 | *Optional:* clusters of 2–4 ion pairs, with and without CO₂, for a few representative liquids (e.g. [BMIM][PF6], [BMIM][Tf2N], [BMPYR][Tf2N] and one poorly predicted liquid). Start from MD snapshots (Phase 4) or CREST cluster searches; r²SCAN-3c optimisation, then single points with the Phase 2 method. | how much binding energies and orbitals change from an isolated pair towards the liquid; a check on whether ion-pair descriptors are good enough |

Every step is checked against `data/literature_reference/` before its results are used.

### Levels of description, from cheap to realistic

| Level | System | Tool | Status |
|---|---|---|---|
| 1 | Isolated ion, gas phase | ORCA | done (thesis) |
| 2 | Ion pair (+ CO₂) | ORCA | Phases 2.1 and 3.5 |
| 3 | Clusters of 2–4 ion pairs (+ CO₂) | ORCA | Phase 3.6, optional |
| 4–5 | Continuum solvent; COSMO-RS (liquid-phase activity, CO₂ Henry constants) | ORCA / openCOSMO-RS | Phase 3.4 |
| 6 | Bulk liquid, periodic box, finite temperature | GROMACS (classical MD) | Phase 4 |
| 6b | Small periodic boxes with ab initio MD, to check the force field | CP2K | optional |

## Phase 4: molecular dynamics (about 1–2 months; can start during Phase 3)

Goal: properties of the **bulk liquid**, because solubility is a liquid property, not an isolated-ion property.

| Step | What | Check against |
|---|---|---|
| 4.1 | Force-field set-up (CL&P, or polarisable CL&Pol) for 3 benchmark liquids: [BMIM][PF6], [BMIM][BF4], [BMIM][Tf2N] | experimental and Doherty 2017 densities and viscosities |
| 4.2 | CO₂ Henry constants from free-energy calculations | Kerlé 2017 (28 bar for [BMIM][Tf2N] at 300 K) and experiment |
| 4.3 | Extend to the liquids the model predicts worst (from the per-IL errors in Phase 1) | experiment |
| 4.4 | MD features for ML: free volume, CO₂ coordination, radial distribution functions | — |

GROMACS 2023.3 is installed on Mahameru. Still to check: whether GPU nodes are available to your account.

## Phase 5: ML v2 (after Phases 3–4, about 3–4 weeks)

- Features: physical QM, COSMO-RS and MD descriptors, plus T and P.
- Validation: leave-IL-out and leave-anion-out only, with errors reported per family.
- Models: Gaussian-process and boosted-tree ensembles that give **error bars**. Graph neural networks if time allows.
- Physics: solubility rises with pressure and follows Henry's law at low pressure.
- **Active learning:** the model's uncertainty picks the next ions to compute (Phases 3–4), and the results go back into the model.

## Phase 6: paper and pipeline

- **Paper:** the honest-validation story (Phase 1), the protocol validated against CCSD(T) (Phase 2), and physical descriptors and MD (Phases 3–5).
- **Pipeline** (internship track 2): the `scripts/` and job templates become a documented QM → MD → ML pipeline for BRIN's computational chemistry team.

---

## Where things stand (25 Nov 2026)

- Phase 1.1–1.4 done: all 9 ions fixed, descriptors re-extracted, corrected training table built, grouped comparison run (decision A answered). Next: 1.5, update the manuscript.
- HPC job **575850** (ion-pair benchmark, Phase 2.1) is queued.

## Working rules

1. Validate every computed output against `data/literature_reference/` before using it.
2. New QM protocols go in their own folders under `calculations/`. Never overwrite `data/raw/` or `data/thesis/`.
3. HPC account limits: at most 2 running and 4 submitted jobs. Pack many small runs into one job.
4. Run frequencies of ions above about 30 atoms on Mahameru; the macOS ORCA crashes on them.
5. Commit after each finished step, so the repo always shows where things stand.

## Who does what

- **You:** submit and monitor HPC jobs, copy results between the laptop and the HPC, make the decisions (A, B, and what goes in the paper), and agree priorities with your mentor and supervisor.
- **Claude:** writes the inputs, scripts and job files; checks and analyses the results; updates the data, the manuscript and the repo.
