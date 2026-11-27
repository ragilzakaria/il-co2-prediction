# Research roadmap: CO₂ solubility in ionic liquids

Last updated: 26 November 2026. Time estimates are rough and assume HPC queue times like today's.

## The goal in one sentence

Predict CO₂ solubility in ionic liquids, **including liquids the model has never seen**, from computed properties of the ions and the liquid. Validate every computed number against the literature, and turn the workflow into a reusable QM → MD → ML pipeline for BRIN's computational chemistry team.

## The big picture

```
Phase 0  Thesis + clean-up             DONE
Phase 1  Fix the thesis-level data     DONE       → honest baseline (paper v1 numbers)
Phase 2  Choose and validate methods   NOW        → documented QM protocol
Phase 3  New QM descriptors            next       → physical descriptors for all ions and pairs
Phase 4  Molecular dynamics            parallel   → bulk-liquid properties and features
Phase 5  ML v2                         after 3–4  → feature ladder, modern baselines, uncertainty, benchmark
Phase 6  Paper + pipeline              end        → publication and pipeline handed to BRIN
```

Phase 1 is done; Phase 2 is under way. Phase 4 can start during Phase 3.

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
| 1.5 | Update the manuscript with the honest numbers | Claude + you | **done** 25 Nov (d5a7281); rawdata.xlsx source identified and all dataset DOIs verified 26 Nov (`data/README.md`) |

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

**Update (26 Nov 2026), on the curated dataset** (12,026 points, 110 ILs; rebuilt from the raw ILThermo records, inconsistent publications removed, atom counts from the corrected structures):

| Features | R² |
|---|---|
| T + P only | 0.77 / 0.80 |
| + thesis DFT descriptors | 0.84 / 0.89 (SD 0.07–0.11) |
| + corrected DFT descriptors | **0.90 / 0.93** (SD ≈ 0.02) |
| + atom counts only (no QM) | 0.87 / 0.91 |
| + corrected DFT + atom counts | 0.89 / 0.92 |

- On the clean data the corrected descriptors are **slightly better than atom counts** (better in 4 of 5 folds for both models, by 0.02 on average). The conclusion stands: the gain from isolated-ion DFT is small, and Phases 3–4 are still the way forward.

## Phase 2: choose and validate the QM method (now, about 1–2 weeks)

Goal: a QM protocol that is **proven** against high-level references before it is used on everything.

| Step | What | Check against | Status |
|---|---|---|---|
| 2.1 | Ion-pair benchmark: 66 published structures, r²SCAN-3c and ωB97X-D3/def2-TZVPD (± counterpoise) | CCSD(T) values from Karu 2018 and IL-2013. Target: mean absolute deviation < 5 kJ/mol | **done** 26 Nov: `calculations/benchmark_ion_pairs/RESULTS.md` |
| 2.2 | Test parallel ORCA on Mahameru (`openmpi4/4.1.4`, `%pal`) | speed-up versus one core | **done** 26 Nov: 3.6–3.8× on 4 cores, 11–13× on 16 |
| 2.3 | Anion orbitals with diffuse functions (ωB97X-D3/def2-TZVPD single points) | published ionisation energies (rank order) | **done** 26 Nov: anion rank correlation 0.96 (r²SCAN-3c 0.82) |
| 2.4 | Write the protocol down (`docs/QM_PROTOCOL.md`) | — | **done** 26 Nov |

**Decision point B:** which method to use for energies and for orbitals.

**Answer (26 Nov 2026): r²SCAN-3c geometries and frequencies, ωB97X-D3/def2-TZVPD single points for energies and orbitals, no counterpoise.** Mean absolute deviation from CCSD(T) over the 66 ion pairs, kJ/mol:

| Method | All 66 | Without [Cl]⁻ | [Cl]⁻ only | Time per pair |
|---|---|---|---|---|
| r²SCAN-3c | 6.0 | 4.1 | 17.9 | ~2 min |
| ωB97X-D3/def2-TZVPD | **3.7** | **3.0** | 8.2 | ~30 min |
| ωB97X-D3/def2-TZVPD + counterpoise | 4.0 | 3.2 | 9.0 | more runs |

- ωB97X-D3/def2-TZVPD meets the target for every anion; counterpoise does not help, so it is dropped.
- r²SCAN-3c is about 35 times cheaper and good enough for geometries, but overbinds chloride by 12–33 kJ/mol.
- All methods rank the pairs almost like the reference (rank correlation 0.98–0.99).
- Not yet tested: r²SCAN-3c geometries instead of the published ones, and the sulfate, phosphate, carboxylate and nitrate anions of the training set.

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

## Phase 5: ML v2 (after Phases 3–4, about 2 months)

Goal: a result that holds beyond this one dataset: which kinds of physical information let a model generalise to ions it has never seen, tested against strong modern baselines.

**Core model**
- Features: physical QM, COSMO-RS and MD descriptors, plus T and P.
- Validation: leave-IL-out and leave-anion-out only, with errors reported per family.
- Models: Gaussian-process and boosted-tree ensembles that give **error bars**.
- Physics: solubility rises with pressure and follows Henry's law at low pressure.
- **Training-range check:** flag any new ion whose descriptors lie outside the range of the training ions (tree models cannot extrapolate; [MDEGSO4] with the ωB97X-D3 LUMO is the worked example).
- **Active learning:** the model's uncertainty picks the next ions to compute (Phases 3–4), and the results go back into the model.

**5.1 Feature ladder with null tests.** Every feature level, from cheap to expensive (atom counts → r²SCAN-3c ions → ωB97X-D3 ions → COSMO-RS σ-moments and K_H → ion pairs and anion–CO₂ complexes → MD), is compared on the same leave-anion-out folds, each against its own shuffled-feature null band and over several seeds (`check_descriptor_noise.py`, `check_cosmors_feature.py`, `check_sigma_moments_feature.py`, `check_seed_spread.py` are the first instances). Written up as a general protocol: does a feature carry physical information, or does it only identify the ion?

**5.2 Modern baselines** on the same folds:
- graph neural networks on the cation and anion graphs (Chemprop or a similar message-passing model);
- pretrained molecular representations (embeddings from a pretrained molecular model, and from a machine-learned interatomic potential such as MACE-OFF);
- a tabular foundation model (TabPFN) on the same feature sets;
- the question: do large pretrained representations beat a handful of physics numbers for unseen anions?

**5.3 Uncertainty that holds for new anions.** Calibrated intervals under group shift (conformal prediction with anion-grouped calibration, compared with the GP and ensemble error bars); report coverage per held-out anion, and whether high uncertainty flags the anions that are actually predicted badly.

**5.4 Physics in the model, not only in the features.** Constraints that hold by construction (solubility monotone in pressure, Henry's-law limit at low pressure, sensible temperature dependence); compare constrained and unconstrained models on held-out anions, including whether the constraints help most where data are scarce.

**5.5 Benchmark release.** The curated dataset (12,026 points, 110 ILs, verified source per point) with fixed split files (point-wise, leave-IL-, leave-cation-, leave-anion-out), every feature level as a ready table, the null-test code, and baseline results. If time allows, widen it with other ILThermo properties or gases measured on the same ILs (for example N₂, CH₄ or H₂ solubility, density, viscosity), so the conclusions are tested on more than one target.

## Phase 6: paper and pipeline

- **Paper:** the honest-validation story (Phase 1), the protocol validated against CCSD(T) (Phase 2), physical descriptors and MD (Phases 3–4), and the feature ladder, baselines, uncertainty and benchmark (Phase 5). Framed around the general question (what information lets a model generalise to unseen components), with CO₂ in ionic liquids as the worked case.
- **Pipeline** (internship track 2): the `scripts/` and job templates become a documented QM → MD → ML pipeline for BRIN's computational chemistry team.

---

## Where things stand (27 Nov 2026)

- **Phase 1 done** (25 Nov): all 9 ions fixed, descriptors re-extracted, corrected training table built, grouped comparison run (decision A answered), manuscript updated.
- **Phase 2.1 done** (26 Nov): ion-pair benchmark finished (HPC job 575850); decision B answered. Protocol drafted in `docs/QM_PROTOCOL.md`.
- Frequency checks for [N1,8,8,8], [P6,6,6,14] and [oleyl-MIM] done on Mahameru (batch02): all 83 computed ions are verified minima.
- **Phase 2 done** (26 Nov, job 576022): parallel ORCA works with `openmpi4/4.1.4` (11–13× on 16 cores); ωB97X-D3/def2-TZVPD single points of all 83 ions (`data/processed/dft_descriptors_wb97xd3.csv`). Anion −HOMO against published ionisation energies: rank correlation 0.96 (r²SCAN-3c 0.82), every anion HOMO bound; cations 0.77 at both levels (weak test: four imidazolium IEs within 0.03 eV). Protocol finished in `docs/QM_PROTOCOL.md`.
- ωB97X-D3 descriptors in the model (`scripts/compare_orbital_levels.py`, `scripts/check_descriptor_noise.py --add-wb97xd3`; leave-one-anion-out MAE in x(CO₂)): gradient boosting 0.0509 (counts), 0.0455 (r²SCAN-3c), **0.0444** (ωB97X-D3); random forest 0.0603, 0.0561, 0.0520. The oxyanion penalty is gone ([MeSO4] +0.035 → +0.009, [TFA] +0.027 → +0.009, now inside the noise band), but they do not become better than noise. Better than noise: [Tf2N], [BETA], [NfO], [PF6], [FSA]. One anion much worse: [MDEGSO4] (+0.008 → +0.107). Cause (`scripts/check_anion_lumo.py`): its ωB97X-D3 LUMO is the lowest of all anions, outside the training range, so the trees treat it like the high-solubility fluorinated anions; the anion LUMO with a diffuse basis mostly follows anion size. Dropping the anion LUMO and gap makes the pooled error worse (0.0468; [PF6], [Tf2N] and [FSA] lose most of their gain, [TFA] is worse than noise again), so all 10 ωB97X-D3 descriptors are kept and Phase 5 needs a training-range check for new ions (`docs/QM_PROTOCOL.md` section 6).
- Step 3.2 job ready: `scripts/make_conformer_reopt.py` + `calculations/conformers_r2scan3c/submit.sh` (r²SCAN-3c optimisation of the current structure and up to 5 CREST conformers per ion, keep the lowest, frequencies). Tested locally on [BMIM] with a short CREST run: a CREST conformer came out 1.0 kcal/mol below the current structure. Run once CREST job 576149 is back. Conformer procedure in `docs/QM_PROTOCOL.md` section 2.
- Per-IL errors on the curated set (`scripts/per_il_errors.py`, 26 Nov; tables `paper/data/per_*_errors_curated.csv`, figure `images/checks/per_il_errors.png`). Gradient boosting, MAE in x(CO₂): leave-one-IL-out 0.024 (counts) and 0.025 (DFT + counts); leave-one-anion-out 0.051 and 0.046. A new anion doubles the error. Null test (`scripts/check_descriptor_noise.py`: the same descriptors shuffled among the ions, 20 times) shows that adding *any* 10 ion-level numbers makes new-anion predictions worse on average (MAE 0.059 against 0.051 for counts), while the real r²SCAN-3c descriptors make them better (0.046, better than all 20 shuffles). Per anion, the real gain lies outside the shuffled range for [Tf2N], [BETA], [TCB], [DCA], [NfO] and [methide]; for the small oxyanions ([MeSO4], [NO3], [MeSO3], [DEPO4]) the descriptors behave like random numbers (no information, not real harm), except [TFA], which is really worse. With ωB97X-D3 orbitals (bullet above) the [TFA] and oxyanion penalty disappears. Worst single ILs ([BBIM][MeSO3], [BMIM][NfO], [HOC2MIM][TFA], [HOC3MIM][NO3]) each come from one publication; error does not follow ion novelty (Spearman |ρ| ≤ 0.15).
- Step 3.4 started (27 Nov): openCOSMO-RS from ORCA 6.1.1 works (`scripts/cosmors_henry.py`; the COSMO step on the laptop, the short openCOSMORS step on the Mahameru login node). CO₂ Henry constants at 25 °C for all 110 training pairs (`calculations/cosmors/henry.csv`). Against Palomar et al. 2011 (18 ILs): rank correlation with experiment 0.89 (COSMOtherm 0.87), same scatter, but about twice too high (COSMOtherm 1.3×). Against P/x of the dataset at low pressure and 293–303 K: rank correlation 0.86 for the 54 ILs with at least 5 points, again about 2× too high. Its temperature dependence has the wrong sign, so only 25 °C values are used. **K_H as a model feature** (`scripts/check_cosmors_feature.py`, gradient boosting, leave-one-anion-out): counts 0.0509 → **0.0428** with ln K_H, ωB97X-D3 0.0444 → **0.0422**; both better than all shuffled-K_H draws (means 0.0521 and 0.0463). One COSMO-RS number per IL helps new anions as much as the ten ωB97X-D3 descriptors, and still adds on top of them; it also fixes most of the [MDEGSO4] problem (0.140 → 0.119). Worse than noise: the small oxyanions [MeSO4], [MeSO3], [NO3] (details in `docs/QM_PROTOCOL.md` section 4). **σ-moments** (7 per ion, `scripts/cosmors_sigma_moments.py`, `scripts/check_sigma_moments_feature.py`): counts 0.0509 → **0.0380**, ωB97X-D3 0.0444 → **0.0371**, ωB97X-D3 + K_H 0.0422 → 0.0377, all better than every shuffled draw; the best QM features so far. Over 5 seeds (`scripts/check_seed_spread.py`, spread 0.0003–0.0007): ωB97X-D3 + σ 0.0376, counts + σ 0.0384 (ωB97X-D3 better with every seed, by about 2 %), and K_H adds nothing on top of σ. Feature choice for Phase 5: σ-moments + ωB97X-D3 + counts; K_H stays a physics baseline. [MDEGSO4] 0.140 → 0.089; [TFA] worse on every base. Next: rerun on the step-3.2 conformers.
- Step 3.5 job ready and tested (27 Nov): `scripts/make_ion_pairs.py` + `calculations/ion_pairs/` (110 ion pairs, 23 anion–CO₂ complexes; GFN2-xTB from about 22 starts, r²SCAN-3c on the 3 lowest distinct, ωB97X-D3 single points). Laptop test on [EMIM][BF4], [EMIM][DCA] and [BF4]⁻·CO₂ (`calculations/ion_pairs/test_laptop/README.md`): −378.2 and −363.0 kJ/mol against CCSD(T) −378.7 and −368.2. The test found a bug in the duplicate check (three copies of one structure went to r²SCAN-3c, giving −371.7 for [EMIM][BF4]); fixed and confirmed. Submitted 27 Nov as job 576441 (short partition, 23 h; resubmit until done, finished runs are skipped). Started 27 Nov 16:44; stage A (2,940 GFN2-xTB runs) took 4.5 h; stage B has 359 r²SCAN-3c optimisations (1.6–5 h each for 53–66 atoms; about 23–32 h on 16 cores in all), then 359 single points. Continuations 576612 (`--dependency=afterany:576441`) and 576613 (`afterany:576612`) queued, for the rest of stage B and stage C (359 ωB97X-D3 single points, roughly 30 h). No stage B for [DCA]⁻·CO₂ and [DEPO4]⁻·CO₂: every GFN2-xTB start bonded the CO₂ to the anion (N–C 1.64 Å, O–C 1.78 Å, O–C–O 141–146°), probably GFN2-xTB overbinding; rerun both at r²SCAN-3c from the starts in the selective rerun. It starts from the current ion structures, not the CREST conformers. **After step 3.2:** list the ions whose re-optimised conformer is noticeably different (more than about 1 kcal/mol lower, or another chain shape), rerun only the pairs and CO₂ complexes containing them (`scripts/make_ion_pairs.py --geometries reopt --only ...`), and report how much the interaction energies change.
- Step 3.1 running: CREST 576149 (26 Nov 17:44, ends at the 23 h limit on 27 Nov 16:44): 4 of 83 ions done after 19 h; the rest needs about 45 h. The follow-up job keeps the same CREST settings (one search depth for all ions) and moves to the 3-day partition medium-small (finished ions are skipped). Submitted as 576440 with `--dependency=afterany:576149`; 576149 timed out on 27 Nov 16:44 with 4 of 83 ions done, and 576440 waits for a free job slot (the ion-pair job took the other one), likely 28 Nov. The ion-pair job (3.5) can use the other slot once 576212 ends. CREST 3.0.2 conformer search of all 83 ions, `calculations/conformers_crest/` (`scripts/make_conformer_search.py`). ORCA's GOAT was tried first and dropped: about 5 h on one core for [BMIM].

## Working rules

1. Validate every computed output against `data/literature_reference/` before using it.
2. New QM protocols go in their own folders under `calculations/`. Never overwrite `data/raw/` or `data/thesis/`.
3. HPC account limits: at most 2 running and 4 submitted jobs. Pack many small runs into one job.
4. Run frequencies of ions above about 30 atoms on Mahameru; the macOS ORCA crashes on them.
5. Commit after each finished step, so the repo always shows where things stand.

## Who does what

- **You:** submit and monitor HPC jobs, copy results between the laptop and the HPC, make the decisions (A, B, and what goes in the paper), and agree priorities with your mentor and supervisor.
- **Claude:** writes the inputs, scripts and job files; checks and analyses the results; updates the data, the manuscript and the repo.
