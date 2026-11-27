# QM protocol

Status: **26 November 2026** (roadmap step 2.4). Decision B is made; roadmap steps 2.2 (parallel
ORCA) and 2.3 (anion orbitals) are done and filled in below.

This is the recipe every new quantum-chemistry number in the project follows: ions, ion pairs,
anion–CO₂ complexes and new candidate ions. Program: ORCA 6.1.1 (Mahameru module
`nuclear/orca/6.1.1`; the macOS build gives identical energies but crashes in frequency steps above
about 30 atoms).

## 1. Two levels

| Purpose | Level | Why |
|---|---|---|
| Geometries and frequencies | **r²SCAN-3c** | Cheap (about 2 min per ion pair on one core), built-in D4 dispersion and gCP; ion-pair interaction energies within 4 kJ/mol of CCSD(T) except for chloride. |
| Energies and orbitals (single points on the r²SCAN-3c geometry) | **ωB97X-D3/def2-TZVPD**, RIJCOSX | Ion-pair interaction energies 3.7 kJ/mol from CCSD(T) (mean absolute deviation, 66 pairs), including chloride; range-separated hybrid with diffuse functions, so anion orbitals are bound. |

Evidence: `calculations/benchmark_ion_pairs/RESULTS.md`.

## 2. Conformers (roadmap steps 3.1–3.2)

Flexible ions (alkyl chains, ether and hydroxy groups) are not used in the shape their SMILES
happened to give. For every ion:

1. **Search:** CREST 3.0.2 (iMTD-GC, GFN2-xTB, gas phase), from the current r²SCAN-3c structure;
   `--quick` for ions above 50 atoms. Job: `calculations/conformers_crest/`, on the 3-day partition
   (16 h for [C12MIM], 49 atoms, on 4 threads; more than 19 h for [P6,6,6,14] and [N1,8,8,8]).
2. **Re-rank:** r²SCAN-3c optimisation of the current structure plus up to 5 CREST conformers within
   3 kcal/mol (GFN2-xTB), duplicates removed (within 0.1 kcal/mol and 0.3 Å heavy-atom RMSD, mirror
   images included). The lowest at the r²SCAN-3c level is kept, so the result can never be worse
   than the starting structure. Job: `calculations/conformers_r2scan3c/`
   (`scripts/make_conformer_reopt.py`).
3. **Confirm:** frequencies of the kept conformer (section 3), then the single points of section 3
   on it.

GOAT (ORCA's own search) was tried and dropped: about 5 h on one core for [BMIM]. Cite CREST:
Pracht, Bohle, Grimme, *Phys. Chem. Chem. Phys.* **22**, 7169 (2020), doi:10.1039/C9CP06869D;
Pracht et al., *J. Chem. Phys.* **160**, 114110 (2024), doi:10.1063/5.0197592.

## 3. ORCA inputs

Geometry optimisation and frequencies:

```
! r2SCAN-3c Opt Freq MiniPrint
%maxcore 6000
* xyzfile <charge> 1 <ion>.xyz
```

- Use `TightOpt DefGrid3` for small, stiff anions with a soft mode (as for [C₃F₇CO₂]⁻) when the
  normal settings leave a small imaginary frequency.
- Run frequencies of species above about 30 atoms on Mahameru.

Single point:

```
! wB97X-D3 def2-TZVPD def2/J RIJCOSX TightSCF
%maxcore 2500
* xyzfile <charge> 1 <ion>_opt.xyz
```

- Charges: cation +1, anion −1, ion pair 0, anion–CO₂ complex −1; all closed-shell singlets.

Running on Mahameru (roadmap step 2.2, job 576022): parallel ORCA needs `module load openmpi4/4.1.4`
next to `nuclear/orca/6.1.1` and `%pal nprocs N end` in the input. `%maxcore` is per core, so
`%maxcore 2500` with 16 cores asks for 40 GB. Single point, wall time
(`calculations/ions_wb97xd3/parallel_timing.csv`):

| Ion | Atoms | 1 core | 4 cores | 16 cores | Speed-up 4 / 16 |
|---|---|---|---|---|---|
| [Tf2N]⁻ | 15 | 798 s | 222 s | 70 s | 3.6 / 11.4 |
| [BMIM]⁺ | 25 | 749 s | 204 s | 65 s | 3.7 / 11.5 |
| [P4,4,4,4]⁺ | 53 | 3119 s | 833 s | 235 s | 3.8 / 13.3 |

- **Many small runs:** serial ORCA, 16 side by side in one 16-core job (no lost efficiency).
- **A few large runs, or the largest of a batch:** `%pal nprocs 16` (efficiency 0.7–0.8), or
  `nprocs 4` (0.9) side by side. In job 576022 the 89 single points ran serially, and the four largest
  ions ([P6,6,6,14]⁺ 3.9 h, [N1,8,8,8]⁺ 2.4 h, [oleyl-mMIM]⁺, [C8F13MIM]⁺) set the 5.6 h wall time; run
  those in parallel first.
- Parallel ORCA does not run on the macOS laptop (its Open MPI 5 does not match the ORCA build).

## 4. Derived quantities

- **Interaction (binding) energy** of a complex A·B: E_int = E(A·B) − E(A) − E(B), all at the
  ωB97X-D3/def2-TZVPD level, monomers **at the complex geometry**, **no counterpoise**
  (it did not improve agreement with CCSD(T): 4.0 against 3.7 kJ/mol).
- **Orbital energies, dipole, charges** from the ωB97X-D3/def2-TZVPD single point
  (`data/processed/dft_descriptors_wb97xd3.csv`, all 83 ions; `scripts/extract_wb97xd3_descriptors.py`).

Check against published ionisation energies (roadmap step 2.3, `scripts/check_descriptors_vs_literature.py`,
`data/processed/checks/descriptors_vs_literature.csv`), −HOMO vs IE:

| Ions | n | r²SCAN-3c rank / r | ωB97X-D3/def2-TZVPD rank / r |
|---|---|---|---|
| Anions | 14 | 0.82 / 0.88 | **0.96 / 0.96** |
| Cations | 9 | 0.77 / 0.72 | 0.77 / 0.78 |

- Every anion HOMO is bound at ωB97X-D3/def2-TZVPD (−HOMO 3.3 eV for [SCN]⁻ to 8.3 eV for [TCB]⁻);
  at r²SCAN-3c they span −0.4 to 4.9 eV, and [NO₃]⁻ is unbound.
- The cation check is weak: four imidazolium IEs lie within 0.03 eV, and [P6,6,6,14]⁺ has the lowest
  −HOMO but not the lowest IE.

**The anion LUMO is not an electron affinity.** Every anion LUMO is positive at both levels (the
extra electron would not be bound). With the diffuse def2-TZVPD basis it is a spread-out state whose
energy mostly follows the size of the anion (rank correlation with the number of atoms −0.82, against
−0.56 at r²SCAN-3c). It is still useful to the model (section 6), but read it as a size-like number,
not as chemistry.

## 5. Checks before a number is used

1. Optimisation converged (`THE OPTIMIZATION HAS CONVERGED`) and `ORCA TERMINATED NORMALLY`.
2. No imaginary frequencies; a mode below 50 cm⁻¹ in a small ion is a warning sign.
3. Closed shell with the expected charge (wrong SMILES once gave doublet radicals).
4. The ion did not fall apart during optimisation (compare bonds with the start structure; [PF₆]⁻
   once lost a fluoride).
5. Compare with `data/literature_reference/` where a value exists (roadmap working rule 1).

## 6. Known limits

- Benchmarked on published (CCSD(T)-quality) geometries. Tested with our own r²SCAN-3c geometries on
  two pairs (`calculations/ion_pairs/test_laptop/`): [EMIM][BF4] −378.2 kJ/mol against −378.7 for
  CCSD(T) (1.6 kJ/mol from the published-geometry value), [EMIM][DCA] −363.0 against −368.2 (0.6 kJ/mol
  from the published-geometry value). The geometry level is not the limiting error; finding the
  right structure is (next point).
- **Ion-pair structure search** (`calculations/ion_pairs/pipeline.py`). GFN2-xTB ranks contact pairs
  poorly: the [EMIM][BF4] structure that is lowest at r²SCAN-3c lies 0.14 kcal/mol *above* the xTB
  lowest, and the xTB lowest is 7.1 kJ/mol higher at r²SCAN-3c. Re-optimising three distinct xTB
  structures catches this, but only if the duplicate check ignores swapped equivalent atoms; the
  first version did not, sent three copies of one structure to r²SCAN-3c, and gave −371.7 kJ/mol.
  Fixed (sorted interatomic distances) before the HPC run.
- The benchmark has no sulfate, phosphate, carboxylate or nitrate anions, which occur in the training
  data and are among the hardest anions for the model.
- r²SCAN-3c overbinds chloride ion pairs by 12–33 kJ/mol; use only the single-point level for
  energies involving chloride.
- Gas-phase species, lowest conformer only (no Boltzmann averaging). Long chains may fold onto
  the charged ring in the gas phase, which they would not do in the liquid; checked after step 3.2.
- **Descriptor values outside the training range.** Tree models cannot extrapolate: an ion whose
  descriptor lies beyond every training ion is treated like the ions at that edge. Example: the
  ωB97X-D3 LUMO of [MDEGSO4]⁻ (3.74 eV, largest anion, 23 atoms) is the lowest of the 23 anions,
  so in leave-one-anion-out it is grouped with the next-lowest ([methide]⁻, [TOS]⁻, [NfO]⁻,
  [BETA]⁻; all but [TOS]⁻ with high CO₂ solubility) and predicted too high: mean x(CO₂) 0.53
  against 0.39 measured, error 0.140 against 0.041 with r²SCAN-3c. Taking only that one value from
  r²SCAN-3c brings it back to 0.038; no other descriptor matters (`scripts/check_anion_lumo.py`,
  `paper/data/anion_lumo_swap.csv`). Dropping the anion LUMO and gap is not the fix: the pooled
  error gets worse (0.0468 against 0.0444), [PF6]⁻, [Tf2N]⁻ and [FSA]⁻ lose most of their gain,
  [TFA]⁻ is worse than noise again, and [MDEGSO4]⁻ is still worse than noise (+0.029)
  (`paper/data/anion_lumo_dropped.csv`). Rule: before trusting a
  prediction for a new ion, check that each descriptor lies inside the training range.
