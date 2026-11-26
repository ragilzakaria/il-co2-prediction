# QM protocol

Status: **draft, 26 November 2026** (roadmap step 2.4). Decision B is made; the parts marked *pending*
are filled in after roadmap steps 2.2 (parallel ORCA) and 2.3 (anion orbitals).

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

## 2. ORCA inputs

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
- *Pending (2.2):* `%pal nprocs` setting and memory per core for Mahameru.

## 3. Derived quantities

- **Interaction (binding) energy** of a complex A·B: E_int = E(A·B) − E(A) − E(B), all at the
  ωB97X-D3/def2-TZVPD level, monomers **at the complex geometry**, **no counterpoise**
  (it did not improve agreement with CCSD(T): 4.0 against 3.7 kJ/mol).
- **Orbital energies, dipole, charges** from the ωB97X-D3/def2-TZVPD single point.
  *Pending (2.3):* check of anion −HOMO against published ionisation energies.

## 4. Checks before a number is used

1. Optimisation converged (`THE OPTIMIZATION HAS CONVERGED`) and `ORCA TERMINATED NORMALLY`.
2. No imaginary frequencies; a mode below 50 cm⁻¹ in a small ion is a warning sign.
3. Closed shell with the expected charge (wrong SMILES once gave doublet radicals).
4. The ion did not fall apart during optimisation (compare bonds with the start structure; [PF₆]⁻
   once lost a fluoride).
5. Compare with `data/literature_reference/` where a value exists (roadmap working rule 1).

## 5. Known limits

- Benchmarked on published (CCSD(T)-quality) geometries; the effect of using r²SCAN-3c geometries
  has not been measured yet.
- The benchmark has no sulfate, phosphate, carboxylate or nitrate anions, which occur in the training
  data and are among the hardest anions for the model.
- r²SCAN-3c overbinds chloride ion pairs by 12–33 kJ/mol; use only the single-point level for
  energies involving chloride.
- Gas-phase species, one conformer each, until the conformer search of roadmap step 3.1.
