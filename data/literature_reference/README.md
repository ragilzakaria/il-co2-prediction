# Literature reference values for validating our calculations

Published **computed** values for the ions and ion pairs in our dataset. Before trusting a new
ORCA, COSMO-RS or MD result (on the HPC or locally), compare it with the matching row here.

Collected 25 November 2026. Each value was read from the source file or table named in the `location`
column. Nothing was typed from memory.

## Files

| File | What it is |
|---|---|
| `reference_values.csv` | All 518 reference values, one per row. Gives system, quantity, value, units, method, conditions, the experimental value when the source reports one, the source, its DOI and the table it came from. |
| `coverage_by_pair.csv` | Lists each of the 109 training pairs with its number of data points and which reference quantities exist for it. |
| `coverage_by_ion.csv` | The same for each of the 84 ions. |
| `raw/` | The values as extracted, one file per source. |
| `build_reference_table.py` | Rebuilds the three CSVs from `raw/` (`python data/literature_reference/build_reference_table.py`). |

## Coverage

- **54 of 109 training pairs** have at least one reference. Those pairs hold **9,564 of 12,297 data points (78%)**, including all of the most-measured pairs ([BMIM][PF6], [BMIM][BF4], [BMIM][Tf2N], [HMIM][Tf2N] and others).
- **Ion-pair quantum chemistry:** 37 pairs.
- **COSMO-RS CO₂ Henry constants:** 28 pairs.
- **MD of the bulk liquid:** 22 pairs.
- **Single-ion orbital or ionisation-energy references:** 28 of the 61 training ions.

## Sources

| Source | What it gives | Level of theory | Our systems |
|---|---|---|---|
| Karu et al. 2018, *Int. J. Quantum Chem.* 118, e25582. [doi:10.1002/qua.25582](https://doi.org/10.1002/qua.25582). Data on GitHub as doublelayer/SCAN. | Ion-pair interaction energy, dipole moment and optimised geometry | DLPNO-CCSD(T)/CBS on B2PLYP-D3/def2-TZVPP geometries | [BMIM], [BMPYR] and [C4PY] with [BF4], [PF6], [Tf2N], [TfO], [FSA], [DCA], [SCN], [C(CN)3], [TCB] and [Cl] |
| Zahn, MacFarlane, Izgorodina 2013, *PCCP* 15, 13664 (the IL-2013 set). [doi:10.1039/C3CP51682B](https://doi.org/10.1039/C3CP51682B) | Ion-pair interaction energies for 236 structures. Ion identities were assigned here by matching structures (102 structures, 36 of our pairs). | CCSD(T)/CBS, counterpoise corrected | Imidazolium and pyrrolidinium cations with [BF4], [PF6], [Tf2N], [DCA], [MeSO3], [TOS] and [Cl] |
| Stavarache, Calero, Vicent-Luna 2026, arXiv:2610.15190. Dataset [doi:10.5281/zenodo.20539350](https://doi.org/10.5281/zenodo.20539350) (CC0). | Ion and ion-pair geometries, binding energies and HOMO–LUMO gaps. 10 conformers per pair. | PBE-D3(BJ)/DZP, no counterpoise correction | 91 of our pairs (13 cations × 7 anions) |
| *PCCP* 2025, 27, 9068, "Accurate prediction of ionic liquid density-of-states…". [doi:10.1039/D5CP00214A](https://doi.org/10.1039/D5CP00214A), ESI Tables S12–S13. | Ionisation energy of each ion in the liquid, roughly equal to −E(HOMO). Validated against valence XPS. | B3LYP-D3(BJ)/6-311+G(d,p) with SMD continuum | 15 of our anions and 9 of our cations |
| Palomar et al. 2011, *Ind. Eng. Chem. Res.* 50, 3452. [doi:10.1021/ie101572m](https://doi.org/10.1021/ie101572m), SI Tables S3–S4. | CO₂ Henry constants from COSMO-RS, next to experiment | COSMO-RS. Ion-pair and separate-ion models, gas-phase geometries (see the paper for the parametrisation). | 28 pairs with imidazolium and pyrrolidinium cations and [BF4], [PF6], [Tf2N], [DCA], [TfO] and [MeSO4] |
| Doherty et al. 2017, *J. Chem. Theory Comput.* 13, 6131. [doi:10.1021/acs.jctc.7b00520](https://doi.org/10.1021/acs.jctc.7b00520), Tables 2 and 6. | Simulated density and viscosity at 298 K, next to experiment | OPLS-2009IL force field, with full charges and with charges scaled by 0.8 | 21 [EMIM], [BMIM] and [OMIM] pairs |
| Kerlé et al. 2017, *PCCP* 19, 1770. [doi:10.1039/C6CP06792A](https://doi.org/10.1039/C6CP06792A), Table 2. | Simulated CO₂ Henry constants at 300–500 K | MD with solvation free energies | [C₂/C₄/C₆/C₈mim][Tf2N] |

## How to use it for validation

1. **Compare like with like.** Absolute numbers only agree when the method is the same. Otherwise compare the **ranking and the trends**: Spearman correlation across ions, the anion ordering, and how a value changes with chain length.
2. **Ion-pair interaction energy** is the best check that the planned QM protocol is right.
   - Compute E_int = E(pair) − E(cation) − E(anion) at the pair geometry, counterpoise corrected, for the Karu and IL-2013 pairs.
   - Compare with the CCSD(T) values. IL-2013 reports a mean absolute deviation of about 3–5 kJ/mol for good dispersion-corrected functionals with diffuse basis sets.
   - Aim for a mean absolute deviation under 5 kJ/mol.
   - Karu and IL-2013 agree with each other: for [BMIM][BF4], Karu gives −375.2 kJ/mol and the two IL-2013 structures give −378.0 and −375.4 kJ/mol.
3. **Single-ion orbital energies:**
   - Check the ordering of our anion HOMOs against the ionisation energies in the 2025 PCCP SI. In that set [SCN]⁻ < [C(CN)3]⁻ < [DCA]⁻ < [NO3]⁻ < [MeSO3]⁻ < … < [Tf2N]⁻ < [FSA]⁻ < [TCB]⁻ < [BF4]⁻ < [PF6]⁻.
   - Gas-phase anion HOMOs depend strongly on diffuse functions. For [SCN]⁻ the HOMO moves from 0.40 eV with 6-31G(d,p) to 1.09 eV with 6-311+G(d,p) (2025 PCCP SI, Table S6).
4. **COSMO-RS:** reproduce the Palomar CO₂ Henry constants for 3–5 liquids, e.g. [BMIM][PF6], [BMIM][BF4] and [EMIM][Tf2N], before screening anything new.
5. **MD:**
   - Simulated densities should be within a few percent of Doherty's values and of experiment.
   - The CO₂ Henry constant at 300 K should be close to Kerlé's 28 ± 3 bar for [BMIM][Tf2N].
   - Viscosity is a harder test: errors of 10–25% are common.

## Lessons for our method choices

- **The current descriptors are r²SCAN-3c.** They come from a gas-phase optimisation and frequency calculation on one conformer per ion (`! r2SCAN-3c Opt Freq`; see `calculations/ions/*.out`). The r²SCAN-3c basis (mTZVPP) has no diffuse functions, so the gas-phase anion HOMO values are sensitive to the basis (the [SCN]⁻ example above). For the orbital descriptors, use a single point with diffuse functions (ωB97X-D3/def2-TZVPD, as planned).
- **Avoid B3LYP without a dispersion correction, and small basis sets such as def2-SVP, for ion-pair energies.** For B3LYP/def2-SVP the mean absolute deviation is 33.6 kJ/mol without counterpoise correction and 15.6 kJ/mol with it (IL-2013, Table 3). Adding D3 and a diffuse or triple-zeta basis brings it to about 4–7 kJ/mol. One stray input, `calculations/ions/BMIM.inp`, is a B3LYP/def2-SVP single point and does not match `BMIM.out`. It was not used for the descriptors.
- **Use a long-range-corrected, dispersion-corrected functional for anion–CO₂ complexes.** García et al. 2015 ([doi:10.1039/C5CP04283F](https://doi.org/10.1039/C5CP04283F)) recommend ωB97X-D. They found anion–CO₂ binding is the main term in IL–CO₂ binding. Their per-liquid values are not in an ESI, so they are not in the table.

## Gaps: no published computed reference found yet

- **Pairs.** The largest pairs with no reference are:
  - [EMIM][EtSO4] (262 points)
  - [HMIM][C(CN)3] (174)
  - [EMIM][DEPO4] (114)
  - [HMIM][BETA] (108)
  - [OMIM][C(CN)3] (100)
  - [C7MIM][C(CN)3] (96)
  - [C5MIM][Tf2N] (81)
  - [HMIM][SCN] (79)

  See `coverage_by_pair.csv` for the full list.
- **Ions with no reference of any kind:**
  - Cations:
    - pyridinium: [C4MPY], [HMPY], [PMPY], [HOEPY]
    - pyrrolidinium: [C7MPYR], [OMPYR], [COCMPYR]
    - functionalised imidazolium: [HOC3MIM], [NMIM], [aMIM], [C5MIM]
    - ammonium: [N1,1,1,2-OH], [N1,4,4,4], [N2,1,1,3], [N6,1,1,1], [deme], [thma]
    - phosphonium: [P1,4,4,4], [P2,4,4,4]
  - Anions: [EtSO4], [DEPO4], [DMPO4], [MDEGSO4], [C3F7CO2], [NfO], [methide]
- **For these, trust comes from internal checks instead:**
  - no imaginary frequencies
  - the correct charge and multiplicity
  - agreement between r²SCAN-3c and ωB97X-D3 values
  - smooth trends within each family (alkyl chain length, sulfate versus sulfonate)

## Sources worth fetching next

- **Ong & Ceder 2010, *Electrochim. Acta* 55, 3804:** B3LYP ionisation and electron affinity values for more than 200 ions. MIT DSpace was blocking access on the day of this search.
- **Sudha & Khanna 2009:** MP2 anion–CO₂ binding energies for [BF4], [PF6], [Tf2N], [TfO] and [NO3].
- **Seymour et al. 2021, *PCCP* (probably doi:10.1039/D1CP02441H; not yet confirmed):** experimental ionisation energies for 60 ionic liquids. This is the experimental anchor for HOMO energies.
- **Seeger & Izgorodina 2022, *J. Comput. Chem.*:** DLPNO-CCSD(T) benchmarks for protic and aprotic ion pairs. Useful for the ammonium cations.
- **The IL195×8 set** (Carter-Fenk & Herbert 2022, *Mol. Phys.* 121, e2055504): IL-2013 pairs at 8 separations.
