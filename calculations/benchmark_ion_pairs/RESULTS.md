# Ion-pair benchmark results (roadmap step 2.1)

66 published ion-pair structures (30 from Karu 2018, DLPNO-CCSD(T)/CBS; 36 from IL-2013, CCSD(T)/CBS),
single points on the published geometries with ORCA 6.1.1 on Mahameru (job 575850, 26 Nov 2026,
528 of 528 runs terminated normally). Interaction energy E_int = E(pair) − E(cation) − E(anion);
the references range from −420 to −298 kJ/mol.

Produced by `python scripts/analyse_ionpair_benchmark.py` (per-pair values in `results.csv`).

## Deviation from the reference (method − reference, kJ/mol)

| Set | Method | MD | MAD | MAXD | MAD % | Rank corr. |
|---|---|---|---|---|---|---|
| all 66 | r²SCAN-3c | −5.6 | 6.0 | 33.0 | 1.6 | 0.988 |
| | ωB97X-D3/def2-TZVPD | +1.1 | **3.7** | 10.4 | 1.0 | 0.986 |
| | ωB97X-D3/def2-TZVPD + CP | +2.7 | 4.0 | 11.2 | 1.1 | 0.986 |
| without [Cl]⁻ (57) | r²SCAN-3c | −3.7 | **4.1** | 8.3 | 1.1 | 0.989 |
| | ωB97X-D3/def2-TZVPD | 0.0 | **3.0** | 9.5 | 0.9 | 0.981 |
| | ωB97X-D3/def2-TZVPD + CP | +1.7 | 3.2 | 9.8 | 0.9 | 0.980 |
| [Cl]⁻ only (9) | r²SCAN-3c | −17.9 | 17.9 | 33.0 | 4.4 | 0.817 |
| | ωB97X-D3/def2-TZVPD | +8.2 | 8.2 | 10.4 | 2.0 | 0.983 |

Target from the roadmap: MAD < 5 kJ/mol.

## Cost (one core each, serial)

| Method | Whole benchmark | One ion pair, median (max) |
|---|---|---|
| r²SCAN-3c | 3.4 core-h | 1.8 min (3.9 min) |
| ωB97X-D3/def2-TZVPD (incl. CP runs) | 121.5 core-h | 32.5 min (91 min) |

## Findings

1. ωB97X-D3/def2-TZVPD meets the target for every anion family (MAD 3.7 kJ/mol, about 1 %).
2. Counterpoise correction does not improve it (4.0 against 3.7 kJ/mol): def2-TZVPD with its diffuse
   functions has little basis-set superposition error here. CP runs can be dropped, which saves
   two of the five calculations per pair.
3. r²SCAN-3c overbinds by a nearly constant 4–6 kJ/mol and meets the target for all anions except
   chloride, which it overbinds by 12–33 kJ/mol (compact anion, delocalisation error, no diffuse
   functions). No IL in the training data has chloride.
4. All three methods order the pairs almost exactly as the reference does (rank correlation 0.98–0.99),
   so for ML features the choice matters less than for absolute energies.
5. r²SCAN-3c is about 35 times cheaper.

## Limits

- Single points on the published (high-level) geometries; in production the geometries will come
  from r²SCAN-3c, which adds a small error not tested here.
- The benchmark anions are [BF4], [PF6], [Tf2N], [DCA], [SCN], [TCB], [C(CN)3], [FSA], [MeSO3],
  [TfO], [TOS] and [Cl]. Training anions not covered: the sulfates, phosphates, carboxylates,
  [NO3], [NfO], [BETA] and [methide].
