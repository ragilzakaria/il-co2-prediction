# Ion-pair pipeline test (laptop, 26–27 Nov 2026)

Three systems through all three stages of `../pipeline.py` on the MacBook (7 runs side by side):
[EMIM][BF4] and [EMIM][DCA], which have CCSD(T) reference values, and [BF4]⁻·CO₂.

| Stage | Runs | Wall time |
|---|---|---|
| A: GFN2-xTB optimisation of every start | 59 | 49 min |
| B: r²SCAN-3c optimisation of 3 picks per system | 9 | 5.9 h (laptop shared with other runs) |
| C: ωB97X-D3/def2-TZVPD single points | 9 | 65 min |

Files: `selected_A.tsv`, `selected_B.tsv`, `results.tsv`, `test.log` as written by the pipeline;
`geometries/` the stage-B structures; `singlepoints/` the stage-C outputs.

**These files come from the pipeline before the duplicate check was fixed.** The old check compared
interatomic distances atom by atom, so the same structure with swapped equivalent atoms (the F of
[BF4]⁻, methyl H) looked 0.6–1.0 Å different. All three stage-B picks of [EMIM][BF4] (starts 8, 0, 2)
were one structure. The fixed check (sorted distances) picks starts 8, 10 and 4 instead; starts 10
and 4 both optimise at r²SCAN-3c to a structure 7.1 kJ/mol lower than start 8, the same as the
CCSD(T) structure of Zahn et al. re-optimised at r²SCAN-3c (within 0.03 kJ/mol).

Interaction energy, kJ/mol:

| System | Pipeline, old check | Pipeline, fixed check | ωB97X-D3 on the CCSD(T) geometry | CCSD(T)/CBS |
|---|---|---|---|---|
| [EMIM][BF4] | −371.7 | **−378.2** | −376.6 | −378.7 (2 structures, to −372.7) |
| [EMIM][DCA] | −363.0 | −363.0 (best pick unchanged) | −362.4 | −368.2 (6 structures, to −361.0) |
| [BF4]⁻·CO₂ | −28.2 | −28.2 | — | no reference in `data/literature_reference/` |

The fixed-check [EMIM][BF4] value is the ωB97X-D3 single point on the r²SCAN-3c re-optimised CCSD(T)
structure (scratch run, same settings as stage C). Using r²SCAN-3c geometries instead of the
published ones changes the interaction energy by 1.6 kJ/mol for [EMIM][BF4]. For [EMIM][DCA] the
pipeline structure differs from the 0.6 kJ/mol lower CCSD(T) structure (r²SCAN-3c energies), and its
interaction energy lies 0.6 kJ/mol from the published-geometry value and inside the CCSD(T) range of
the six structures; the remaining 5 kJ/mol to the best CCSD(T) value is the method error measured in
the benchmark (5.8 kJ/mol for this pair).
