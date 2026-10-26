# Research findings: Part 2B

All 360 planned fits completed; 1,080 split-level B metric rows report no unexpected failure or unsupported sample in these generated cohorts. Independent checks replay saved predictions, verify the exhaustive objective oracle and unchanged train/validation/test IDs, and pass the numerical gates. These are synthetic feasibility findings. The original `metrics.csv` is currently in the main checkout, while the committed [provenance derivative](results/metrics_provenance.csv) preserves its 1,080 numeric/timing rows and adds units and dependency hash. The [sidecar](results/metrics_provenance.json) records the original source SHA-256 `e24079a4cfeb13e42ecc2ba900284ed85b68800796b406ab3c63261437749245` and derivative SHA-256 `d969e506275ffafa3906360dcc63c8aa0656401db4a3f5b6e84d31efe403135a`.

Means below use each of three 20-row test splits, with exhaustive order search for the B columns. All RMSEs are whole-sample fractions, rounded to six decimals. A is the **predeclared** `cls_eq` comparator except `missing_phase` (`cls_le`). B row and distribution objectives use the same allocator.

| Scenario | A assay | A mineral | B row assay | B row mineral | B distribution assay | B distribution mineral |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| easy | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 |
| ambiguous | 0.000000 | 0.083206 | 0.000000 | 0.086836 | 0.000000 | 0.086836 |
| noisy | 0.004316 | 0.008752 | 0.004978 | 0.009772 | 0.004978 | 0.009772 |
| missing_phase | 0.064301 | 0.122035 | 0.064301 | 0.122035 | 0.064301 | 0.122035 |
| wrong_chemistry | 0.005288 | 0.014450 | 0.012803 | 0.033874 | 0.012803 | 0.032428 |
| infeasible | 0.045520 | 0.081103 | 0.063043 | 0.116787 | 0.097596 | 0.080003 |

No method/scenario achieved mineral RMSE reduction ≥`0.01` against the predeclared A comparator in at least two of three matched replicates. In the infeasible case, moving from B row to distribution objective reduced mineral RMSE by about `0.032–0.055` per replicate (method dependent) but raised assay RMSE by about `0.034–0.036`, far above the `0.001` allowance. It is therefore an assay/mineral tradeoff, not a passing benefit. The ambiguous case remains mineral-inaccurate despite exact assays; the 180 paired Anatase/Rutile swaps have identical observed assay rows. Matching even a complete assay law would not resolve that chemical null direction.

At seed 101, GA matched the exhaustive training objective within `1e-10` in 18/18 row pairs and 17/18 distribution pairs; the one exception was wrong chemistry replicate 0 by `0.000408165`. Joint score/fewer evaluations/lower-time counts versus exhaustive were 15/18 and 14/18. Versus **beam**, GA used fewer evaluations in 0/18 pairs for either objective, so the predeclared GA efficiency claim fails. Search times are milliseconds (GA median about 0.015 s, beam about 0.018 s); these measurements do not support a durable speed claim. Exhaustive is a finite training-objective oracle here, not a universal mineral optimum.

Fresh external validation remains pending. Earlier Fortescue 14/5/4 work was already exposed and its failure is not a newly replicated finding. The numerical record retains early Part A failures, historical v2's missing snapshot reference (55/56) and the complete current freeze (55/55). The findings support a narrower conclusion: chemical information and independently verified phase data matter more than order-search complexity for this pilot.
