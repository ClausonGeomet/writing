# Research findings: Part 2A

The current freeze passes its prospective numerical gates: 5,400 baseline fits have status `ok`; constrained positivity and closure errors are below `1e-8`; forward replay is below `1e-10`; easy-case maximum absolute phase error across methods is `4.17e-8`, below `1e-6`. These are synthetic correctness checks, not field validation. The historical v1 13 SLSQP failures and v2 incomplete 55/56 snapshot remain visible; current v3 verifies 55/55.

Mean test errors below pool three 20-specimen replicate tests. Values are whole-sample fractions and rounded from [`results/summary.csv`](results/summary.csv). The comparator was fixed before B results: `cls_eq` except `cls_le` for `missing_phase`.

| Scenario | Assay RMSE | Mineral RMSE | Interpretation |
| --- | ---: | ---: | --- |
| easy | <0.000001 | <0.000001 | Full-rank exact recovery |
| ambiguous | <0.000001 | 0.083206 | Exact chemistry leaves phase split unresolved |
| noisy | 0.004316 | 0.008752 | Error propagates into mineral fit |
| missing_phase | 0.064301 | 0.122035 | Omitted mineral leaves large unexplained assay |
| wrong_chemistry | 0.005288 | 0.014450 | Fixed-response misspecification matters |
| infeasible | 0.045520 | 0.081103 | No exact closed fit exists |

In the ambiguous three-phase response, chemical rank and rank after adding closure are both 2. Exact LP Anatase bounds span widths 0.213–0.617 across generated specimens. Swapping Anatase and Rutile truth in 180 paired rows leaves each assay exactly unchanged while phase-specific truth differs. This is stronger than matching first moments or marginal histograms: the entire empirical assay rows coincide. Mineral accuracy cannot be read from assay residual alone.

The [ILR summary](results/summary.csv) contains 270 rows, of which 126 have finite `ilr_rms` and 144 have zero eligible rows/NaN. Across 900 row instances per scenario, eligible counts are ambiguous 540, easy 889, infeasible 334, missing phase 0, noisy 268 and wrong chemistry 360. Support restrictions are part of the result, not missing output. No uncertainty intervals were produced. The [dataset screen](results/dataset_screen.csv) admits no independent paired field cohort; previous Fortescue 14/5/4 exposure and failure are historical context, not a new test.
