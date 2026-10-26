# R2 saved-results comparison audit

This is a read-only analysis of the saved Part A `summary.csv` and Part B `metrics.csv` **test** rows. No fits, search runs, seed selection, or baseline rewrites were performed. Source evidence: [Part A summary](../emc-part2a/results/summary.csv) and the saved Part B `metrics.csv` test rows, whose original fields are retained in the public [provenance derivative](results/metrics_provenance.csv).

## Complete machine-readable comparisons

- `results/all_seed_comparisons.csv`: 1,296 paired rows. It includes every scenario, replicate 0/1/2, objective, method, and saved search seed. Families: B versus predeclared A 360; search versus fixed with the **same objective** 324; distribution versus row with the **same method and seed** 180; GA efficiency versus exhaustive, greedy, beam, and **equal-seed random** 432. Each row carries candidate/reference scores, mineral and assay errors, actual evaluations, measured fit times, signed changes, orders, and criterion booleans.
- `results/all_seed_medians.csv`: 648 groups at scenario × replicate × objective × method × reference granularity, including median signed change across all available saved seeds and `seed_count`. For random and GA each median pools exactly seeds 101/102/103; deterministic methods have one saved seed 0. The 648 groups comprise 216 B/A, 180 search/fixed, 108 distribution/row, and 144 GA efficiency groups.
- `results/all_seed_benefit.csv`: all 120 scenario × objective × method × seed B/A evaluations, each with three replicate pairs, mineral-gain count, mean assay change, median mineral change, and the prospective benefit verdict.
- `results/all_seed_ga_efficiency.csv`: all 144 scenario × objective × comparator × GA-seed groups, each with exactly three matched replicate pairs, component counts, joint count, medians, and the unchanged ≥2/3 verdict.

All signed deltas are **candidate minus reference** in original whole-sample mass-fraction units for RMSE. Negative mineral or assay delta means lower error. Search scores are training scores; `delta_train_score` compares only the same objective. A baseline has no order or search score. Fit time and evaluations are per fitted cohort and are repeated on its test row; they are not per-test-prediction costs. All 360 B test rows and all 90 A test-summary rows are present. Comparison keys are unique; there is no best-seed filter.

## Scientific findings

**Search versus fixed order.** Search has substantial held-out mineral and assay gains over the lexical fixed order in easy, missing-phase, noisy, and wrong-chemistry cases for exhaustive, beam, random, and GA. In the easy case, median search-minus-fixed mineral RMSE is about −0.173455 and assay RMSE about −0.073394; in missing phase, about −0.048894 and −0.009680. Greedy selected the fixed result in all saved cases. Ambiguous orders give identical held-out errors. Infeasible distribution search reduces mineral error versus fixed by about 0.123 median, while *increasing* assay error by about 0.0093; row search lowers both errors in this case. These are gains over a weak fixed allocation order, not proof that B beats the predeclared constrained baseline.

**Distribution versus row objective.** For the same search and seed, most saved held-out errors are identical. In infeasible cohorts, distribution search lowers mineral RMSE by roughly 0.03–0.054 per replicate, but raises assay RMSE by about 0.034. That assay increase exceeds the unchanged **≤0.001** mean-increase rule. Wrong chemistry has smaller method-specific mineral changes, including a beam regression. Fixed and greedy are identical between objectives in held-out error. The detailed CSV retains every nonzero and near-zero difference, including all three random and GA seeds.

**B versus predeclared A.** The reference is `cls_eq` in every scenario except `missing_phase`, whose reference is `cls_le`. Across all 120 B method/objective/seed groups, **0/120** meet the prospective rule: mineral RMSE reduction **≥0.01** in **≥2/3** matched replicates *and* mean assay RMSE increase **≤0.001**. No B/A pair reaches the mineral threshold. Fixed/greedy B is often much worse than A; optimized B is near A in easy and missing phase but still lacks the required mineral gain. Infeasible distribution can slightly lower mineral RMSE relative to A in some replicates, but raises assay RMSE substantially. This corrects the earlier broad conclusion based on a held-out-selected A comparator. No scientific acceptance threshold was changed.

**GA efficiency.** For each scenario × replicate × objective × GA seed, the comparator is exhaustive, greedy, or beam seed 0, or random at the *same* seed. Each objective/comparator/GA-seed denominator is **18 matched scenario–replicate pairs**; each scenario/objective/comparator/seed denominator is **3 matched replicates**. The joint criterion is training score no worse than comparator by more than **1e-10**, **strictly fewer actual evaluations**, and **strictly lower measured fit time** in **≥2/3** matched replicates. Counts below are joint successes out of 18, followed by objective-only successes out of 18. Both objectives and all seeds are shown.

| Objective | Comparator | Seed 101 joint; score | Seed 102 joint; score | Seed 103 joint; score |
| --- | --- | ---: | ---: | ---: |
| row | exhaustive | 15; 18 | 12; 15 | 18; 18 |
| row | greedy | 0; 18 | 0; 18 | 0; 18 |
| row | beam | 0; 18 | 9; 15 | 9; 18 |
| row | equal-seed random | 15; 18 | 12; 15 | 18; 18 |
| distribution | exhaustive | 14; 17 | 11; 14 | 17; 17 |
| distribution | greedy | 0; 18 | 0; 18 | 0; 18 |
| distribution | beam | 0; 17 | 8; 14 | 6; 17 |
| distribution | equal-seed random | 13; 17 | 11; 14 | 17; 17 |

At the scenario/objective/comparator/seed level, the joint **≥2/3** rule passes for row exhaustive **15/18** groups, greedy **0/18**, beam **6/18**, equal-seed random **15/18**; for distribution, exhaustive **15/18**, greedy **0/18**, beam **5/18**, equal-seed random **15/18**. Each of the 18 groups per objective and comparator has three replicate pairs and its own recorded verdict in `ALL_SEED_GA_EFFICIENCY.csv`.

The prior seed-101-only statement of **0/36 GA-versus-beam joint wins** is not the full saved evidence: across all seeds GA has **18/54 row** and **14/54 distribution** joint wins versus beam. Against greedy it has **0/54** for each objective; against exhaustive **45/54 row** and **42/54 distribution**; against equal-seed random **45/54 row** and **41/54 distribution**. These pooled counts describe the saved pairs, not a prospective ≥2/3 verdict for every three-replicate scenario. GA fit times range from about 0.0065 to 0.0262 seconds (median about 0.0152 seconds); seed-dependent ranking at this scale does not support a confident general speed advantage. Per-scenario three-replicate outcomes and seed medians are reconstructible directly from the complete CSVs.

## Comparison conventions and limits

The exact prospective thresholds used here equal `test-spec.md`: mineral gain **≥0.01**, mean assay increase **≤0.001**, matched benefit count **≥2/3**, and GA objective tolerance **1e-10** with strict evaluation/runtime improvement and **≥2/3** pairs. No tolerance was relaxed. No interval estimates were produced, so coverage, size, power, and SE-ratio are not applicable. The fixed-order gains are algorithmic order gains within B; A superiority is assessed separately with the predeclared baseline. All claims are confined to these saved synthetic test cohorts and measured fit times.
