# Prior-guided EMC experiment plan

**Prospective only. Checkpoint NOT ACCEPTED; campaigns NOT RUN; external validity unestablished.** Current evidence consists of a proof and tiny exact arithmetic checks: independent tester 171/171 exact assertions passed; five builder unit tests passed. The aid handles up to four phases and a diagonal score; it is not a covariance estimator, optimizer, prior validation or Monte Carlo result.

## Comparator matrix

| Arm | Estimand/method | Present status |
| --- | --- | --- |
| A0 | Prior-free constrained least squares on declared noisy feasible Q | NOT RUN |
| A1 | External-prior-only posterior predictive mean under squared mineral loss, respecting support/closure; no target assay | NOT RUN |
| A2 | Assay plus external prior, same Q and likelihood; primary comparator | NOT RUN |
| B0 | Frozen prior-free greedy sequential rule, historical reference only | Not rerun here |
| B1 | Prior-guided fixed order, free per-sample fractions and common J on P | NOT RUN; redundant order parameterization |
| B2 | Exhaustive order with common J on P | NOT RUN; redundant order parameterization |
| B3 | Conditional GA over justified discrete rule/support and fractions | NOT RUN / NOT JUSTIFIED by current model |

A0/A1/A2 isolate assay and prior contribution; B1/B2/B3 at fixed prior isolate search. Primary prospective comparison is B3 versus A2 **only if** a scientifically meaningful discrete decision remains after review. Give both arms the same calibration, target evidence, chemistry uncertainty, support, closure, tuning budget and final-state objective. Also report A2-on-P as a matched auxiliary comparator for the decoder's one-sided resource cap; do not silently change primary A2's noisy Q. If a support mixture is used, A2 receives it too, with enumeration or an appropriate mixed/discrete baseline. Compare finite search to exhaustive search in the tiny three/four-phase setting and to simple random/multistart methods if discrete search is justified. A lower training objective is not necessarily a lower held-out mineral error.

## Analytical gate and future data design

Current quartz–anatase–rutile witnesses cover exact nullspace, strong/weak accurate/wrong prior centers, greedy versus free fractions, noisy P⊂Q, signed negative assays, absent phase and incomplete inventory. They are derivations/checks, not empirical prior performance. The theorem makes an order-plus-free-fractions B3 a no-go. Reopen only after an independently motivated discrete physical rule or cohort parameter restriction is specified, and demonstrate that it changes the **common** feasible model rather than only an optimizer path. No production GA or main campaign should start at this checkpoint.

After user acceptance, design a separate calibrated study using simCoDaR revision `d9a0b27676a0cb3a4acb300a3f6b2ebc845d205d`, pinned for **future** use only. Predeclare four-phase cohorts, DGP, assay/noise and chemistry laws, domain/support distribution, independent calibration groups, validation/test groups, and a seed/version manifest before generation. Include accurate and biased prior means, weak/strong precision, calibration sizes 10/30/100 independent units (sensitivity levels, not adequacy claims), domain misclassification, QXRD error and censoring, uncertain C, absent and missing phases, and signed negative assays. Separate oracle-prior diagnostics from finite-calibration analysis. Keep all unsupported cases and solver failures; define their primary-loss treatment or restricted estimand before testing. Do not compare different surviving subsets.

## Frozen future decision rule to preregister

Primary paired unit is an independently generated cohort, not a row or optimizer restart. Define g=RMSE_m(A2)−RMSE_m(B3) and d=RMSE_a(B3)−RMSE_a(A2), both on original whole-sample fractions. The historical .01 mineral gain and .001 maximum assay degradation are **prospective practical margins**, not achieved results. Strong confirmatory success would require lower one-sided 95% bound for mean g>.01 **and** upper one-sided 95% bound for mean d<.001 in every declared primary scenario. Predeclare multiplicity adjustment if making an any-scenario claim; never select the best test scenario. Report paired confidence intervals, MCSE=s/√R, failures and calibration uncertainty. A2 versus A0 prior benefit is a separate secondary contrast.

Propose at least R=200 independent cohort pairs per scenario, then increase before confirmatory execution using separate pilot/conservative variance evidence. For desired two-sided 95% half-width .002 on g, plan R≥ceil((1.96 s_g/.002)²); for half-width .0002 on d, R≥ceil((1.96 s_d/.0002)²). R=200 alone does not guarantee precision. If B3 is stochastic, predeclare 20 search seeds within each cohort, report mean loss over seeds as primary and spread/failure as secondary; those seeds are not independent geological replicates. Interval coverage, if later attempted, needs at least 2000 independent designated replicates per scenario (MCSE near .95 ≈.0049), explicit binomial intervals and calibration uncertainty. No coverage claim follows from three old seeds or the exact checks. Freeze all tuning on calibration/validation only. No repeated looking and favorable stopping without a separate rule.

Next authorized checkpoint actions are independent scientific review and a stacked draft PR, followed by stop awaiting user acceptance. No simulation activation, package changes, article render or publication is part of this stage.
