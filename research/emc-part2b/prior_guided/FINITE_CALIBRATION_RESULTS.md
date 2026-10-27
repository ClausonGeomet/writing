# Bounded finite-calibration QAR pilot: results

**27 October 2026 · numerical audit PASS; final internal independent agent review PASS; STOP after this pilot.** The initial root attempt is **INCOMPLETE**. The completed evidence is the retained `correction-001` attempt. These descriptive synthetic results are neither a confirmatory test nor evidence of field performance or general empirical superiority.

## Design and estimators

The [frozen plan](FINITE_CALIBRATION_PLAN.md) defines 10 cells, each with 50 independently generated calibration/evaluation cohorts and 40 held-out evaluation specimens per cohort (500 fits; 20,000 evaluation specimens). Primary cells 1–6 use calibration sizes 10, 30, 100 for fixed or variable Quartz. Stress cells 7–10 use n=30 only. Calibration and evaluation populations are generated independently; the mismatch cell changes the calibration source center. The known synthetic laboratory covariance is supplied with signed, unclosed QXRD-like calibration measurements. No known DGP center, latent truth, held-out assay, or evaluation truth enters prior fitting.

A0 minimizes assay loss on the common closed feasible set, with an explicit minimum-Euclidean-norm tie for identical Ti chemistry. A1 is the projected **calibrated center**, without a target assay; it is not the earlier proposed posterior predictive mean. A2 globally solves the convex assay-plus-learned-quadratic-prior objective on the same feasible set. The fitted tangent covariance subtracts known measurement covariance, applies an eigenvalue floor of 0.0001 and 0.25 shrinkage, and adds an approximate mean-estimation covariance to the plug-in penalty. The mean projection, floor, and shrinkage can bias the result. Covariance-parameter uncertainty is **not integrated**. The chemical response has identical Anatase/Rutile columns, leaving their contrast exactly unobserved by bulk assay.

## Held-out cohort metrics

Each entry below is the mean of 50 cohort metrics, rounded for display. Mass RMSE is over the three phases; assay RMSE uses two channels. `Obs` is reconstruction versus noisy assay, while `Latent` is reconstruction versus noiseless assay. `Ti` is total Ti-bearing phase mass RMSE. Phase columns give Quartz/Anatase/Rutile RMSE and signed bias, in that order. The [full precision cohort CSV](finite_calibration/attempts/correction-001/results/cohort_metrics.csv) is authoritative.

| Cell | Arm | Mass | Assay Obs | Assay Latent | Phase RMSE Q/A/R | Phase bias Q/A/R | Ti |
| ---: | --- | ---: | ---: | ---: | --- | --- | ---: |
| 1 | A0 | 0.08375 | 0.00994 | 1.71e-16 | 0 / 0.10257 / 0.10257 | 0 / 0.09580 / -0.09580 | 2.41e-16 |
| 1 | A1 | 0.03115 | 0.00994 | 2.84e-16 | 2.91e-16 / 0.03815 / 0.03815 | 2.91e-16 / 0.00083 / -0.00083 | 2.6e-16 |
| 1 | A2 | 0.03115 | 0.00994 | 1.35e-15 | 1.65e-15 / 0.03815 / 0.03815 | 8.66e-18 / 0.00083 / -0.00083 | 9.2e-16 |
| 2 | A0 | 0.08361 | 0.00975 | 1.63e-16 | 0 / 0.10239 / 0.10239 | 0 / 0.09582 / -0.09582 | 2.28e-16 |
| 2 | A1 | 0.02994 | 0.00975 | 2.81e-16 | 3e-16 / 0.03666 / 0.03666 | 3e-16 / -0.00045 / 0.00045 | 2.53e-16 |
| 2 | A2 | 0.02994 | 0.00975 | 1.35e-15 | 1.65e-15 / 0.03666 / 0.03666 | 3.82e-17 / -0.00045 / 0.00045 | 9.1e-16 |
| 3 | A0 | 0.08413 | 0.00977 | 1.72e-16 | 0 / 0.10304 / 0.10304 | 0 / 0.09615 / -0.09615 | 2.43e-16 |
| 3 | A1 | 0.03053 | 0.00977 | 2.9e-16 | 3.04e-16 / 0.03739 / 0.03739 | 3.04e-16 / -0.00018 / 0.00018 | 2.65e-16 |
| 3 | A2 | 0.03053 | 0.00977 | 1.07e-15 | 1.3e-15 / 0.03739 / 0.03739 | 4.72e-18 / -0.00018 / 0.00018 | 7.43e-16 |
| 4 | A0 | 0.09106 | 0.00715 | 0.00704 | 0.00704 / 0.11140 / 0.11141 | -7.74e-05 / 0.09961 / -0.09953 | 0.00704 |
| 4 | A1 | 0.08663 | 0.10625 | 0.10588 | 0.10588 / 0.04264 / 0.09697 | 1.68e-05 / -0.00110 / 0.00108 | 0.10588 |
| 4 | A2 | 0.03441 | 0.00718 | 0.00704 | 0.00704 / 0.04161 / 0.04207 | -6.28e-05 / -0.00049 / 0.00055 | 0.00704 |
| 5 | A0 | 0.09227 | 0.00713 | 0.00729 | 0.00729 / 0.11303 / 0.11275 | -0.00021 / 0.10098 / -0.10076 | 0.00729 |
| 5 | A1 | 0.08454 | 0.10398 | 0.10344 | 0.10344 / 0.04096 / 0.09489 | 0.00381 / -0.00046 / -0.00335 | 0.10344 |
| 5 | A2 | 0.03246 | 0.00715 | 0.00727 | 0.00727 / 0.03932 / 0.03950 | -0.00019 / 0.00062 / -0.00043 | 0.00727 |
| 6 | A0 | 0.09169 | 0.00700 | 0.00686 | 0.00686 / 0.11201 / 0.11238 | 0.00030 / 0.09921 / -0.09951 | 0.00686 |
| 6 | A1 | 0.08585 | 0.10512 | 0.10486 | 0.10486 / 0.04153 / 0.09659 | 0.00176 / -0.00113 / -0.00063 | 0.10486 |
| 6 | A2 | 0.03242 | 0.00702 | 0.00687 | 0.00687 / 0.03907 / 0.03972 | 0.00031 / -0.00070 / 0.00040 | 0.00687 |
| 7 | A0 | 0.09338 | 0.00684 | 0.00696 | 0.00696 / 0.11419 / 0.11434 | 0.00015 / 0.10154 / -0.10169 | 0.00696 |
| 7 | A1 | 0.08613 | 0.10568 | 0.10508 | 0.10508 / 0.04198 / 0.09682 | 0.00464 / -0.00107 / -0.00357 | 0.10508 |
| 7 | A2 | 0.03271 | 0.00684 | 0.00696 | 0.00696 / 0.03951 / 0.03997 | 0.00015 / -0.00024 / 9.73e-05 | 0.00696 |
| 8 | A0 | 0.09208 | 0.00706 | 0.00705 | 0.00705 / 0.11252 / 0.11280 | 0.00025 / 0.10075 / -0.10099 | 0.00705 |
| 8 | A1 | 0.18499 | 0.10588 | 0.10536 | 0.10536 / 0.20446 / 0.22229 | 0.00183 / 0.20015 / -0.20198 | 0.10536 |
| 8 | A2 | 0.17443 | 0.00708 | 0.00704 | 0.00704 / 0.21339 / 0.21375 | 0.00026 / 0.20113 / -0.20139 | 0.00704 |
| 9 | A0 | 0.09202 | 0.00708 | 0.00688 | 0.00688 / 0.11257 / 0.11261 | -2.76e-05 / 0.10033 / -0.10031 | 0.00688 |
| 9 | A1 | 0.08702 | 0.10709 | 0.10666 | 0.10666 / 0.04300 / 0.09702 | 0.00025 / 0.00155 / -0.00180 | 0.10666 |
| 9 | A2 | 0.03418 | 0.00711 | 0.00686 | 0.00686 / 0.04125 / 0.04189 | -3.42e-05 / 0.00194 / -0.00191 | 0.00686 |
| 10 | A0 | 0.09267 | 0.00711 | 0.00698 | 0.00698 / 0.11321 / 0.11355 | -0.00019 / -0.10555 / 0.10575 | 0.00698 |
| 10 | A1 | 0.11353 | 0.11660 | 0.11617 | 0.11617 / 0.11321 / 0.10977 | 0.05086 / -0.10555 / 0.05469 | 0.11617 |
| 10 | A2 | 0.09250 | 0.00715 | 0.00698 | 0.00698 / 0.11321 / 0.11315 | 0.00013 / -0.10555 / 0.10542 | 0.00698 |

## Paired cohort contrasts

Entries are mean A2 minus comparator ± descriptive SE, computed from the 50 independent *cohort* differences in each cell (planned=complete=50 for every entry). Negative means A2 has lower error. These are not specimen SEs, confidence intervals, or hypothesis tests. [Full precision paired summary](finite_calibration/attempts/correction-001/results/paired_summary.csv).

| Cell | Scenario | A2−A0 mass | A2−A1 mass | A2−A0 assay Obs | A2−A1 assay Obs | A2−A0 assay Latent | A2−A1 assay Latent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | fixed, n=10 | -0.05260 ± 0.00081 | 2.15e-15 ± 2.8e-15 | -4.14e-16 ± 7.27e-17 | -4.14e-16 ± 7.25e-17 | 1.18e-15 ± 1.89e-16 | 1.06e-15 ± 1.9e-16 |
| 2 | fixed, n=30 | -0.05367 ± 0.00076 | -1.02e-16 ± 2.57e-15 | -3.97e-16 ± 6.02e-17 | -3.95e-16 ± 5.99e-17 | 1.18e-15 ± 1.62e-16 | 1.07e-15 ± 1.62e-16 |
| 3 | fixed, n=100 | -0.05361 ± 0.00080 | 9.87e-16 ± 3.35e-15 | -3.1e-16 ± 5.08e-17 | -3.1e-16 ± 5.1e-17 | 9e-16 ± 1.25e-16 | 7.83e-16 ± 1.25e-16 |
| 4 | variable, n=10 | -0.05664 ± 0.00101 | -0.05222 ± 0.00145 | 2.69e-05 ± 4.19e-06 | -0.09907 ± 0.00190 | -3.69e-06 ± 1.42e-05 | -0.09885 ± 0.00193 |
| 5 | variable, n=30 | -0.05981 ± 0.00107 | -0.05208 ± 0.00115 | 2.39e-05 ± 2.45e-06 | -0.09683 ± 0.00158 | -2.16e-05 ± 1.46e-05 | -0.09617 ± 0.00158 |
| 6 | variable, n=100 | -0.05927 ± 0.00121 | -0.05342 ± 0.00105 | 2.24e-05 ± 1.04e-06 | -0.09810 ± 0.00143 | 5.64e-06 ± 1.27e-05 | -0.09800 ± 0.00145 |
| 7 | weak, n=30 | -0.06068 ± 0.00103 | -0.05342 ± 0.00106 | 2.22e-09 ± 1.57e-10 | -0.09884 ± 0.00142 | -3.78e-07 ± 1.04e-07 | -0.09811 ± 0.00145 |
| 8 | mismatch, n=30 | 0.08235 ± 0.00099 | -0.01057 ± 0.00076 | 2.56e-05 ± 2.08e-06 | -0.09879 ± 0.00135 | -7.02e-06 ± 1.64e-05 | -0.09832 ± 0.00139 |
| 9 | noisy_qxrd, n=30 | -0.05783 ± 0.00103 | -0.05284 ± 0.00130 | 2.89e-05 ± 4.66e-06 | -0.09999 ± 0.00170 | -1.25e-05 ± 1.42e-05 | -0.09980 ± 0.00171 |
| 10 | wrong_occurrence, n=30 | -0.00016 ± 9.61e-06 | -0.02103 ± 0.00099 | 4.25e-05 ± 3.53e-06 | -0.10944 ± 0.00213 | -3.21e-06 ± 1.7e-05 | -0.10919 ± 0.00209 |

The correct-model synthetic design gives the calibrated, deliberately skewed Anatase/Rutile center an advantage over A0's equal-Ti minimum-norm tie in cells 1–6. The advantage does not remove irreducible specimen-level variation along the exact assay null direction. In fixed-Quartz cells 1–3, A1 and A2 coincide within 1.37×10⁻¹² in mass: the assay cannot learn the Ti split. In variable-Quartz primary cells 4–6, A2 improves on A1 mainly by fitting total Ti/Quartz, not by identifying the Anatase/Rutile contrast. The observed-assay and latent-assay columns must remain distinct; A0's close noisy-assay fit is an objective consequence, not a recovery guarantee.

## Validation, provenance, and limits

The tester audit (`runs/REQ-20261027-finite-calibration/audit.md` in the workflow record) records nine unit checks and 84 independent analytic checks at first execution, then 56/56 result and replay checks after correction. All 500 fits and 60,000 predictions completed; failures were 0/500 fits and 0/20,000 predictions per arm, including 0/50 fits and 0/2,000 predictions per arm in every cell. Maximum constraint violation was 1.175×10⁻¹³ (gate 10⁻⁹), maximum scaled KKT residual 4.012×10⁻¹⁴ (gate 10⁻⁸). Fixed-Q A1/A2 mass difference was ≤1.368×10⁻¹² (gate 10⁻⁸). There were 210 negative calibration entries retained as signed measurements, 144 cohorts with an eigenvalue floor, and 3,500 wrong-support latent contradictions in cell 10. The pilot's failure policy kept full denominators; no incomplete cohort was silently omitted.

The first authorized `--run` stopped on its first public R truth draw: a list-backed `mineral_fraction` data frame reached `is.finite`. It produced header-only CSVs and no comparison. Its source freeze, files, and failure log remain at the root `finite_calibration/` path. Before comparative outcomes existed, the simulator recorded [the correction amendment](finite_calibration/manifests/correction_amendment.json), converted the selected public columns to a double matrix, and added a guarded `--attempt correction-001` path. The seed schedule, grid, generation laws, estimator, and public plan/config stayed fixed. `correction-001` has a separate source/environment freeze and COMPLETE execution manifest. This was a corrected completion of one pilot, not a second independent study or a clean one-pass run. The tester's older runtime notation note was corrected from `(a_i-r_i)/C_ij` to `r_i/C_ij` with a backup; accepted proof and historical tests were unchanged.

The [execution manifest](finite_calibration/attempts/correction-001/manifests/execution.json), [freeze](finite_calibration/attempts/correction-001/manifests/freeze.json), [diagnostics](finite_calibration/attempts/correction-001/results/diagnostics.json), [predictions](finite_calibration/attempts/correction-001/results/predictions.csv), [prior fits](finite_calibration/attempts/correction-001/results/priors.jsonl), and [data directory](finite_calibration/attempts/correction-001/data/) give the complete machine record. Figures: [mass RMSE](finite_calibration/attempts/correction-001/figures/mass_rmse.png), [calibration size](finite_calibration/attempts/correction-001/figures/calibration_size.png), and [stress](finite_calibration/attempts/correction-001/figures/stress.png). These paths resolve in accepted commit `07b3488dfbc6513ad9172655bebe2566fbeaa205`. At the time of the earlier internal review, the generated outputs were still uncommitted in the main checkout; the shipper subsequently committed the exact reviewed outputs.

Reproduce from a fresh checkout only after installing the pinned simCoDaR source and required R/Python environment recorded in the manifests. Retain the initial failure tree and corrected source commit/history; do not overwrite root outputs. From the repository root, run `python research/emc-part2b/prior_guided/calibration_experiment.py --config research/emc-part2b/prior_guided/finite_calibration/preregistration.json --output research/emc-part2b/prior_guided/finite_calibration --attempt correction-001 --check`, inspect its freeze, then run the same command with `--run`. The guarded correction path requires the retained initial freeze/failure artifacts and original source hash reference; a fresh checkout must include those prerequisites. It refuses a second run over existing outputs. The original root attempt's separate initial `--check`/`--run` procedure and failure are documented in the manifests and simulation correction record (`runs/REQ-20261027-finite-calibration/simulation-correction.md` in the workflow record).

No oracle prior, predictive intervals, interval coverage, confirmatory claim, field validation, or four-phase/multidomain campaign was performed. These ten frozen cells do not establish a trend beyond themselves or general superiority. A phase-sensitive measurement would be required to identify the Ti split from target measurements; any such experiment needs separate authorization and a new frozen plan. The final internal independent agent review passed on 27 October 2026; this is not external peer review or field validation.

The initial independent review STOP is retained in the workflow record. Its R1 missing deterministic/interface checks were addressed by a focused 55/55 follow-up; R2 corrected two audit diagnostic column labels without changing values. The final internal independent agent review on 27 October 2026 cleared both findings. That review inspected commit `614c770a3f27b57d702c152c361f24d34d0138a2`; accepted commit `07b3488dfbc6513ad9172655bebe2566fbeaa205` subsequently included the reviewed pilot outputs. This status update changes no numerical evidence and is not external peer review or authorization for the later study.

## Five stop answers

1. **Does independent calibration help relative to A0?** In all six matched-domain primary cells, A2 had lower held-out mass RMSE than A0, with paired mean differences from −0.05260 to −0.05981. The synthetic source center is skewed while A0 splits tied Ti mass equally, so this comparison favors a correct learned prior by design. It does not establish field transfer.
2. **When does prior information hurt?** The mismatched calibration domain in cell 8 raises A2 mass RMSE by 0.082348 versus A0 (descriptive cohort SE 0.000994). Cell 10's false support excludes true Anatase for all arms and leaves an A2−A0 gain of only −0.000165. Weak strength (7) and high QXRD noise (9) did not reverse the mass advantage in their single frozen cells. None of these cells proves a trend outside the grid.
3. **Does target assay add to the learned prior?** Variable-Quartz primary cells 4–6 show A2−A1 mass gains of about −0.052 to −0.053 and observed/latent assay RMSE gains of roughly −0.096 to −0.099. Fixed-Quartz cells 1–3 have A1=A2 within numerical tolerance; assay supplies no Anatase/Rutile contrast information.
4. **How much finite-calibration uncertainty appears?** Across-cohort fitted Anatase-center SD falls from 0.01123 to 0.00374 across fixed-Quartz n=10 to 100 and from 0.01272 to 0.00362 across variable-Quartz n=10 to 100. The [sensitivity report](PRIOR_SENSITIVITY.md) separates center variability, estimated specimen covariance T, and approximate mean covariance M. Covariance-estimation uncertainty is unintegrated, so these are not predictive-coverage results.
5. **What next?** The pilot justifies proposing a separately authorized phase-sensitive measurement test for the exact Ti null direction, with mismatch and support robustness considered in a new frozen design. The pilot stopped here after final internal independent agent review PASS; it did not itself authorize an automatic larger experiment. A later, separately authorized bounded phase-sensitive study has its own frozen design and independent pre-execution review gate.
