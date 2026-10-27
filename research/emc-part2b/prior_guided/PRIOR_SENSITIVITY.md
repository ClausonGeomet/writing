# Prior sensitivity in the bounded QAR pilot

**27 October 2026 · exploratory, descriptive; independent scientific review pending.** This companion to the [results](FINITE_CALIBRATION_RESULTS.md) uses the ten frozen cells only. Each cell has 50 independent calibration/evaluation cohorts, with 40 held-out evaluation specimens per cohort. The full precision [diagnostics](finite_calibration/attempts/correction-001/results/diagnostics.json), [cohort metrics](finite_calibration/attempts/correction-001/results/cohort_metrics.csv), and [paired contrasts](finite_calibration/attempts/correction-001/results/paired_summary.csv) are authoritative.

## Frozen factors and cohort behavior

A2−A0 and A2−A1 below are mean paired mass RMSE differences ± descriptive SE across 50 cohort differences. Negative favors A2. `SD center A` is the across-cohort SD of the fitted Anatase center. `T_A` is the mean Anatase diagonal of the estimated physical within-population covariance; `M_A` is the mean Anatase diagonal of the approximate calibration-mean covariance. These are different quantities, and neither incorporates covariance-parameter uncertainty.

| Cell | Factor | n | QXRD SD | Strength | A2−A0 mass | A2−A1 mass | SD center A | Mean T_A | Mean M_A | Floor / 50 |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | fixed | 10 | 0.02 | 1 | -0.05260 ± 0.00081 | 2.15e-15 ± 2.8e-15 | 0.01123 | 0.00099 | 0.00013 | 39 |
| 2 | fixed | 30 | 0.02 | 1 | -0.05367 ± 0.00076 | -1.02e-16 ± 2.57e-15 | 0.00672 | 0.00133 | 5.34e-05 | 44 |
| 3 | fixed | 100 | 0.02 | 1 | -0.05361 ± 0.00080 | 9.87e-16 ± 3.35e-15 | 0.00374 | 0.00126 | 1.53e-05 | 48 |
| 4 | variable | 10 | 0.02 | 1 | -0.05664 ± 0.00101 | -0.05222 ± 0.00145 | 0.01272 | 0.00319 | 0.00035 | 1 |
| 5 | variable | 30 | 0.02 | 1 | -0.05981 ± 0.00107 | -0.05208 ± 0.00115 | 0.00790 | 0.00303 | 0.00011 | 0 |
| 6 | variable | 100 | 0.02 | 1 | -0.05927 ± 0.00121 | -0.05342 ± 0.00105 | 0.00362 | 0.00303 | 3.29e-05 | 0 |
| 7 | weak | 30 | 0.02 | 0.01 | -0.06068 ± 0.00103 | -0.05342 ± 0.00106 | 0.00855 | 0.00314 | 0.00011 | 0 |
| 8 | mismatch | 30 | 0.02 | 1 | 0.08235 ± 0.00099 | -0.01057 ± 0.00076 | 0.01663 | 0.00797 | 0.00027 | 0 |
| 9 | noisy_qxrd | 30 | 0.08 | 1 | -0.05783 ± 0.00103 | -0.05284 ± 0.00130 | 0.01344 | 0.00317 | 0.00025 | 12 |
| 10 | wrong_occurrence | 30 | 0.02 | 1 | -0.00016 ± 9.61e-06 | -0.02103 ± 0.00099 | 5.46e-17 | 0.00306 | 0.00011 | 0 |

## Interpretation by frozen comparison

- **Calibration size.** In fixed-Quartz primary cells 1–3, center SD for Anatase is 0.01123, 0.00672, 0.00374 at n=10, 30, 100. In variable-Quartz primary cells 4–6 it is 0.01272, 0.00790, 0.00362. These independent cells suggest less center variation with larger n in this synthetic setup; they do not prove a general monotone trend. The estimated T reflects individual specimen heterogeneity and remains distinct from M, the approximate center-mean uncertainty.
- **Fixed Quartz and exact null space.** Cells 1–3 have A1=A2 within 1.37×10⁻¹², because their feasible assay reconstruction is constant. The bulk Ti assay has no Anatase/Rutile split signal in any cell. Weakening the quadratic penalty in cell 7 changes how much variable Quartz is fitted; it cannot identify an exact Ti-contrast null direction. A positive penalty chooses a solution using calibration information.
- **Weak strength.** Cell 7 uses 0.01 against cell 5’s 1. Its A2−A0 mass contrast is -0.06068 ± 0.00103 versus cell 5’s -0.05981 ± 0.00107. This single stress comparison does not show weak-prior harm and is not a strength-response curve.
- **Wrong mean/domain.** Cell 8 changes the calibration generator center to (.6,.3,.1) while evaluation remains (.6,.1,.3). A2−A0 mass is 0.08235 ± 0.00099; the mean difference is +0.082348, an adverse result. The mean fitted Anatase center minus held-out cohort truth mean is 0.20015. This cohort-mean comparison is retrospective, noisy, and not a known population bias parameter. It was not used to fit the prior.
- **QXRD noise.** Cell 9 uses known lab SD .08 against cell 5’s .02; 189 signed calibration entries are negative, while 12/50 cohorts use an eigenvalue floor. Its A2−A0 mass contrast is -0.05783 ± 0.00103. The higher noise makes finite calibration and deconvolution less stable, but this one cell does not establish a trend beyond the frozen settings.
- **Wrong occurrence.** Cell 10 imposes Anatase=0 on A0, A1, and A2, despite 3,500 calibration/evaluation latent contradictions. A2−A0 mass is only -0.00016 ± 9.61e-06. The failure is a shared false support assumption; it should not be attributed to the calibrated prior alone.

## What the variance calculation does and does not say

The fitted center is the Euclidean projection of the raw noisy calibration mean. In raw-mass tangent coordinates, the estimator deconvolves known QXRD measurement covariance, floors eigenvalues at 0.0001, shrinks by 0.25 toward isotropy, and calls the resulting within-population covariance T. It sets M=(T+R_z)/n for approximate uncertainty of the *unprojected* mean, then uses T+M to form the plug-in precision. Projection, floor and shrinkage can alter center and variance estimates. The 144 floor-affected cohorts, especially 39/44/48 in fixed-Quartz n=10/30/100, show these choices were active, not merely formal. The calibration covariance itself is estimated from finite data and its uncertainty is not integrated. No predictive intervals, coverage, posterior probability, or field error rates follow from T or M.

The QAR assay has response matrix `[[1,0,0],[0,1,1]]`: the vector (0,1,−1) lies in its null space even under closure. No stronger penalty or more precise bulk Ti assay can create information about this direction. There is also irreducible specimen variation around a learned population center. The correct-model synthetic generator deliberately favors the skewed calibrated center over A0’s equal-Ti minimum-norm split; this is a design-specific result, without general empirical or field validity.

## Decision boundary

The five stop answers are in the [results report](FINITE_CALIBRATION_RESULTS.md): calibrated A2 helps relative to A0 in the six matched primary cells; domain mismatch hurts; the assay helps A2 beyond A1 only for total Ti/Quartz in variable-Quartz cells; finite center variation decreases across the declared sizes while covariance uncertainty remains unintegrated; and a **separately authorized phase-sensitive measurement study** is the sensible next scientific test. A larger experiment would require a new freeze and independent review. This bounded pilot has ended; no automatic expansion, GA, four-phase, multidomain or confirmatory campaign follows.
