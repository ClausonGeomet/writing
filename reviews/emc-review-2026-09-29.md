# EMC article review

Date: 29 September 2026  
Article: *From elements to minerals: what an assay can tell us*  
Repository: `ClausonGeomet/writing`  
Original article blob: `6ecace6641376a76b2ae349db696c849b461c96b`  
Revised article blob: `e121b697677b0dfdc416793df8c1ab8477943490`

## Overall assessment

The strongest part of the article is its practical question: which mineral distinctions does an assay support, and what measurement would change the answer? The quartz–anatase–rutile example is a good way to explain this without hiding behind a large dataset or an algorithm.

The original core linear algebra and conditional grid calculation were sound. The main problems were attribution, a few over-broad interpretations, an optional prediction-code bug and too much repeated explanation. This revision is a teaching article informed by the literature, not a reproduction of each paper's method or an independently validated operational EMC model.

## Editorial changes

The opening now gets to the assay/mineralogy gap once rather than twice. The literature discussion is organised around contributions to the argument, rather than mainly by publication date. The recurring “Question / What this means / What it does not mean” pattern has been removed; it made otherwise useful qualifications feel repetitive.

The worked mixture, plain explanations and phrases such as “This is bookkeeping, not inference” remain. First-person interpretation is used sparingly, without adding claims about professional experiences that were not in the source. The conclusion returns to the mining decision rather than ending with another catalogue of algorithms.

The long genetic-search demonstration, decorative calibration/network graphic, prior-density gallery and duplicate sampler examples have been replaced with shorter explanations. The two main plots now each have a distinct job: showing the null direction and showing an update from genuinely different information. A sensitivity table adds more evidence than another general warning about priors.

This is a substantial proposed edit, not merely a spelling pass. The original is retained in Git history and the changes are on a review branch, not merged into the live site.

## Technical findings and changes

### 1. Optional NumPyro posterior prediction: corrected

The original model defaulted `y` to the observed value, and its `Predictive` call did not override that default. That keeps the response observed instead of drawing a new response. The model now defaults to `y=None`; fitting explicitly supplies `y=0.28`, while prediction explicitly passes `y=None`.

The standalone example also imports ArviZ itself, rather than relying on an import in a different, non-executed PyMC cell. These changes follow the [NumPyro Predictive documentation](https://num.pyro.ai/en/stable/utilities.html).

This bug was in an optional `eval: false` snippet. It does not invalidate the original grid-based posterior-predictive demonstration. The corrected sampler has been syntax-checked and its observation/prediction call contract inspected programmatically, but it has not been run in this review environment.

### 2. Conditional prior versus original prior: clarified and derived

For `(quartz, anatase, rutile) ~ Dirichlet(6,2,2)`, define `S = anatase + rutile` and `u = anatase / S`. Then `S ~ Beta(4,6)` and `u ~ Beta(2,2)` are independent. Conditioning on the exact assay fixes `S=0.40`, giving `anatase = 0.40u`.

The density on that fixed-total line is a scaled Beta(2,2), not the original marginal Beta(2,8) distribution of anatase. The article now includes the Jacobian factorisation and scale factor. “Posterior equals prior” applies to the unresolved share in this conditional example, not to the complete mineral vector before using the assay.

The numerical grid is checked against independent quadrature and a finer grid. This improves the explanation without suggesting that the original grid's basic calculation was wrong.

### 3. Closure and uniqueness: qualified

In this toy, closure is redundant: its row is the sum of the two assay rows. Elsewhere it can increase the rank. Non-negativity can also leave one feasible boundary solution even when the equality matrix is rank-deficient.

The revision therefore uses the augmented system `[C; 1^T]` and discusses feasible null directions, rather than treating the number of assays alone as a complete test. The checker includes examples where closure adds information and where a rank-deficient system has a unique non-negative boundary solution.

### 4. Allocation order: the original comparison changed two things

The original two recipes changed both analyte order and titanium tie-break. In this particular matrix, reversing order while keeping the recipient fixed does not change the final mineral composition.

The revised implementation crosses both orders with both titanium recipients. The four results make it clear that the phase-selection rule, not order alone, drives the difference here. This does not imply that order is irrelevant in more complicated sequential allocations.

### 5. Calibration claims: made consistent with the implementation

The original 70/30 weights were hard-coded; there were no calibration labels or weight-fitting procedure in the code. They are now explicitly described as illustrative assigned weights. No fitted labels, hold-out evaluation or literature-method reproduction is claimed.

A genuine calibration extension would need paired observations, a defined objective, a training/tuning/evaluation split and mineral-error results. A mass-balance check alone would not establish that the mineral predictions are right.

### 6. Optimisation and uncertainty: separated

The pseudoinverse is the minimum-norm least-squares solution, not “usually” such a solution. NNLS adds non-negativity but does not enforce closure; the toy happens to close. LP ranges are coupled feasibility bounds, not probability intervals.

The residual chart magnified by 10^15 has been removed. Floating-point differences between effectively exact fits do not distinguish mineralogical quality. A compact table is enough.

The fixed linear, simplex-constrained least-squares example is convex. A genetic algorithm is unnecessary for this example, and repeated seeds are not posterior samples. The explanation retains a role for nonlinear or discrete search without implying that search identifies an unobserved phase distinction.

### 7. Assay basis and uncertain mineral chemistry: made more explicit

Mineral closure is distinguished from normalising an incomplete or oxide-equivalent assay suite. A normalised-response formula is included, along with the need to handle potentially singular covariance on closed assays.

The article keeps matching chemistry and abundance posterior draws together. It now shows the covariance term in `E[Cm]`, explaining why a product of separate posterior means generally differs from the mean reconstructed assay.

### 8. Reproducibility: improved, with limits recorded

Main code is folded. `execute.error` is now false, so an execution error stops a Quarto render instead of being printed into the article. A reusable checker compiles all cells, executes the seven main cells in order, draws the two plots, and explicitly skips the optional sampler.

No repository dependency or lockfile changes were needed. Local validation used the installed environment, not a fresh `uv sync --frozen` environment. That distinction should remain visible until the repository environment and complete site render have been checked.

## Literature review: evidence and limits

| Source | Use in the revision | Verification boundary |
| --- | --- | --- |
| Johnson, Chu and Hussey (1985) | Component-property uncertainty and simultaneous-equation workflow | Publisher abstract; not a reproduction of its soil calculations |
| Whiten (2008) | SVD, non-uniqueness and pseudoinverse interpretation | Paper abstract/bibliographic record, supported by direct algebra and NumPy documentation; no full example replication |
| Berry, Hunt and McKnight (2011) | QXRD-supported calibration and error-weighted combination | AusIMM proceedings abstract; no claim to reproduce exact LP preferences |
| Lund, Lamberg and Lindberg (2013) | Operational bridge between chemistry, diagnostic measurements and actual mineralogy | Publisher abstract/preview; original Malmberget data not reanalysed |
| Parian and colleagues (2015) | Combined EMC and quantitative XRD | Publisher abstract/preview; removed the earlier unverified “constrained non-negative factorisation” attribution |
| Mena Silva and colleagues (2018) | LS-XRD versus regression-based R-XRD in a particular deposit | Publisher article text/search-accessible sections; no independent reanalysis or universal performance claim |
| Rodrigues and colleagues (2023), Pedras | Bayesian direction | Conference abstract only; not enough to implement or validate its exact model |

Public references replace the original “local background paper” placeholders. The post now distinguishes model-based modal estimation from a prescribed normative allocation. It does not claim that the generic examples reproduce HSC software.

Key primary records: [Johnson](https://doi.org/10.1346/CCMN.1985.0330204), [Whiten](https://doi.org/10.1080/08827500701257860), [Berry](https://www.ausimm.com/publications/conference-proceedings/first-ausimm-international-geometallurgy-conference-geomet-2011/estimating-mineralogy-in-bulk-samples/), [Lund](https://doi.org/10.1016/j.mineng.2013.04.005), [Parian](https://doi.org/10.1016/j.mineng.2015.04.023), [Mena Silva](https://doi.org/10.3390/min8080325), [Pedras](https://doi.org/10.48380/2r85-w395).

## Numerical checks completed

All seven main cells passed, all eight Python cells compiled, and both Plotnine figures were drawn and visually inspected. Independent checks also covered the additional phase-sensitive row, informative closure, non-negative boundary uniqueness and conditional-prior distinction.

| Calculation | Result |
| --- | --- |
| Pseudoinverse mineral composition | 60/20/20 wt% |
| NNLS result in the tested environment | 60/40/0 wt%; its particular tie-break is not mineralogical evidence |
| Closed non-negative fit | 60/20/20 wt% |
| LP bounds | Quartz 60–60; anatase 0–40; rutile 0–40 wt%, with coupled Ti total 40 |
| Base-case posterior anatase mean | 27.20 wt% |
| Base-case central 90% credible interval | 21.01–33.25 wt% |
| Beta(8,8) share prior, same observation/error | Posterior mean 24.85 wt% |
| Beta(6,2) share prior, same observation/error | Posterior mean 29.10 wt% |
| Simulated central 90% response-predictive interval | 18.01–36.11 wt%, fixed seed; not a mineral credible interval |

The posterior mean agrees with adaptive quadrature within 10^-7 on the fraction scale. The 4,000-bin and 8,000-bin grid mean/interval summaries agree within 2×10^-6 on that scale. These are numerical checks of a stipulated toy model, not evidence of real-data EMC accuracy.

Environment: Python 3.13.5; NumPy 2.3.5; SciPy 1.17.0; pandas 2.2.3; Plotnine 0.15.3; IPython 9.14.0.

## Before publication

Run the checker in the locked repository environment, then render and inspect the complete Quarto page, particularly equations, callouts, narrow-screen tables and figure legends. Quarto was unavailable here. The optional NumPyro cell also needs an actual runtime check and comparison against the grid before its validation disclaimer is removed.

The existing `rendered-report.html` is intentionally unchanged: it is a historical pre-review snapshot, not a render of the revised source. The accompanying README now makes that distinction explicit. No claim is made that this review rendered or deployed the live page.

My editorial recommendation is to keep this as a focused foundational tutorial. A later deposit-specific paper replication or measured-data calibration would be a stronger next article than restoring a large algorithm catalogue to this one.
