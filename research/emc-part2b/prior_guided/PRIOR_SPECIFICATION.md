# Independent geological priors: proposed specification

**Checkpoint NOT ACCEPTED; calibration and campaigns NOT RUN; external validity unestablished.** These are candidate models, not fitted geological distributions. The prior must be estimated without evaluation mineral truth. A synthetic generator's known parameters are an oracle diagnostic, not evidence of independently learned priors.

## Evidence ledger to create before fitting

| Evidence | Intended role | Required record and boundary |
| --- | --- | --- |
| Domain and petrography | Phase occurrence/support and domain assignment | Source, site/drillhole, sample IDs, date, classifier version, ambiguity; freeze classifier on training data |
| QXRD/QEMSCAN | Noisy phase abundance calibration or target auxiliary likelihood | Method, detection limit/censoring, laboratory bias/precision; never count one measurement in both prior training and target likelihood |
| EPMA/mineral standards | C and chemistry uncertainty | Phase, method, revision, within/between-domain variability; propagate identically to A2/B |
| Bulk assay | Target likelihood | Units, detection limit, covariance, missingness and signed errors; evaluation mineral labels excluded |
| Split manifest | Leakage control | Group correlated site/drillhole samples; disjoint calibration, tuning and test groups |

An unsupplied geological calibration dataset cannot be replaced by held-out mineral masses, their empirical covariance, or a test-derived domain label. A target's auxiliary measurement enters **once** as a likelihood, or is withheld; cross-fitting requires grouped splits.

For domain d, distinguish population variation Σ_pop,d=Var(M_new|θ_d), uncertainty in learned mean Var(μ_d|D_cal), uncertainty in covariance/hyperparameters, and laboratory error in observed Y_cal. Model Y_cal|M_cal with its measurement error and censoring. The predictive variance obeys Var(M_new|D_cal)=E[Var(M_new|θ_d,D_cal)|D_cal]+Var(E[M_new|θ_d,D_cal]|D_cal). Raw QXRD sample variance mixes geological and analytical variation. Sparse calibration calls for stated shrinkage/low-dimensional hierarchy and sensitivity, not an unrestricted inverse covariance. Integrate over θ or use calibration-only grouped bootstrap in a future study.

## Candidate models and constraints

Reference quadratic: μ_d from independent calibration, W assay precision, K_d PSD prior precision/regularizer, λ frozen on calibration/validation only; J is defined in MATHEMATICAL_DESIGN. For closed compositions, raw covariance is singular. Let UᵀU=I,Uᵀ1=0, tangent covariance V≻0; K=UV⁻¹Uᵀ is a declared tangent regularizer. For incomplete inventory, represent total known-phase mass s≤b separately from conditional composition x=m/s. Handle s=0 explicitly; never silently set s=b or assign zero chemistry to missing phases.

On active support S, Dirichlet α_j>0 gives μ_j=α_j/α₀, Var(x_j)=μ_j(1−μ_j)/(α₀+1), Cov(x_j,x_l)=−μ_jμ_l/(α₀+1). One concentration α₀ cannot encode arbitrary covariance; positive α gives no structural-zero point mass. Its negative-log penalty is convex only when all α_j≥1. A logistic-normal alternative uses z=Uᵀlog x∼N(η,V), x=closure(exp(Uz)); it requires x_j>0 on S. V is logratio, not raw-mass covariance. A quadratic in z is a coordinate regularizer. The negative log density on simplex Lebesgue measure also contains Σ_j log x_j (up to constants); declare the base measure and use the same one in A2/B. MAP and posterior mean are different estimands.

Use a support mixture P(S|d,evidence), conditional positive composition on S, and exactly zero masses outside. Single-phase support is a point mass. Nondetection is censored evidence, not proof of absence and not a reason for a universal pseudocount. Compare wrong-support and false-positive costs. Domain uncertainty is a weighted mixture Σ_d P(d|available evidence)p(m|d,D_cal); report misclassification sensitivity. Correlated occurrence may require joint support probabilities rather than independent Bernoulli factors.

## Prospective sensitivity ledger

Vary accurate versus biased mean, weak versus strong precision, small calibration n, covariance shrinkage, QXRD error/censoring, incorrect domain, incorrect support/structural zero, missing phase, and uncertain C. Treat true abundance variability, estimated-mean uncertainty and measurement error as separate axes. For fixed C and additive ε, Cov(a)=CΣ_m Cᵀ+Σ_ε+C Cov(m,ε)+Cov(ε,m)Cᵀ; cross terms vanish only under an independence or zero-correlation assumption. With random C, E[Cm] generally differs from E[C]E[m]. The same chemistry hierarchy and auxiliary data must be available to every informed arm.

Decision gate: independent calibration and its provenance, frozen prior family/hyperparameters, grouped split, measurement model and domain/support sensitivity must be reviewable before empirical claims. None has yet been executed.
