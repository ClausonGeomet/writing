# R1 approved contract amendment: frozen legacy baseline

Date: 2026-10-26. Approved under the explicit focused planner dispatch selecting the legacy-acceptance option in `review-first.md`, R1. This is a contract amendment for downstream verification and disclosure, not an executed test result or final reviewer clearance. Preserve the original specifications as historical requirements. This amendment supersedes only the baseline residual convention and clarifies the reported objective semantics; all other requirements and tolerances remain unchanged.

## Original requirement, observed behavior, approved convention

| Surface | Original requirement / prior description | Observed behavior reported by reviewer | Approved legacy convention |
| --- | --- | --- | --- |
| Part A baseline residual, all methods | Original test contract requires `a - reconstructed`, where reconstruction is `C@m`. | Frozen baseline returns `reconstructed - a`. | Formally accept and disclose `reconstructed - a` for frozen Part A. |
| Part B allocation/prediction residual | `a - reconstructed`. | Uses the requested sign. | Retain `a - reconstructed`, opposite to Part A by design under this amendment. |
| Part A `objective` field | Tester contract specifies squared weighted residual without a factor of one half; ordinary-fit labels and provenance descriptions can imply the minimized unweighted objective for `pinv`/`nnls`. | Field always uses sigma weighting, including at ordinary `pinv`/`nnls` solutions. | Disclose a sigma-weighted squared-residual diagnostic evaluated at the returned solution, without a factor of one half. For ordinary `pinv`/`nnls`, this is not the minimized unweighted objective. Sigma defaults to ones; fitting semantics remain unchanged. |

**The old tester oracle was wrong against the original contract.** It asserted the legacy sign, and its historical PASS does not establish compliance with the original `a - reconstructed` requirement. Do not retrospectively describe it as an original-contract PASS. A new focused PASS can establish compliance with this explicitly amended contract only.

Residual sign reversal does not affect mineral estimates, reconstructions, squared errors, squared objectives, or reported benefit conclusions. Clarifying the objective field does not change any fitting objective, estimates, or frozen numeric values.

## Focused nonzero test witness

Use the baseline black-box interface with `C=[[1],[1]]`, `a=[0,1]`, `sigma=[1,2]`, and `method="wnnls"`.

- Expected `m=[0.2]`, with the existing estimate tolerance `1e-6`.
- Expected reconstruction `[0.2,0.2]`, with existing replay tolerance `1e-10`.
- Original-contract residual: `[-0.2,0.8]`.
- Approved frozen legacy Part A residual: `[0.2,-0.8]`. Check the nonzero returned residual against `reconstructed - a` with the unchanged `1e-10` replay tolerance.
- Part B retains `a - reconstructed`. At the same reconstruction its residual would be `[-0.2,0.8]`. This compares conventions; it does not assert that Part B allocation returns the weighted-baseline solution on this input.

These are expected values, not executed evidence. Tester must run the focused check and record its actual outcome. No new formulas, numerical thresholds, or broad validation campaign are introduced.

## Preservation and downstream handoff

Preserve original `spec.md`, `test-spec.md`, and `sim-spec.md`, frozen v3 artifacts and hashes, historical records, and original timings. Do not regenerate fits or replace frozen records. No target source/header edits or corrected derivative interface are necessary for this chosen resolution: explicit legacy-semantics disclosure is sufficient.

Tester consumes this amendment with its own permitted contract, corrects the audit expectation, and provides a named Before/After/deviation table separating the original requirement, historical mismatched oracle, and approved legacy convention. Keep the new focused result separate from historical passes.

Scriber discloses both residual conventions and the meaning of the baseline objective field. Ordinary `pinv`/`nnls` fitting remains ordinary fitting even when the reported diagnostic uses sigma weighting. Explain historical labels without silently rewriting frozen values.

This role-neutral amendment contains shared observable semantics and the explicitly authorized witness only, with no implementation algorithm or other pipeline specification. Tester and scriber need not read another pipeline's specification. R2 and final reviewer clearance remain separate; this amendment alone does not clear the review STOP.
