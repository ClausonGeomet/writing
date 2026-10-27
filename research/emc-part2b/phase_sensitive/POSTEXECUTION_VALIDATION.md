# Post-execution validation lineage

**Status (27 October 2026): numerical results complete; final internal scientific review STOP pending independent saved-result validation and re-review.** This page records a correction to the *validation harness* after the bounded synthetic experiment. It does not amend the preregistered design, generate new observations, change the reported estimates, or claim a final PASS. The original [STOP review](REVIEW_HISTORY/FINAL_REVIEW_STOP.md) remains the decision in force.

## Why an additive check was needed

The focused draw-free follow-up to the final review found three defects in the original `measurement_experiment.py` harness:

1. The resource gate could admit a sampling-capable subprocess after its time or logical-output cap was exhausted.
2. A failed matched block suppressed common error metrics but left own-channel error and envelope-derived fields available.
3. The historical `--verify` could accept a feasible but nonoptimal saved mass when a copied output hash and reported successful diagnostics were forged. File integrity alone did not prove the convex fit's optimality.

The first additive correction, commit `bf1b508`, introduced `measurement_experiment_corrected.py` and `test_experiment_corrected.py`. It provides a draw-free saved-input verifier, reconstructs permitted observations and priors, checks every saved mathematical certificate independently of reported KKT values, completes failed-block masking, and exercises the resource gate through a required harmless callback. An independent tester audit found one further defect: when a required prediction row was *absent*, the corrected scorer counted 29 successes of 30 planned but zero failures. Commit `807a0b8` changed that count to one failure and added a regression case. The tester's full certificate for this second source identity and final independent re-review are pending; an earlier valid-input PASS for `bf1b508` does not clear its missing-row failure or certify `807a0b8`.

## Separate identities and preservation

The correction-002 execution remains bound to its original five computational sources: `measurement_experiment.py`, `generation.R`, `test_experiment.py`, `phase_models.py`, and `test_phase_models.py`. The six frozen design files, original `run/` tree and all 116 archived regular files remain byte-exact. The source-file index SHA256 is `2876fedc6ebc72a3c9f83874ab5f717f51c4259eea6d1824e03b93a9bd9973dd`; the reviewed design digest is `7f011812d4ecfee33f8184bab59ea0591d98580215782d9dc509e547bfacabad`. The new verifier must make its own later source, environment, input and evidence record outside `run/`. It cannot become the source that produced the historical data or backdate a pre-draw certificate.

The historical `measurement_experiment.py --verify` remains an integrity check under its exact original source and environment. Its recorded PASS is weaker than the additive mathematical certificate: it did not independently reconstruct and certify every saved fit. The independent tester's audit of `bf1b508` reported a valid-input certificate covering 350,250 estimates, 7,750 priors, 42,750 bounds, 26,250 information records, 600 acquisition rows, 21,250 cohort metrics and 3,050 paired records, plus an independent all-row mass/likelihood oracle. Its missing-row adversarial check BLOCKED that source identity. The `807a0b8` source needs a fresh independent certificate and negative probe under its own hash before this lineage can support a new review decision.

The required detailed prediction/prior/update key inventory was absent before the first draw; only aggregate planned counts existed then. Later exhaustive reconstruction of keys is **retrospective**. It can check saved-row completeness but does not cure or erase this prospective process deviation. No pre-draw timestamp or prospective PASS is assigned to it.

## Draw-free saved-result validation

First restore and byte-check the six archives using [the archive recipe](artifacts/README.md). From the **original main-checkout root** at `/home/matt/Documents/GitHub/writing`, with the recorded Python/R environment, original absolute paths, accepted predecessor files, pinned simCoDaR installation, historical seed specifications, and runtime `/tmp/emc-part2-workspace/writing/runs/REQ-20261027-phase-sensitive/simulator-correction-spec.md` available at their recorded paths, run the additive verifier with a fresh output directory outside `run/`:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 \
/tmp/bif-v3-world1-validation-env-20261023/bin/python \
  research/emc-part2b/phase_sensitive/measurement_experiment_corrected.py \
  --input-attempt /home/matt/Documents/GitHub/writing/research/emc-part2b/phase_sensitive/run/attempts/correction-002 \
  --evidence-output /tmp/emc-phase-postexecution-validation-fresh \
  --verify-saved
```

Choose an unused `--evidence-output` path for each attempt. The command reads saved observations/results and writes separate post-execution evidence; it has no generation option and authorizes no new draws. It should fail if original-path, source/input identity, environment or other required provenance is unavailable. The correction specifications define the certificate and unchanged tolerances; the checker output and independent tester audit, once complete, determine its actual result. Do not infer a PASS from this command's presence in documentation.

Later delivery may package the correction specification, audits, scripts and validation results in a separate validation-review archive. That package is not yet part of this study's committed `artifacts/` set, so this page gives no archive link or completed verdict for it. The bounded numerical estimates in [RESULTS.md](RESULTS.md) remain the reported scientific results; their final clearance awaits independent review of this additive evidence and the disclosed process deviation.
