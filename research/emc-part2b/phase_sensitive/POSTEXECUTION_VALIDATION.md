# Post-execution validation lineage

Final internal independent scientific convergence re-review: PASS WITH NOTE for bounded exploratory research artifact delivery in a stacked draft PR only. R1 validation and R2 report reproduction are cleared. The missing prospective detailed-key inventory remains a disclosed historical process deviation; exhaustive reconstruction and certification are retrospective. The original STOP, pre-execution review, failed attempts and unchanged numerical results remain preserved. This is not confirmatory evidence, external peer review, field validation, or authorization for merge, article publication, new draws or expanded work. See the exact [final scientific re-review](FINAL_SCIENTIFIC_REVIEW.md), [historical STOP](REVIEW_HISTORY/FINAL_REVIEW_STOP.md), and [pre-execution review](PREEXECUTION_REVIEW.md).

This page records a correction to the *validation harness* after the bounded synthetic experiment. It does not amend the preregistered design, generate new observations, or change the reported estimates. The original [STOP review](REVIEW_HISTORY/FINAL_REVIEW_STOP.md) remains preserved as history.

## Why an additive check was needed

The focused draw-free follow-up to the final review found three defects in the original `measurement_experiment.py` harness:

1. The resource gate could admit a sampling-capable subprocess after its time or logical-output cap was exhausted.
2. A failed matched block suppressed common error metrics but left own-channel error and envelope-derived fields available.
3. The historical `--verify` could accept a feasible but nonoptimal saved mass when a copied output hash and reported successful diagnostics were forged. File integrity alone did not prove the convex fit's optimality.

The first additive correction, commit `bf1b508`, introduced `measurement_experiment_corrected.py` and `test_experiment_corrected.py`. It provides a draw-free saved-input verifier, reconstructs permitted observations and priors, checks every saved mathematical certificate independently of reported KKT values, completes failed-block masking, and exercises the resource gate through a required harmless callback. An independent tester audit found one further defect: when a required prediction row was *absent*, the corrected scorer counted 29 successes of 30 planned but zero failures. Commit `807a0b8` changed that count to one failure and added a regression case. The independent tester passed all 33 deterministic tests and issued a fresh full draw-free certificate for commit `807a0b8dcdc2913e07033385e9d996ec99503352`. The earlier valid-input PASS for `bf1b508` remains historical and did not clear its missing-row failure. The [final scientific re-review](FINAL_SCIENTIFIC_REVIEW.md) subsequently gave PASS WITH NOTE for bounded exploratory draft delivery only.

## Separate identities and preservation

The correction-002 execution remains bound to its original five computational sources: `measurement_experiment.py`, `generation.R`, `test_experiment.py`, `phase_models.py`, and `test_phase_models.py`. The six frozen design files, original `run/` tree and all 116 archived regular files remain byte-exact. The source-file index SHA256 is `2876fedc6ebc72a3c9f83874ab5f717f51c4259eea6d1824e03b93a9bd9973dd`; the reviewed design digest is `7f011812d4ecfee33f8184bab59ea0591d98580215782d9dc509e547bfacabad`. The new verifier wrote its later source, environment, input and evidence record outside `run/`. It did not become the source that produced the historical data or backdate a pre-draw certificate. Final source-record SHA256: `e382c97d4fcfe5a93141f574f2a75521cd22d8db655b3485f3deab1724d01da2`; validation-record SHA256: `79a3814c5c5c447e46c73b8c53536f96bdce71577ba33fdb70e871ba3cf0b316`.

The historical `measurement_experiment.py --verify` remains an integrity check under its exact original source and environment. Its recorded PASS is weaker than the additive mathematical certificate: it did not independently reconstruct and certify every saved fit. The independent tester's audit of `bf1b508` reported a valid-input certificate covering 350,250 estimates, 7,750 priors, 42,750 bounds, 26,250 information records, 600 acquisition rows, 21,250 cohort metrics and 3,050 paired records, plus an independent all-row mass/likelihood oracle. Its missing-row adversarial check BLOCKED that source identity. The `807a0b8` source received its own independent certificate and missing-row negative probe. The probe produced 30 planned, 29 successful and one failed prediction; all 13 matched arm rows and 90 affected contrasts were masked, while unaffected rows stayed identical. The full certificate passed 350,250 saved estimates, 7,750 priors, 42,750 bounds, 26,250 information records, 600 acquisition rows, 21,250 cohort metrics, 3,050 paired records and one diagnostic record with zero draws and zero generator calls. The independent tester found all 116 original run files and six archives byte-identical. See runtime `audit-final-validation.md` and `final-validation-evidence/` for exact checks and hashes. This numerical PASS was considered separately in the subsequent [scientific re-review](FINAL_SCIENTIFIC_REVIEW.md).

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

Choose an unused `--evidence-output` path for each attempt. The command reads saved observations/results and writes separate post-execution evidence; it has no generation option and authorizes no new draws. It should fail if original-path, source/input identity, environment or other required provenance is unavailable. The correction specifications define the certificate and unchanged tolerances. The independent tester ran this command against the preserved correction-002 attempt and recorded PASS in runtime `audit-final-validation.md`; another invocation must be judged by its own output.

The existing [supplemental validation-review evidence package](artifacts/VALIDATION_REVIEW_README.md) documents the archive and its exclusions; it predates the separately published final scientific re-review and is not amended to include it. The bounded numerical estimates in [RESULTS.md](RESULTS.md) remain the reported scientific results, cleared only for exploratory draft delivery with the disclosed process deviation. Ruff and mypy were exploratory, nonconfigured static checks: both exited 1 with pre-existing findings (six Ruff findings matching parent `bf1b508` by rule/message, five mypy errors in unchanged `phase_models.py`). Exact outputs are retained in runtime `final-validation-evidence/`; these are separate from the 33 passing deterministic tests.
