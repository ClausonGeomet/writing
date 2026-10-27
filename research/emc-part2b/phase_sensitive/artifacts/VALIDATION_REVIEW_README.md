# Post-execution validation and review evidence

`validation-review-evidence.tar.gz` is an additive, draw-free evidence package. It restores under `validation_review/`, separate from the immutable `run/` tree. It contains 121 exact runtime files: the correction bridge and both correction specifications; every runtime audit and review version (including the initial scientific STOP and current tester PASS); source-correction handoffs; focused follow-up scripts and results; corrected-harness scripts, results and its earlier full certificate; and final-validation scripts, 33-test result, evidence index and new full certificate. The runtime report plotting script is included. The later [final independent scientific re-review](../FINAL_SCIENTIFIC_REVIEW.md) gave PASS WITH NOTE for bounded exploratory stacked draft delivery only; the [initial STOP](../REVIEW_HISTORY/FINAL_REVIEW_STOP.md) remains preserved. This archive predates that re-review and does not contain it. The tester's PASS is a validation finding distinct from the scientific verdict.

The gzip/tar output is deterministic: lexical members, gzip timestamp zero and no filename, tar timestamps and numeric ownership zero, empty owner/group names, regular files only, mode 0644. The archive is 217,141 bytes (under 90 MiB), SHA256 `90879684b7cf162bf533ba46aca785cecdc1db1dbf47c8eb58a5d88f5f094b4e`. `validation-review-files.sha256` is the governing exact-byte index for included members (SHA256 `55af5cd9797d980921bea0f72dee6553715da6be70ee48cdf34a025b735f03a6`). `validation-review-exclusions.json` names all 38 omitted files or fixture symlinks individually, with a reason and recreation route (SHA256 `265cc51a40c0fabfc61590ef16e9db74b1d07c41d71778f6e45bb64dbf1ed091`). The copied fixtures and raw datasets are available from the six original run archives plus the retained probe scripts; no fixture symlink or raw dataset is embedded here. Tool caches have no evidentiary role.

## Restore and check

From the repository root, in a checkout with the six original archives available:

```bash
set -euo pipefail
study=research/emc-part2b/phase_sensitive
sha256sum "$study/artifacts/validation-review-evidence.tar.gz"
tar -tzf "$study/artifacts/validation-review-evidence.tar.gz" | head
tar -xzf "$study/artifacts/validation-review-evidence.tar.gz" -C "$study"
(cd "$study" && sha256sum -c artifacts/validation-review-files.sha256)
```

Compare the archive digest to the value above. The extraction creates `validation_review/`; it never targets `run/`. Extract into a fresh directory if `validation_review/` already exists. The index must report 121 `OK` checks. `validation_review/final-validation-evidence/evidence-index.json` covers its original 33-file subtree, but **is not a complete index for this package or a self-contained study reproduction manifest**. Its saved-input and fixture references may point outside this archive. Use `validation-review-files.sha256` for this archive's contents and `validation-review-exclusions.json` for every omitted fixture. The original `artifacts/source-files.sha256` and `artifacts/archive-files.sha256` remain unchanged provenance for the 116 original run files and six archives; their paths refer to files **excluded from this validation-review archive**. Restore those six archives using the existing [run archive instructions](README.md) before any saved-input certification.

## Corrected verifier's source identity

The included `validation_review/simulator-correction-spec.md` has SHA256 `dd865b1452a6309ea77439e4cdaeb5842aa03bedaa350111423aeb07cd097bd1`. The corrected verifier's current source identity reads it at the **original runtime path** `/tmp/emc-part2-workspace/writing/runs/REQ-20261027-phase-sensitive/simulator-correction-spec.md`. A package restore elsewhere alone does not satisfy that source identity. Keep the already-existing original file there, or restore this exact copy to that explicit path and confirm its SHA256 before running the corrected CLI. Do not change the frozen source hash or backdate a freeze to accommodate relocation.

With the original study `run/` restored at its recorded checkout path, the historical environment and pinned dependencies present, and that correction specification at its required path, the draw-free command is:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/tmp/bif-v3-world1-validation-env-20261023/bin/python \
  research/emc-part2b/phase_sensitive/measurement_experiment_corrected.py \
  --input-attempt research/emc-part2b/phase_sensitive/run/attempts/correction-002 \
  --evidence-output /tmp/emc-phase-postexecution-validation-fresh \
  --verify-saved
```

Use a fresh evidence-output directory. This post-execution verification generates no new observations and writes no files into the immutable historical `run/` tree. Its certificate has a later source identity and timestamp; it does not retroactively certify a pre-draw step. The missing prospective detailed-key inventory remains a disclosed historical process deviation in the final scientific re-review.
