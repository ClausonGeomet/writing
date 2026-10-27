# Reproducing and inspecting the bounded study

The audited output was generated in the **main checkout** at `research/emc-part2b/phase_sensitive/run/attempts/correction-002/`. It was uncommitted when the initial report was prepared. The complete `run/` tree, including earlier failed attempts, is now preserved in committed [archives with a restoration and byte-check recipe](artifacts/README.md). Restore those archives before following the saved-run instructions below; no raw run files are staged directly in Git. Preserve the root failed run and `attempts/correction-001/` byte-for-byte. The eight correction-002 data CSVs are byte-identical to correction-001; compare `audit/correction001_data_before.json` with `audit/data_replay_hashes.json`. `manifests/execution.json` lists SHA256 for every generated output; `manifests/correction.json`, `source.json`, `environment.json`, `seeds.json`, `freeze.json`, `design_freeze.json`, `deviations.json`, and `preservation.json` bind the lineage. `audit/postrun_audit.json`, `postchecks_summary.json`, `selective_replay.json`, and the runtime `audit.md` carry the original numerical audit. [POSTEXECUTION_VALIDATION.md](POSTEXECUTION_VALIDATION.md) records the later additive validation lineage.

## Read-only inspection of the saved run

From the main checkout root, set `STUDY=research/emc-part2b/phase_sensitive` and `ATTEMPT=$STUDY/run/attempts/correction-002`. Inspect `$ATTEMPT/results/paired_summary.csv` for all preregistered cohort contrasts, `$ATTEMPT/results/cohort_metrics.csv` for phase/assay/totalTi errors, `$ATTEMPT/results/bounds.csv` and `information.jsonl` for feasible sets and rank, and `$ATTEMPT/manifests/execution.json` for counts/hashes. The primary rows have `stage=2`, `budget=5`, `stratum=unmeasured`, `contrast=A2_P-A2_R`, `metric=mass_RMSE`; keep cell/support separate. `mean_difference` is the mean of 50 cohort differences, and `descriptive_se` is their sample SD divided by √50. `--verify` reads existing outputs without a new draw, but it also checks the frozen environment and source; use the original environment described below. Do not rerun `--run` as a routine report check.

The exact final public command, from the main checkout root, was:

```sh
export PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR=/tmp/emc-phase-mpl
/tmp/bif-v3-world1-validation-env-20261023/bin/python research/emc-part2b/phase_sensitive/measurement_experiment.py \
  --config research/emc-part2b/phase_sensitive/preregistration.json \
  --output research/emc-part2b/phase_sensitive/run \
  --review research/emc-part2b/phase_sensitive/PREEXECUTION_REVIEW.json \
  --attempt correction-002 --verify
```

The tester used the same base command with `--check`, then `--correction-reason '<recorded source-only reason>' --freeze`, then `--run`, and finally `--verify`; only `--run` generated observations. Historical `--verify` checks integrity under executed sources; it is weaker than the later independent mathematical certificate. The final freeze preceded the first draw. The original failed root and correction-001 attempts have their own freeze/execution records; they are part of reproduction history, not alternative complete results. Exact source SHA256s are in correction-002 `manifests/source.json`; notably `generation.R=16432637b16c4a72eb07ec340d08cc96910d2d3172285d503b18e96428066610`, `measurement_experiment.py=9577c4799d5ff8fbb8bbd14098cac23cafe8567fab91521ec46f0320c092d222`, and `phase_models.py=2910f0f43c7effb392a8d71edf64c87e4b314712774b076738b127d324ead2e1`. Frozen design digest is `7f011812d4ecfee33f8184bab59ea0591d98580215782d9dc509e547bfacabad`.

The six frozen source SHA256s (study-relative) are: `EXPERIMENT_PLAN.md 9f3cbe1e916261789232448d5d474fe9a9f550cb2aac16b6ec593fd65ee2054d`; `LITERATURE_REVIEW.md 49692fd953093789f1c714037a8645582d710ac831624bb1fcf7cb084eea8c3b`; `MATHEMATICAL_DESIGN.md 6c4d09ac7962a0611c4c8366c01cf66590ab299e77225ab32400b676c293558c`; `PREREGISTRATION.md 24137254834f434978b6f4fc7b7740c39b9ab47d954128c75a925d0ad3b63ac1`; `RESEARCH_DECISIONS.md f573570271ca52395a9b693665847e86cd74357e856f5f15891b5afb030b60d2`; `preregistration.json 8fd9e88fde929b4ca6b74a7ce15c81a6f8e1d300182446d274301d29151907f4`. `design_freeze.json` gives the digest rule and review linkage. New results prose and figures are outside this frozen inventory.

## Environment and fresh generation limits

The frozen environment used Python 3.12.3, NumPy 2.5.3, SciPy 1.18.1, matplotlib 3.11.2, pytest 9.1.1, R 4.3.3, jsonlite 2.0.0, and simCoDaR 0.1.0.9000 from public commit `d9a0b27676a0cb3a4acb300a3f6b2ebc845d205d`; exact R/BLAS/installation provenance and thread settings are in `manifests/environment.json` and `source.json`. The historical-seed gate inventories 11 earlier manifest paths and 2,009 historical integers as well as prior master/search seeds. Its path hashes and collision proof are in runtime `historical-seed-inventory.json` and correction-002 `manifests/seeds.json` (SHA256 `c759ec08ef7b2f7bb6ce497000dfabe783ec1d3c94ee880baadb16d2dfd1a195`). A clean fresh clone alone does **not** recreate that historical runtime inventory or the old blocked attempts. Fresh generation requires reconstructing and independently certifying those legacy inputs and the pinned environment, retaining prior attempts, then using the same unchanged 2,000 seed schedule. Do not claim a clean-clone one-command rerun.

## Redraw the two report figures without generating observations

After [restoring and checking the archives](artifacts/README.md), run this from the repository root. The input is the saved correction-002 summary; choose an **explicit, separate output directory** so the tracked report figures and run evidence are not overwritten:

```sh
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/emc-phase-mpl \
/tmp/bif-v3-world1-validation-env-20261023/bin/python \
  research/emc-part2b/phase_sensitive/render_figures.py \
  --paired-summary research/emc-part2b/phase_sensitive/run/attempts/correction-002/results/paired_summary.csv \
  --output-dir /tmp/emc-phase-report-figures
```

This writes `transfer_summary.png` and `direct_vs_transfer.png` into the chosen directory using only `paired_summary.csv`. The script is a report renderer outside the frozen six-file design and five computational sources; it does not run the estimator, generate observations, or change saved evidence. The original PNGs remain in `figures/`. Matplotlib renderer or metadata differences may change PNG bytes across environments without changing the plotted summaries; use the pinned environment above for an exact comparison. In the pinned environment this renderer produced PNGs byte-identical to the tracked figures. Numerical audit PASS and historical pre-execution PASS remain distinct from the initial final scientific **STOP**. R2 reproduction is addressed; R1 additive saved-result certification now passes for `807a0b8`; independent scientific re-review remains pending.

## Additive post-execution check

From the original main-checkout root `/home/matt/Documents/GitHub/writing`, after restoring the saved attempt, run this draw-free command with a fresh evidence directory:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 \
/tmp/bif-v3-world1-validation-env-20261023/bin/python \
  research/emc-part2b/phase_sensitive/measurement_experiment_corrected.py \
  --input-attempt /home/matt/Documents/GitHub/writing/research/emc-part2b/phase_sensitive/run/attempts/correction-002 \
  --evidence-output /tmp/emc-phase-postexecution-validation-fresh \
  --verify-saved
```

The output path must be unused and separate from `run/` and the source tree. [POSTEXECUTION_VALIDATION.md](POSTEXECUTION_VALIDATION.md) explains its scope. The [archive recipe](artifacts/README.md) restores bytes but does not supply the original absolute checkout, recorded environment, pinned simCoDaR package, historical seed specifications or runtime `/tmp/emc-part2-workspace/writing/runs/REQ-20261027-phase-sensitive/simulator-correction-spec.md` required for source identity. The independent tester issued a fresh draw-free certificate for `807a0b8`: 33 deterministic tests passed; 350,250 saved estimates, 7,750 priors, 42,750 bounds, 26,250 information records, 600 acquisition rows, 21,250 cohort metrics and 3,050 paired records certified; zero draws and generator calls. Its source-record SHA256 is `e382c97d4fcfe5a93141f574f2a75521cd22d8db655b3485f3deab1724d01da2`, and validation-record SHA256 is `79a3814c5c5c447e46c73b8c53536f96bdce71577ba33fdb70e871ba3cf0b316`. All 116 original run files and six archives remained byte-identical. Independent scientific re-review still decides whether the existing STOP can be lifted. The prospective detailed-key inventory was absent before draw; any full comparison now is retrospective.
