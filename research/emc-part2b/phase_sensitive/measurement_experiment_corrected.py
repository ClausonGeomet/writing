#!/usr/bin/env python3
"""Draw-free, post-execution validation of the frozen phase study.

The historical attempt is an input. This module never calls its generation path.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import os
from pathlib import Path
import platform
import sys
import time

import numpy as np
import scipy

import measurement_experiment as old
import phase_models as model

HERE = Path(__file__).resolve().parent
MAX_SECONDS = 5400
MAX_BYTES = 500_000_000
METADATA = frozenset(("cell_id", "cohort_id", "support_group", "stage", "budget",
                      "stratum", "arm", "status", "reason", "n_planned", "n_success", "n_failure"))
NUMERIC_TOL = 1e-12


def insist(condition, reason):
    if not condition:
        raise RuntimeError(reason)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for part in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(part)
    return h.hexdigest()


def atomic_json(path, value):
    path = Path(path)
    insist(not path.exists(), f"evidence overwrite refused: {path}")
    data = (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    insist(not tmp.exists(), "temporary evidence path occupied")
    with tmp.open("xb") as stream:
        stream.write(data)
    tmp.replace(path)


def file_bytes(root):
    root = Path(root)
    insist(root.is_dir() and not root.is_symlink(), "resource tree missing or linked")
    total = 0
    for base, dirs, names in os.walk(root, followlinks=False):
        for name in dirs + names:
            path = Path(base) / name
            insist(not path.is_symlink(), f"linked file/tree refused: {path}")
            if path.is_file():
                size = path.stat().st_size
                insist(isinstance(size, int) and size >= 0, "invalid resource measurement")
                total += size
    return total


class ResourceGate:
    def __init__(self, output, *, clock=time.monotonic, seconds=MAX_SECONDS,
                 byte_cap=MAX_BYTES):
        self.output = Path(output)
        self.clock = clock
        self.deadline = clock() + seconds
        self.byte_cap = byte_cap

    def check(self):
        now = self.clock()
        insist(np.isfinite(now) and now < self.deadline, "wall-time cap reached")
        size = file_bytes(self.output)
        insist(size <= self.byte_cap, "logical output byte cap reached")
        return self.deadline - now

    def write_json(self, name, value):
        self.check()
        path = self.output / name
        payload = (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()
        insist(file_bytes(self.output) + len(payload) <= self.byte_cap,
               "evidence write would cross byte cap")
        atomic_json(path, value)
        self.check()


def non_generating_boundary(output, callback, *, preflight, clock=time.monotonic,
                            seconds=MAX_SECONDS, byte_cap=MAX_BYTES, workers=1):
    """Inject a harmless callback; source-level pre-invocation repair only."""
    insist(callable(callback) and callable(preflight), "callback and preflight required")
    insist(type(workers) is int and 1 <= workers <= 2, "two-worker cap")
    insist(all(os.environ.get(k) == "1" for k in
               ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")),
           "one BLAS thread per worker required")
    gate = ResourceGate(output, clock=clock, seconds=seconds, byte_cap=byte_cap)
    gate.check()
    preflight(gate)
    gate.check()  # after all preflight writes, at the actual invocation boundary
    result = callback()
    gate.check()
    return result


def masked_record(keys, rows, n, reason="matched arm failed or missing"):
    """Only accounting fields survive a failed matched block."""
    success = sum(row.get("status") == "ok" for row in rows)
    failure = n - success  # Required rows absent from the block are failures too.
    return {**keys, "status": "unavailable_matched_failure", "reason": reason,
            "n_planned": n, "n_success": success, "n_failure": failure}


def score_block_corrected(rows, truth, observed, bounds):
    """Independent formulas for an entirely successful arm block."""
    n = len(rows)
    if not n:
        return {"status": "empty_stratum", "reason": "no purchased specimens",
                "n_planned": 0, "n_success": 0, "n_failure": 0}
    insist(all(row["status"] == "ok" for row in rows), "partial block passed to scorer")
    estimates = np.array([[float(row[f"{phase}_hat"]) for phase in "QAR"] for row in rows])
    truths = np.array([[float(truth[row["sample_id"]][phase]) for phase in "QAR"] for row in rows])
    bulk = np.array([[float(observed[row["sample_id"]][assay]) for assay in ("SiO2", "TiO2")]
                     for row in rows])
    error = estimates - truths
    result = {"status": "ok", "reason": "", "n_planned": n, "n_success": n, "n_failure": 0,
              "mass_RMSE": float(np.sqrt(np.mean(error ** 2))),
              "totalTi_RMSE": float(np.sqrt(np.mean(np.sum(error[:, 1:3], axis=1) ** 2))),
              "bulk_assay_RMSE_observed": float(np.sqrt(np.mean((estimates @ old.C.T - bulk) ** 2))),
              "bulk_assay_RMSE_latent": float(np.sqrt(np.mean(
                  (estimates @ old.C.T - truths @ old.C.T) ** 2)))}
    for j, phase in enumerate("QAR"):
        result[f"phase_RMSE_{phase}"] = float(np.sqrt(np.mean(error[:, j] ** 2)))
        result[f"phase_bias_{phase}"] = float(np.mean(error[:, j]))
    residuals = [row for row in rows if row["measurement_design_n"] not in ("", None)]
    result["measurement_design_RMSE_observed"] = (
        float(np.sqrt(sum(float(row["measurement_design_sse"]) for row in residuals) /
                      sum(int(row["measurement_design_n"]) for row in residuals)))
        if residuals else None)
    chosen = {row["observation_design_id"] for row in rows}
    relevant = [bound for bound in bounds if bound["observation_design_id"] in chosen]
    valid = [bound for bound in relevant if bound["status"] == "ok"]
    result["n_envelope_nonempty"] = len(valid)
    result["n_envelope_infeasible"] = sum(bound["status"] == "infeasible_envelope"
                                          for bound in relevant)
    result["n_envelope_numerical_failure"] = sum(bound["status"] not in ("ok", "infeasible_envelope")
                                               for bound in relevant)
    for phase in "QAR":
        result[f"mean_width_{phase}"] = (float(np.mean([float(bound[f"width_{phase}"])
                                                     for bound in valid])) if valid else None)
    return result


def numeric_equal(actual, expected, *, tolerance=NUMERIC_TOL):
    if actual is None or expected is None:
        return actual is None and expected is None
    if isinstance(expected, (list, tuple)):
        return isinstance(actual, (list, tuple)) and len(actual) == len(expected) and all(
            numeric_equal(a, b, tolerance=tolerance) for a, b in zip(actual, expected))
    if isinstance(expected, dict):
        return isinstance(actual, dict) and actual.keys() == expected.keys() and all(
            numeric_equal(actual[k], v, tolerance=tolerance) for k, v in expected.items())
    if isinstance(expected, (int, float)) and not isinstance(expected, bool):
        try:
            return np.isfinite(float(actual)) and abs(float(actual) - float(expected)) <= tolerance
        except (TypeError, ValueError):
            return False
    return actual == expected


def csv_equal(saved, expected, *, fields=None):
    fields = fields or expected.keys()
    for key in fields:
        value = saved.get(key, "")
        reference = expected.get(key)
        if reference is None:
            insist(value in ("", None), f"{key}: expected null")
        elif isinstance(reference, (int, float)) and not isinstance(reference, bool):
            insist(numeric_equal(value, reference), f"{key}: numerical discrepancy")
        else:
            insist(str(value) == str(reference), f"{key}: value discrepancy")


def kkt_certificate(m, H, f, support):
    """Simplex KKT with hard zeros treated as equalities."""
    m = np.asarray(m, float)
    allowed = np.asarray(support, bool)
    insist(m.shape == (3,) and np.isfinite(m).all(), "malformed candidate mass")
    feasibility = max(abs(np.sum(m) - 1), max(0., -np.min(m)),
                      np.max(np.abs(m[~allowed])) if (~allowed).any() else 0.)
    insist(feasibility <= 1e-9, "candidate infeasible")
    g = H @ m - f
    free = np.flatnonzero(allowed)
    positive = [j for j in free if m[j] > 1e-9]
    insist(positive, "candidate has no positive coordinate")
    multiplier = float(np.mean(g[positive]))
    stationarity = max(abs(g[j] - multiplier) for j in positive)
    inactive = max((max(0., multiplier - g[j]) for j in free if j not in positive), default=0.)
    complementarity = max((abs(m[j] * max(0., g[j] - multiplier)) for j in free), default=0.)
    scale = max(1., np.linalg.norm(H, ord=np.inf) * max(1., np.linalg.norm(m, ord=np.inf)) +
                np.linalg.norm(f, ord=np.inf))
    residual = max(stationarity, inactive, complementarity) / scale
    insist(residual <= 1e-8, "independent KKT certificate failed")
    return feasibility, residual


def estimate_certificate(row, D, y, R, prior, support, strength):
    insist(row["status"] == "ok", "required estimate failed")
    mass = old.array(row, ("Q_hat", "A_hat", "R_hat"))
    method = row["arm"][:2]
    if method == "A1":
        insist(prior is not None, "A1 prior missing")
        insist(np.max(np.abs(mass - np.asarray(prior["center"]))) <= 1e-8, "A1 center mismatch")
        insist(abs(np.sum(mass) - 1) <= 1e-9 and np.min(mass) >= -1e-9 and
               (not np.any(~np.asarray(support, bool)) or
                np.max(np.abs(mass[~np.asarray(support, bool)])) <= 1e-9),
               "A1 mass infeasible")
        for key in ("assay_loss", "kkt_residual", "measurement_design_sse", "measurement_design_n"):
            insist(row[key] == "", f"A1 {key} must be null")
        for key in ("prior_loss", "objective"):
            insist(numeric_equal(row[key], 0.), f"A1 {key} must be zero")
        insist(float(row["constraint_violation"]) <= 1e-9, "A1 feasibility diagnostic")
        return
    inverse = np.linalg.inv(R)
    K = np.zeros((3, 3)) if method == "A0" else np.asarray(prior["precision"], float)
    mu = np.zeros(3) if method == "A0" else np.asarray(prior["center"], float)
    lam = 0. if method == "A0" else strength
    H = D.T @ inverse @ D + lam * K
    f = D.T @ inverse @ y + lam * K @ mu
    feasibility, residual = kkt_certificate(mass, H, f, support)
    signed = D @ mass - y
    assay = float(signed @ inverse @ signed / 2)
    penalty = float(lam * ((mass - mu) @ K @ (mass - mu)) / 2)
    expected = {"assay_loss": assay, "prior_loss": penalty, "objective": assay + penalty,
                "measurement_design_sse": float(signed @ signed),
                "measurement_design_n": len(y)}
    for key, value in expected.items():
        insist(numeric_equal(row[key], value), f"{key} mismatches actual candidate math")
    insist(float(row["constraint_violation"]) <= 1e-9 and
           float(row["kkt_residual"]) <= 1e-8, "reported diagnostics outside gate")
    for key, value in zip(("bulk_SiO2_hat", "bulk_TiO2_hat"), old.C @ mass):
        insist(numeric_equal(row[key], value), f"{key} mismatch")
    replay = model.estimate(y, D, R, None if method == "A0" else prior,
                            method=method, strength=strength, support=support)
    insist(replay["status"] == "ok" and
           np.max(np.abs(mass - np.asarray(replay["mass"]))) <= 1e-8,
           "minimum-norm replay mass mismatch")
    return feasibility, residual


def safe_paths(input_attempt, evidence_output):
    source = HERE.resolve()
    inp = Path(input_attempt).resolve(strict=True)
    out = Path(evidence_output).resolve(strict=False)
    historical = Path("/home/matt/Documents/GitHub/writing/research/emc-part2b/phase_sensitive/run").resolve()
    insist(inp == historical / "attempts" / "correction-002" and inp.is_dir(),
           "complete historical correction-002 input required")
    for forbidden in (historical, source, inp):
        insist(out != forbidden and not out.is_relative_to(forbidden) and
               not forbidden.is_relative_to(out), "input/output or source/output overlap")
    insist(not out.exists(), "fresh evidence destination required")
    return inp, out


def historical_gate(inp):
    index = HERE / "artifacts" / "source-files.sha256"
    expected = {}
    for line in index.read_text().splitlines():
        hash_value, rel = line.split("  ", 1)
        expected[rel] = hash_value
    root = inp.parents[1]
    for rel, hash_value in expected.items():
        path = root.parent / rel
        insist(path.is_file() and not path.is_symlink() and sha(path) == hash_value,
               f"historical archive inventory mismatch: {rel}")
    execution = old.read(inp / "manifests/execution.json")
    insist(execution["status"] == "complete" and execution["schema"] == "qar-phase-execution-v1",
           "historical attempt incomplete")
    cfg, review, freeze = old.design_gate(HERE / "preregistration.json", HERE / "PREEXECUTION_REVIEW.json")
    manifests = inp / "manifests"
    frozen = old.read(manifests / "freeze.json")
    insist(frozen["design_digest"] == freeze["design_digest"] and
           frozen["review_sha256"] == sha(HERE / "PREEXECUTION_REVIEW.json") and
           frozen["source_sha256"] == sha(manifests / "source.json") and
           frozen["environment_sha256"] == sha(manifests / "environment.json") and
           frozen["seeds_sha256"] == sha(manifests / "seeds.json"), "historical freeze chain")
    insist(old.read(manifests / "source.json") == old.complete_source_manifest(cfg),
           "historical five-source/pinned-package identity mismatch")
    saved_seeds = old.read(manifests / "seeds.json")
    expected_seeds = old.seed_manifest()
    insist(saved_seeds["streams"] == expected_seeds["streams"] and
           saved_seeds["historical_seed_set"] == expected_seeds["historical_seed_set"] and
           saved_seeds["collision_count"] == 0, "historical seed schedule mismatch")
    insist(len(saved_seeds["historical_manifests"]) == len(expected_seeds["historical_manifests"]),
           "historical seed inventory count mismatch")
    for saved, current in zip(saved_seeds["historical_manifests"],
                              expected_seeds["historical_manifests"]):
        insist(saved["sha256"] == current["sha256"] and
               saved["explicit_seed_values"] == current["explicit_seed_values"] and
               Path(saved["path"]).is_file() and sha(saved["path"]) == saved["sha256"],
               "historical seed inventory file mismatch")
    insist(execution["freeze_sha256"] == sha(inp / "manifests/freeze.json") and
           execution["review_sha256"] == sha(HERE / "PREEXECUTION_REVIEW.json"),
           "historical freeze/review chain mismatch")
    for rel, value in execution["output_hashes"].items():
        insist(sha(inp / rel) == value, f"historical output hash mismatch: {rel}")
    access = old.read(manifests / "access.json")
    insist(access["truth_opened_only_after_prediction_sha256"] ==
           sha(inp / "results/predictions.csv"), "historical prediction/access lineage mismatch")
    for rel, value in access["data_sha256"].items():
        insist(sha(inp / rel) == value, f"historical data access hash mismatch: {rel}")
    return cfg, execution


def source_record(inp):
    sources = (HERE / "measurement_experiment_corrected.py", HERE / "test_experiment_corrected.py")
    reused = old.source_paths()
    return {"schema": "qar-phase-postexecution-source-v1",
            "purpose": "post-execution saved-input validation",
            "recorded_utc": old.utc(),
            "additive_sources": {str(p): sha(p) for p in sources},
            "historical_sources": {str(p): sha(p) for p in reused},
            "correction_spec_sha256": sha(Path("/tmp/emc-part2-workspace/writing/runs/REQ-20261027-phase-sensitive/simulator-correction-spec.md")),
            "historical_manifests": {p.name: sha(p) for p in (inp / "manifests").glob("*.json")},
            "saved_inputs_results": {str(p.relative_to(inp)): sha(p) for section in ("data", "results")
                                     for p in (inp / section).glob("*") if p.is_file()},
            "validation_environment": {"python": sys.version, "executable": sys.executable,
                                       "numpy": np.__version__, "scipy": scipy.__version__,
                                       "platform": platform.platform(),
                                       "thread_environment": {k: os.environ.get(k) for k in
                                           ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")}},
            "historical_environment_sha256": sha(inp / "manifests/environment.json"),
            "design_digest": "7f011812d4ecfee33f8184bab59ea0591d98580215782d9dc509e547bfacabad"}


def cohort_key(row):
    return int(row["cell_id"]), int(row["cohort_id"])


def records(path, *, jsonlines=False):
    with Path(path).open(newline="") as stream:
        rows = (json.loads(line) for line in stream) if jsonlines else csv.DictReader(stream)
        yield from rows


def grouped(path, *, jsonlines=False):
    for key, group in itertools.groupby(records(path, jsonlines=jsonlines), cohort_key):
        yield key, list(group)


def group_at(iterator, key):
    actual, rows = next(iterator)
    insist(actual == key, f"saved cohort order/key mismatch: {actual} != {key}")
    return rows


def expected_priors(cell, cohort, observed, saved):
    s, r = cohort
    pool = observed["pool"]
    cal = observed["calibration"]
    insist(len(cal) == 20 and len(pool) == 10, "calibration/pool cohort size")
    cal_y = np.stack([old.array(x, ("QXRD_Q", "QXRD_A", "QXRD_R")) for x in cal])
    expected = {}
    for group in cell["support_groups"]:
        sup = old.support(group)
        initial_id = f"ps1_s{s}_r{r}_{group}_initial"
        initial = model.fit_initial(cal_y, old.RQ, center_support=tuple(cell["center_support"]))
        expected[initial_id] = {"prior": initial, "update": "initial", "budget": 0,
                                "excluded_unit": None, "support_group": group}
        for budget in (0, 2, 5):
            for mode in ("B", "R", "P"):
                units = []
                for row in pool[:budget]:
                    D, y, cov = old.design(mode, row, True)
                    units.append({"unit_id": row["sample_id"], "D": D, "y": y, "cov": cov})
                for excluded in [None] + [unit["unit_id"] for unit in units]:
                    pid = f"ps1_s{s}_r{r}_{group}_{mode}_k{budget}_{excluded or 'all'}"
                    prior = model.update_mean(initial, units, support=sup, exclude_unit=excluded)
                    expected[pid] = {"prior": prior, "update": mode, "budget": budget,
                                     "excluded_unit": excluded, "support_group": group}
    insist(len(saved) == len(expected) and len({x["prior_id"] for x in saved}) == len(saved),
           "prior key count or uniqueness")
    for row in saved:
        pid = row["prior_id"]
        insist(pid in expected, "unexpected prior ID")
        e = expected[pid]
        for key in ("update", "budget", "excluded_unit", "support_group"):
            insist(row[key] == e[key], f"prior provenance mismatch: {pid}/{key}")
        insist(numeric_equal(row["prior"], e["prior"]), f"prior fitted arrays/diagnostics mismatch: {pid}")
        insist(row["prior"]["status"] == "ok", f"failed required prior: {pid}")
    return {pid: item["prior"] for pid, item in expected.items()}


def expected_prediction(cell, cohort, row, observed, priors):
    s, r = cohort
    group = row["support_group"]
    insist(group in cell["support_groups"], "prediction support group")
    stage, budget, stratum, arm = int(row["stage"]), int(row["budget"]), row["stratum"], row["arm"]
    sample_id = row["sample_id"]
    if stage == 1:
        insist(budget == -1 and stratum == "direct" and arm in old.ARMS1, "stage-one key")
        candidates = observed["direct"]
        mode = "B" if arm == "A1" else arm[-1]
        indices = {"B": [0, 1], "T": [0, 1, 2], "H": [0, 1, 3]}[mode]
        sample = next((x for x in candidates if x["sample_id"] == sample_id), None)
        insist(sample is not None, "unselected direct ID")
        D = old.D4[indices]
        y = old.array(sample, ("SiO2", "TiO2", "redundant_Ti", "H_A"))[indices]
        cov = old.R4[np.ix_(indices, indices)]
        pid = f"ps1_s{s}_r{r}_{group}_initial"
    else:
        insist(stage == 2 and budget in (0, 2, 5) and
               stratum in ("measured", "unmeasured") and arm in old.ARMS2,
               "stage-two key")
        candidates = observed["pool"][:budget] if stratum == "measured" else observed["evaluation"]
        sample = next((x for x in candidates if x["sample_id"] == sample_id), None)
        insist(sample is not None, "unselected pool/evaluation ID")
        tag = arm.split("_", 1)[1]
        update = "N" if tag == "N" or tag.endswith("0") else tag
        mode = "B" if stratum == "unmeasured" or tag == "N" else tag[0]
        D, y, cov = old.design(mode, sample, stratum == "measured")
        pid = (f"ps1_s{s}_r{r}_{group}_initial" if update == "N" else
               f"ps1_s{s}_r{r}_{group}_{update}_k{budget}_" +
               (sample_id if stratum == "measured" else "all"))
    did = f"{mode}_{sample_id}_{group}"
    method = arm[:2]
    actual_pid = "" if method == "A0" else pid
    insist(row["prior_id"] == actual_pid and row["observation_design_id"] == did,
           "prediction prior/design identity or whole-unit exclusion mismatch")
    prior = None if method == "A0" else priors[pid]
    return sample, D, y, cov, prior


def expected_prediction_keys(cell, cohort, observed):
    keys = set()
    for group in cell["support_groups"]:
        for x in observed["direct"]:
            for arm in old.ARMS1:
                keys.add((group, 1, -1, "direct", arm, x["sample_id"]))
        for budget in (0, 2, 5):
            for stratum, samples in (("measured", observed["pool"][:budget]),
                                     ("unmeasured", observed["evaluation"])):
                for x in samples:
                    for arm in old.ARMS2:
                        keys.add((group, 2, budget, stratum, arm, x["sample_id"]))
    return keys


def prediction_key(row):
    return (row["support_group"], int(row["stage"]), int(row["budget"]),
            row["stratum"], row["arm"], row["sample_id"])


def verify_predictions(inp, cfg, gate):
    obs_iters = {kind: grouped(inp / "data" / f"{kind}_observed.csv")
                 for kind in ("calibration", "direct", "pool", "evaluation")}
    prior_iter = grouped(inp / "results/priors.jsonl", jsonlines=True)
    pred_iter = grouped(inp / "results/predictions.csv")
    count = 0
    prior_count = 0
    for cell in cfg["cells"]:
        for r in range(1, 51):
            gate.check()
            cohort = cell["cell_id"], r
            observed = {kind: group_at(iterator, cohort) for kind, iterator in obs_iters.items()}
            insist([len(observed[x]) for x in ("calibration", "direct", "pool", "evaluation")] ==
                   [20, 20, 10, 30], "observed schedule row counts")
            for kind, rows in observed.items():
                insist(len({x["sample_id"] for x in rows}) == len(rows), f"duplicate {kind} ID")
            saved_priors = group_at(prior_iter, cohort)
            priors = expected_priors(cell, cohort, observed, saved_priors)
            prior_count += len(saved_priors)
            rows = group_at(pred_iter, cohort)
            keys = expected_prediction_keys(cell, cohort, observed)
            insist(len(rows) == len(keys) and {prediction_key(x) for x in rows} == keys,
                   "prediction full expected key inventory mismatch")
            for row in rows:
                sample, D, y, cov, prior = expected_prediction(cell, cohort, row, observed, priors)
                estimate_certificate(row, D, y, cov, prior, old.support(row["support_group"]),
                                     cell["strength"])
                mass = old.array(row, ("Q_hat", "A_hat", "R_hat"))
                bulk = old.array(sample, ("SiO2", "TiO2"))
                for key, value in zip(("bulk_resid_SiO2", "bulk_resid_TiO2"), old.C @ mass - bulk):
                    insist(numeric_equal(row[key], value), f"{key} mismatch")
            count += len(rows)
    for iterator in (*obs_iters.values(), prior_iter, pred_iter):
        insist(next(iterator, None) is None, "unexpected trailing cohort")
    insist(count == 350250 and prior_count == 7750, "full prediction/prior count")
    return {"estimates": count, "priors": prior_count}


def verify_auxiliary(inp, cfg, gate):
    """Verify saved acquisition, information and LP records from observations."""
    obs = {kind: list(records(inp / "data" / f"{kind}_observed.csv"))
           for kind in ("direct", "pool", "evaluation")}
    by_id = {x["sample_id"]: x for rows in obs.values() for x in rows}
    acquisition = list(records(inp / "results/acquisition.csv"))
    insist(len(acquisition) == 600, "acquisition row count")
    seen_acquisition = set()
    for row in acquisition:
        gate.check()
        s, r, k = int(row["cell_id"]), int(row["cohort_id"]), int(row["budget"])
        key = s, r, k
        insist(key not in seen_acquisition and k in (0, 2, 5), "acquisition key")
        seen_acquisition.add(key)
        pool = [x for x in obs["pool"] if cohort_key(x) == (s, r)]
        csv_equal(row, {"purchased_ids": "|".join(x["sample_id"] for x in pool[:k]),
                        "n_calibration_QXRD": 20, "n_pool_base_bulk": k,
                        "n_replicate_bulk": k, "n_pool_QXRD": k,
                        "coverage_pool": k / 10})
    bounds = list(records(inp / "results/bounds.csv"))
    info = list(records(inp / "results/information.jsonl", jsonlines=True))
    expected_info = {}
    bound_keys = set()
    for row in bounds:
        gate.check()
        sample = by_id.get(row["sample_id"])
        insist(sample is not None and cohort_key(sample) == cohort_key(row), "bound sample provenance")
        group = row["support_group"]
        insist(group in next(x for x in cfg["cells"] if x["cell_id"] == int(row["cell_id"]))["support_groups"],
               "bound support provenance")
        s, k, stratum = int(row["stage"]), int(row["budget"]), row["stratum"]
        did = row["observation_design_id"]
        mode = did.split("_", 1)[0]
        insist(did == f"{mode}_{sample['sample_id']}_{group}", "bound design identity")
        if s == 1:
            insist(k == -1 and stratum == "direct" and sample in obs["direct"] and
                   mode in ("B", "T", "H"), "direct bound key")
            ix = {"B": [0, 1], "T": [0, 1, 2], "H": [0, 1, 3]}[mode]
            D = old.D4[ix]
            y = old.array(sample, ("SiO2", "TiO2", "redundant_Ti", "H_A"))[ix]
            cov = old.R4[np.ix_(ix, ix)]
        else:
            insist(s == 2 and k in (0, 2, 5) and stratum in ("measured", "unmeasured"),
                   "stage-two bound key")
            if stratum == "measured":
                pool = [x for x in obs["pool"] if cohort_key(x) == cohort_key(row)]
                insist(sample in pool[:k] and mode in ("B", "R", "P"), "unselected measured bound")
            else:
                insist(sample in obs["evaluation"] and mode == "B", "unmeasured bound must be bulk")
            D, y, cov = old.design(mode, sample, stratum == "measured")
        key = (cohort_key(row), group, s, k, stratum, sample["sample_id"], mode)
        insist(key not in bound_keys, "duplicate bound key")
        bound_keys.add(key)
        b = model.abundance_bounds(y, D, cov, support=old.support(group))
        insist(row["status"] == b["status"] and row["reason"] == b["reason"], "LP status mismatch")
        insist(b["status"] in ("ok", "infeasible_envelope"), "numerical LP failure retained")
        for prefix, field in (("lo", "lower"), ("hi", "upper"), ("width", "width")):
            for phase, value in zip("QAR", b[field] if b[field] is not None else (None,) * 3):
                csv_equal(row, {f"{prefix}_{phase}": value})
        expected_info.setdefault(did, (D, cov, old.support(group)))
    expected_bounds = set()
    for cell in cfg["cells"]:
        for r in range(1, 51):
            cohort = cell["cell_id"], r
            direct = [x for x in obs["direct"] if cohort_key(x) == cohort]
            pool = [x for x in obs["pool"] if cohort_key(x) == cohort]
            evaluation = [x for x in obs["evaluation"] if cohort_key(x) == cohort]
            for group in cell["support_groups"]:
                for sample in direct:
                    for mode in ("B", "T", "H"):
                        expected_bounds.add((cohort, group, 1, -1, "direct", sample["sample_id"], mode))
                for budget in (0, 2, 5):
                    for sample in evaluation:
                        expected_bounds.add((cohort, group, 2, budget, "unmeasured",
                                             sample["sample_id"], "B"))
                    for sample in pool[:budget]:
                        for mode in ("B", "R", "P"):
                            expected_bounds.add((cohort, group, 2, budget, "measured",
                                                 sample["sample_id"], mode))
    insist(len(bounds) == len(expected_bounds) == 42750 and bound_keys == expected_bounds,
           "bound full expected key inventory mismatch")
    insist(len(info) == len(expected_info), "information design count")
    info_ids = set()
    for row in info:
        did = row["observation_design_id"]
        insist(did in expected_info and did not in info_ids, "information key uniqueness")
        info_ids.add(did)
        D, cov, support = expected_info[did]
        insist(row["support"] == list(support) and numeric_equal(
            row["information"], model.information(D, cov, support=support)),
            "information rank/singular value mismatch")
    return {"acquisition": len(acquisition), "bounds": len(bounds), "information": len(info)}


def score_saved(inp, cfg, gate):
    """Truth becomes readable only after all estimate certificates pass."""
    truth = {}
    observed = {}
    for kind in ("direct", "pool", "evaluation"):
        truth.update({x["sample_id"]: x for x in records(inp / "data" / f"{kind}_truth.csv")})
        observed.update({x["sample_id"]: x for x in records(inp / "data" / f"{kind}_observed.csv")})
    bounds = {}
    for row in records(inp / "results/bounds.csv"):
        key = (int(row["cell_id"]), int(row["cohort_id"]), row["support_group"],
               int(row["stage"]), int(row["budget"]), row["stratum"])
        bounds.setdefault(key, []).append(row)
    metrics = {}
    for row in records(inp / "results/cohort_metrics.csv"):
        key = (int(row["cell_id"]), int(row["cohort_id"]), row["support_group"],
               int(row["stage"]), int(row["budget"]), row["stratum"], row["arm"])
        insist(key not in metrics, "duplicate metric key")
        metrics[key] = row
    insist(len(metrics) == 21250, "metric count")
    pred_iter = grouped(inp / "results/predictions.csv")
    compared = 0
    for cell in cfg["cells"]:
        for r in range(1, 51):
            gate.check()
            cohort = cell["cell_id"], r
            rows = group_at(pred_iter, cohort)
            groups = {}
            for p in rows:
                block = (cohort[0], cohort[1], p["support_group"], int(p["stage"]),
                         int(p["budget"]), p["stratum"], p["arm"])
                groups.setdefault(block, []).append(p)
            for sup in cell["support_groups"]:
                for stage, budget, stratum, arms, n in (
                    (1, -1, "direct", old.ARMS1, 20),
                    *((2, k, st, old.ARMS2, 30 if st == "unmeasured" else k)
                      for k in (0, 2, 5) for st in ("unmeasured", "measured"))):
                    block = (*cohort, sup, stage, budget, stratum)
                    bad = n > 0 and any(len(groups.get((*block, arm), [])) != n or
                                         any(x["status"] != "ok" for x in groups.get((*block, arm), []))
                                         for arm in arms)
                    for arm in arms:
                        key = (*block, arm)
                        selected = groups.get(key, [])
                        if bad:
                            expected = masked_record({}, selected, n)
                        else:
                            expected = score_block_corrected(selected, truth, observed,
                                                             bounds.get(block, []))
                        csv_equal(metrics[key], expected, fields=expected.keys())
                        if bad:
                            for field in metrics[key]:
                                if field not in METADATA:
                                    insist(metrics[key][field] == "", f"failed block derived field leaked: {field}")
                        compared += 1
    insist(next(pred_iter, None) is None and compared == 21250, "metric loop incomplete")
    paired = {}
    for row in records(inp / "results/paired_summary.csv"):
        key = (int(row["cell_id"]), row["support_group"], int(row["stage"]),
               int(row["budget"]), row["stratum"], row["contrast"], row["metric"])
        insist(key not in paired, "duplicate paired key")
        paired[key] = row
    expected_paired = 0
    for cell in cfg["cells"]:
        s = cell["cell_id"]
        for sup in cell["support_groups"]:
            for stage, budget, stratum, contrasts in (
                (1, -1, "direct", old.CONTRASTS1),
                *((2, k, st, old.CONTRASTS2) for k in (0, 2, 5)
                  for st in ("unmeasured", "measured"))):
                for contrast in contrasts:
                    a, b = contrast.split("-")
                    for metric in old.METRICS:
                        key = s, sup, stage, budget, stratum, contrast, metric
                        values = []
                        for r in range(1, 51):
                            x = metrics[(s, r, sup, stage, budget, stratum, a)]
                            y = metrics[(s, r, sup, stage, budget, stratum, b)]
                            if x["status"] == y["status"] == "ok":
                                values.append(float(x[metric]) - float(y[metric]))
                        expected = {"n_planned": 50, "n_complete": len(values),
                                    "n_failed": 50 - len(values),
                                    "mean_difference": float(np.mean(values)) if len(values) == 50 else None,
                                    "descriptive_se": float(np.std(values, ddof=1) / np.sqrt(50))
                                    if len(values) == 50 else None}
                        csv_equal(paired[key], expected)
                        expected_paired += 1
    insist(len(paired) == expected_paired == 3050, "paired full key count")
    return {"cohort_metrics": compared, "paired_summary": expected_paired}


def verify_saved(input_attempt, evidence_output):
    inp, out = safe_paths(input_attempt, evidence_output)
    out.mkdir(parents=True, exist_ok=False)
    gate = ResourceGate(out)
    gate.check()
    insist(all(os.environ.get(k) == "1" for k in
               ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")),
           "one BLAS thread required for validation")
    record = source_record(inp)
    gate.write_json("source.json", record)
    evidence = {"schema": "qar-phase-postexecution-evidence-v1",
                "purpose": "post-execution saved-input validation", "started_utc": old.utc(),
                "command": [sys.executable, str(Path(__file__).resolve()), "--input-attempt",
                            str(inp), "--evidence-output", str(out), "--verify-saved"],
                "historical_execution_utc": old.read(inp / "manifests/execution.json").get("ended_utc"),
                "source_sha256": sha(out / "source.json"), "input_attempt": str(inp),
                "generator_calls": 0, "draw_count": 0, "verdict": "INCOMPLETE",
                "process_deviation": "Detailed prospective key inventory was absent before draw; full keys reconstructed retrospectively.",
                "tolerances": {"saved_numeric_atol": 1e-12, "saved_numeric_rtol": 0,
                               "mass_oracle_atol": 1e-8, "feasibility": 1e-9,
                               "normalized_kkt": 1e-8}}
    gate.write_json("start.json", {"schema": "qar-phase-postexecution-start-v1",
                                  "source_sha256": evidence["source_sha256"],
                                  "started_utc": evidence["started_utc"],
                                  "purpose": evidence["purpose"],
                                  "verdict": "INCOMPLETE", "generator_calls": 0,
                                  "draw_count": 0})
    try:
        cfg, execution = historical_gate(inp)
        gate.check()
        evidence["historical_execution_sha256"] = sha(inp / "manifests/execution.json")
        evidence["counts"] = verify_predictions(inp, cfg, gate)
        evidence["counts"].update(verify_auxiliary(inp, cfg, gate))
        evidence["counts"].update(score_saved(inp, cfg, gate))
        diagnostic = old.study_diagnostics(inp, cfg, execution["counts"],
                                           sha(inp / "results/predictions.csv"))
        insist(numeric_equal(old.read(inp / "results/diagnostics.json"), diagnostic),
               "saved study diagnostics mismatch")
        evidence["counts"]["diagnostics"] = 1
        insist(source_record_identity(record, inp), "validation source/input changed during operation")
        evidence["verdict"] = "PASS"
    except Exception as exc:
        evidence["verdict"] = "BLOCK"
        evidence["reason"] = str(exc)
        raise
    finally:
        evidence["ended_utc"] = old.utc()
        try:
            gate.write_json("validation.json", evidence)
            gate.write_json("evidence-hashes.json", {"schema": "qar-phase-evidence-hashes-v1",
                "source_sha256": sha(out / "source.json"),
                "start_sha256": sha(out / "start.json"),
                "validation_sha256": sha(out / "validation.json")})
        except Exception:
            # Preserve a reason when an exhausted resource gate prevents the normal record.
            if not (out / "validation.json").exists():
                fallback = {**evidence, "verdict": "INCOMPLETE",
                            "reason": evidence.get("reason", "resource gate blocked evidence")}
                raw = (json.dumps(fallback, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()
                if file_bytes(out) + len(raw) <= gate.byte_cap:
                    atomic_json(out / "validation.json", fallback)
    return evidence


def source_record_identity(record, inp):
    current = source_record(inp)
    for field in ("additive_sources", "historical_sources", "correction_spec_sha256",
                  "historical_manifests", "saved_inputs_results", "historical_environment_sha256"):
        if current[field] != record[field]:
            return False
    return True


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-attempt", type=Path, required=True)
    parser.add_argument("--evidence-output", type=Path, required=True)
    parser.add_argument("--verify-saved", action="store_true", required=True)
    args = parser.parse_args(argv)
    try:
        result = verify_saved(args.input_attempt, args.evidence_output)
        print(json.dumps(result, sort_keys=True))
        return 0
    except Exception as exc:
        print(f"BLOCK: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())


