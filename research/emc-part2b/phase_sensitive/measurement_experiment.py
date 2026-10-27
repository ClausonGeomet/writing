#!/usr/bin/env python3
"""Frozen, gated QAR phase measurement study. Only --run invokes generation.R."""
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import time
from datetime import datetime, timezone

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PIN = "d9a0b27676a0cb3a4acb300a3f6b2ebc845d205d"
ROLES = ("calibration_truth", "direct_truth", "pool_truth", "evaluation_truth",
         "calibration_noise", "direct_noise", "pool_bulk_noise", "evaluation_bulk_noise",
         "pool_replicate_noise", "pool_qxrd_noise")
ARMS1 = ("A0_B", "A1", "A2_B", "A0_T", "A2_T", "A0_H", "A2_H")
ARMS2 = ("A0_B", "A0_R", "A0_P", "A1_N", "A2_N", "A1_B", "A2_B",
         "A1_R", "A2_R", "A1_P", "A2_P", "A2_R0", "A2_P0")
CONTRASTS1 = ("A0_H-A0_T", "A0_H-A0_B", "A2_H-A2_T", "A2_H-A2_B",
              "A2_H-A0_H", "A2_B-A0_B", "A2_B-A1")
CONTRASTS2 = ("A2_P-A2_R", "A2_P-A2_N", "A2_P-A2_B", "A2_P-A1_P",
              "A2_P-A0_P", "A2_R-A2_B", "A2_P-A2_P0", "A2_R-A2_R0", "A1_P-A1_N")
METRICS = ("mass_RMSE", "phase_RMSE_Q", "phase_RMSE_A", "phase_RMSE_R",
           "phase_bias_Q", "phase_bias_A", "phase_bias_R", "totalTi_RMSE",
           "bulk_assay_RMSE_observed", "bulk_assay_RMSE_latent")
C = np.array([[1., 0., 0.], [0., 1., 1.]])
H = np.array([[0., 1., 0.]])
RB = np.array([[1e-4, 3e-5], [3e-5, 1e-4]])
RQ = np.eye(3) * 4e-4
SD4 = np.diag([.01, .01, .01, .02])
CORR4 = np.array([[1., .3, 0., 0.], [.3, 1., .5, 0.],
                  [0., .5, 1., 0.], [0., 0., 0., 1.]])
R4 = SD4 @ CORR4 @ SD4
D4 = np.vstack([C, C[1], H])


def utc():
    return datetime.now(timezone.utc).isoformat()


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n")
    temp.replace(path)


def read(path):
    return json.loads(Path(path).read_text())


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def design_gate(config, review):
    require(config.resolve() == HERE / "preregistration.json", "config must be public frozen preregistration")
    require(review.resolve() == HERE / "PREEXECUTION_REVIEW.json", "review must be public gate")
    f = read(HERE / "design_freeze.json")
    q = read(review)
    files = {p: digest(HERE / p) for p in f["design_files"]}
    require(files == f["design_files"], "frozen design byte mismatch")
    lines = "".join(f"{p} {files[p]}\n" for p in sorted(files))
    actual = hashlib.sha256(lines.encode("utf-8")).hexdigest()
    require(actual == f["design_digest"], "actual-LF design digest mismatch")
    require(q["schema"] == "qar-preexecution-v1" and q["verdict"] == "PASS", "strict independent PASS missing")
    require("reviewer" in q and "reviewer" in q["reviewer"] and "scriber" not in q["reviewer"]
            and "simulator" not in q["reviewer"] and "planner" not in q["reviewer"], "reviewer identity invalid")
    require(q["design_digest"] == actual and q["design_files"] == files, "reviewed design mismatch")
    require(q["preregistration_sha256"] == files["preregistration.json"], "reviewed prereg mismatch")
    require(q["review_markdown_sha256"] == digest(HERE / "PREEXECUTION_REVIEW.md"), "review MD mismatch")
    require(q["label_commit"] == f["label_commit"] and q["label_diff_verdict"] == "PASS", "label review gate")
    stamp = datetime.fromisoformat(q["reviewed_utc"])
    require(stamp.tzinfo is not None and stamp <= datetime.now(timezone.utc), "future/naive review timestamp")
    cfg = read(config)
    require(cfg["schema"] == "qar-phase-prereg-v1" and cfg["replicates"] == 50,
            "unexpected frozen schedule")
    require(cfg["budgets"] == [0, 2, 5] and len(cfg["cells"]) == 4,
            "unexpected frozen grid")
    require(cfg["dgp"]["public_pin"] == PIN, "public source pin mismatch")
    return cfg, q, f


def source_paths():
    return [HERE / p for p in ("measurement_experiment.py", "generation.R", "test_experiment.py",
                               "phase_models.py", "test_phase_models.py")]


def source_manifest():
    paths = source_paths()
    require(all(p.exists() for p in paths), "merged builder source or simulator test absent")
    return {str(p.relative_to(REPO)): digest(p) for p in paths}


def parent_source_provenance():
    # The parent is a read-only comparator, not an imported implementation.
    parent = REPO / "research/emc-part2b/prior_guided/calibration_models.py"
    require(parent.is_file(), "frozen parent source absent")
    for path in source_paths():
        tree = ast.parse(path.read_text()) if path.suffix == ".py" else None
        if tree is not None:
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    require(all(not alias.name.endswith("calibration_models") for alias in node.names),
                            "parent model import requires source freeze")
                if isinstance(node, ast.ImportFrom):
                    require(not (node.module or "").endswith("calibration_models"),
                            "parent model import requires source freeze")
    return {"path": str(parent.relative_to(REPO)), "sha256": digest(parent),
            "use": "read-only reference; no import or call by study code"}


def complete_source_manifest(cfg):
    return {"files": source_manifest(), "parent_reference": parent_source_provenance(),
            "public_package": pinned_package(cfg)}


def package_tree(root):
    require(root.is_dir(), f"missing package tree {root}")
    files = sorted(p for p in root.rglob("*") if p.is_file() and ".git" not in p.parts)
    require(files, "empty package tree")
    h = hashlib.sha256()
    for p in files:
        h.update(str(p.relative_to(root)).encode() + b" " + digest(p).encode() + b"\n")
    return {"root": str(root), "tree_sha256": h.hexdigest(), "files": len(files),
            "provenance_sha256": digest(root / "PROVENANCE.md") if (root / "PROVENANCE.md").exists() else None}


def pinned_package(cfg):
    lib = Path(cfg["dgp"]["library"])
    src = lib.parent / "src"
    require(lib == Path("/tmp/emc-simcodar-d9/lib"), "unapproved package library")
    upstream = Path("/home/matt/Documents/GitHub/simCoDaR")
    require((upstream / ".git").exists(), "public source Git object store absent")
    require(subprocess.check_output(["git", "-C", str(upstream), "cat-file", "-t", PIN],
                                    text=True).strip() == "commit", "public pin absent")
    listing = subprocess.check_output(["git", "-C", str(upstream), "ls-tree", "-r", PIN],
                                      text=True).splitlines()
    listed = {}
    for line in listing:
        meta, name = line.split("\t", 1)
        listed[name] = meta.split()[2]
    actual = {str(p.relative_to(src)): subprocess.check_output(
        ["git", "-C", str(upstream), "hash-object", str(p)], text=True).strip()
        for p in src.rglob("*") if p.is_file()}
    require(actual == listed, "isolated source export diverges from public pinned commit")
    installed = lib / "simCoDaR"
    a = package_tree(src)
    b = package_tree(installed)
    require(digest(src / "inst" / "PROVENANCE.md") == digest(installed / "PROVENANCE.md"),
            "installed/source provenance divergence")
    return {"pin": PIN, "source": a, "installed": b, "git_blobs": actual,
            "library": str(lib), "installed_description_sha256": digest(installed / "DESCRIPTION")}


def historical_seeds():
    # Fixed accepted PR4/PR5 evidence paths. Never scan current study fixtures/config.
    repo_files = (
        "research/emc-part2a/fixtures/source_manifest.json",
        "research/emc-part2a/fixtures/manifest.json",
        "research/emc-part2a/freeze_v1/fixture_manifest.json",
        "research/emc-part2a/freeze_v1/baseline_freeze.json",
        "research/emc-part2b/results/manifest.json",
        "research/emc-part2b/prior_guided/finite_calibration/preregistration.json",
        "research/emc-part2b/prior_guided/finite_calibration/manifests/seeds.json",
        "research/emc-part2b/prior_guided/finite_calibration/attempts/correction-001/manifests/seeds.json",
    )
    run_root = Path("/tmp/emc-part2-workspace/writing/runs")
    run_files = (
        "REQ-20261026-emc-part2/sim-spec.md",
        "REQ-20261027-finite-calibration/sim-spec.md",
        "REQ-20261027-finite-calibration/preregistration.json",
    )
    records, values = [], set()
    def record(path, scope, found):
        require(path.is_file(), f"accepted historical seed evidence absent: {path}")
        found = set(found)
        values.update(found)
        records.append({"path": str(path), "sha256": digest(path), "scope": scope,
                        "explicit_seed_values": sorted(found)})

    for rel in repo_files:
        p = REPO / rel
        require(p.is_file(), f"accepted historical seed evidence absent: {p}")
        obj = read(p)
        found = set()
        def walk(v, key=""):
            if isinstance(v, dict):
                for k, x in v.items():
                    walk(x, k)
            elif isinstance(v, list):
                for x in v:
                    walk(x, key)
            elif isinstance(v, int) and not isinstance(v, bool) and "seed" in key.lower():
                found.add(v)
        walk(obj)
        if rel.endswith("results/manifest.json"):
            found.update(obj["search_seeds"])
        if rel.endswith("finite_calibration/preregistration.json"):
            master = obj["rng"]["master_seed"]
            require(master == 20261027 and obj["replicates"] == 50 and len(obj["cells"]) == 10,
                    "accepted pilot seed schedule changed")
            found.add(master)
            found.update(master + 100000*s + 100*r + j
                         for s in range(1, 11) for r in range(1, 51) for j in range(1, 5))
        if rel.endswith("seeds.json"):
            found.update(int(x["seed"]) for x in obj["streams"])
        record(p, "accepted repository PR4/PR5 evidence: " + rel, found)
    for rel in run_files:
        p = run_root / rel
        source = p.read_text()
        found = set()
        if rel.startswith("REQ-20261026"):
            require("20261026+100*r" in source and "20261027+100*r" in source and
                    "101,102,103" in source, "PR4 run schedule changed")
            found.update((20261026, 20261027, 101, 102, 103))
            found.update(base + 100*r for base in (20261026, 20261027) for r in range(3))
        else:
            require("20261027+100000*s+100*r+k" in source or p.suffix == ".json",
                    "PR5 run schedule changed")
            found.add(20261027)
            found.update(20261027 + 100000*s + 100*r + j
                         for s in range(1, 11) for r in range(1, 51) for j in range(1, 5))
        record(p, "accepted prior workflow run schedule: " + rel, found)
    return records, values


def seed_manifest():
    old, old_values = historical_seeds()
    require(old, "historical manifest inventory empty")
    streams = [{"cell_id": s, "cohort_id": r, "role": role,
                "seed": 731900000 + 100000*s + 100*r + j}
               for s in range(1, 5) for r in range(1, 51)
               for j, role in enumerate(ROLES, 1)]
    fresh = [x["seed"] for x in streams]
    require(len(fresh) == len(set(fresh)) == 2000, "new seed schedule collision")
    collision = sorted(set(fresh) & old_values)
    require(not collision, f"historical seed collision: {collision[:10]}")
    return {"schema": "qar-phase-seeds-v1", "historical_manifests": old,
            "historical_seed_set": sorted(old_values), "streams": streams,
            "collision_count": 0, "historical_inventory_count": len(old)}


def environment(cfg):
    # Nonstochastic R probe is frozen before a generation process exists.
    script = (".libPaths(c('" + cfg["dgp"]["library"] + "',.libPaths()));"
              "library(simCoDaR);library(jsonlite);"
              "cat(jsonlite::toJSON(list(R=R.version.string,RNGkind=RNGkind(),"
              "packages=as.list(sapply(loadedNamespaces(),function(x) as.character(packageVersion(x)))),"
              "libPaths=.libPaths(),BLAS=extSoftVersion()[['BLAS']],LAPACK=La_library()),auto_unbox=TRUE))")
    r = json.loads(subprocess.check_output(["Rscript", "-e", script], text=True))
    versions = {}
    imports = {"numpy", "scipy", "matplotlib", "pytest"}
    for path in (p for p in source_paths() if p.suffix == ".py"):
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Import):
                imports.update(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imports.add(node.module.split(".")[0])
    local = {p.stem for p in source_paths()}
    for name in sorted(imports - local):
        if name in sys.stdlib_module_names or name == "__future__":
            continue
        distributions = importlib.metadata.packages_distributions().get(name, [name])
        versions[name] = {dist: importlib.metadata.version(dist) for dist in distributions}
    require(os.environ.get("OPENBLAS_NUM_THREADS") == "1" and
            os.environ.get("OMP_NUM_THREADS") == "1" and
            os.environ.get("MKL_NUM_THREADS") == "1", "one BLAS thread must be exported")
    return {"python": sys.version, "executable": sys.executable, "packages": versions,
            "R": r, "platform": platform.platform(), "machine": platform.machine(),
            "thread_environment": {k: os.environ.get(k) for k in
                ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "R_PARALLEL_NUM_THREADS")},
            "numpy_config": np.show_config(mode="dicts")}


def frozen(output, cfg, q, f):
    manifests = output / "manifests"
    freeze = read(manifests / "freeze.json")
    require(freeze["design_digest"] == f["design_digest"], "freeze design mismatch")
    require(freeze["review_sha256"] == digest(HERE / "PREEXECUTION_REVIEW.json"), "freeze review mismatch")
    require(freeze["source_sha256"] == digest(manifests / "source.json"), "freeze source manifest mismatch")
    require(freeze["environment_sha256"] == digest(manifests / "environment.json"), "freeze environment manifest mismatch")
    require(freeze["seeds_sha256"] == digest(manifests / "seeds.json"), "freeze seeds manifest mismatch")
    require(read(manifests / "source.json") == complete_source_manifest(cfg),
            "source changed after freeze")
    require(read(manifests / "environment.json") == environment(cfg), "environment changed after freeze")
    require(read(manifests / "seeds.json") == seed_manifest(), "seed inventory changed after freeze")
    return freeze


PRE_CASES = {f"D{i:02d}" for i in (*range(1, 13), 16)}
POST_CASES = {"D13", "D14", "D15"}


def certification_gate(output, *, final=False):
    path = output / "audit" / "checks.json"
    require(path.is_file(), "tester deterministic certification missing")
    cases = read(path).get("cases", [])
    seen = set()
    for case in cases:
        check_id = case.get("check_id")
        if check_id not in PRE_CASES | POST_CASES:
            continue
        require(check_id not in seen, f"duplicate certification {check_id}")
        seen.add(check_id)
        if check_id in POST_CASES and not final:
            require(case.get("verdict") in ("PENDING", "NOT_EXECUTED") and
                    case.get("actual") is None,
                    f"post-run certification executed before run {check_id}")
            continue
        require(case.get("verdict") == "PASS" and case.get("actual") is not None and
                case.get("expected") is not None and case.get("command") and
                case.get("evidence_path"), f"unexecuted/failed certification {check_id}")
        evidence = (output / case["evidence_path"]).resolve()
        require(evidence.is_relative_to(output.resolve()) and evidence.is_file(),
                f"certification evidence unavailable {check_id}")
    required = PRE_CASES | POST_CASES if final else PRE_CASES
    require(required <= seen, f"tester certification incomplete: {sorted(required-seen)}")
    require(POST_CASES <= seen, f"post-run pending rows missing: {sorted(POST_CASES-seen)}")
    return digest(path)


def planned_keys(cfg):
    return {"cohorts": 200, "truth_calls": 800, "specimens": 16000,
            "predictions": 350250, "cohort_metrics": 21250,
            "analysis_blocks": 250, "stage1_predictions": 35000,
            "stage2_predictions": 315250,
            "cells": [x["cell_id"] for x in cfg["cells"]],
            "budgets": cfg["budgets"], "stage1_arms": list(ARMS1), "stage2_arms": list(ARMS2)}


def check(cfg, q, f):
    require(planned_keys(cfg)["predictions"] == 350250, "planned key count changed")
    require(np.linalg.matrix_rank(C) == 2 and np.linalg.matrix_rank(np.vstack([C, H])) == 3,
            "chemical rank invariant")
    require(all(np.linalg.eigvalsh(x).min() > 0 for x in (RB, RQ, R4)), "covariance not SPD")
    seed_manifest()
    source_manifest()  # expected to block before builder merge; never copies builder code
    parent_source_provenance()
    return {"schema": "qar-phase-check-v1", "checked_utc": utc(),
            "design_digest": f["design_digest"], "review_sha256": digest(HERE / "PREEXECUTION_REVIEW.json"),
            "planned": planned_keys(cfg), "verdict": "PASS"}


def freeze(output, cfg, q, f, correction_reason=None, base=None):
    require(not (output / "data").exists() and not (output / "results").exists(),
            "cannot freeze after data/results")
    manifests = output / "manifests"
    require(not (manifests / "freeze.json").exists(), "freeze already exists; preserve record")
    check_result = check(cfg, q, f)
    source = complete_source_manifest(cfg)
    manifests.mkdir(parents=True, exist_ok=True)
    dump(manifests / "design_freeze.json", f)
    dump(manifests / "source.json", source)
    dump(manifests / "environment.json", environment(cfg))
    dump(manifests / "seeds.json", seed_manifest())
    if base is not None:
        require(correction_reason and correction_reason.strip(), "correction reason required")
        old = read(base / "manifests" / "freeze.json")
        require(old["design_digest"] == f["design_digest"], "correction changed frozen design")
        require(read(base / "manifests" / "seeds.json") == read(manifests / "seeds.json"),
                "correction changed seed inventory")
        dump(manifests / "correction.json", {"schema":"qar-phase-correction-v1",
             "reason": correction_reason.strip(), "base_execution_sha256": digest(base / "manifests" / "execution.json"),
             "base_freeze_sha256": digest(base / "manifests" / "freeze.json"),
             "design_digest": f["design_digest"], "seeds_sha256": digest(manifests / "seeds.json"),
             "recorded_utc": utc(), "automatic_retry": False})
    dump(output / "audit" / "checks.json", {"cases": [{"check_id": "simulator-preflight",
         "command": "--check", "expected": "PASS", "actual": "PASS", "tolerance": None,
         "verdict": "PASS", "evidence_path": "manifests/freeze.json"}] +
         [{"check_id": check_id, "command": "after authorized run", "expected": "PASS",
           "actual": None, "tolerance": None, "verdict": "PENDING", "evidence_path": None}
          for check_id in sorted(POST_CASES)], "planned": planned_keys(cfg)})
    dump(manifests / "freeze.json", {"schema": "qar-phase-freeze-v1", "frozen_utc": utc(),
        "design_digest": f["design_digest"], "review_sha256": digest(HERE / "PREEXECUTION_REVIEW.json"),
        "source_sha256": digest(manifests / "source.json"),
        "environment_sha256": digest(manifests / "environment.json"),
        "seeds_sha256": digest(manifests / "seeds.json"), "check": check_result})


def writer(path, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    f = path.open("w", newline="")
    w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
    w.writeheader()
    return f, w


def csv_rows(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def jsonable(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, dict):
        return {k: jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [jsonable(v) for v in obj]
    return obj


def jsonline(f, value):
    f.write(json.dumps(jsonable(value), allow_nan=False, sort_keys=True) + "\n")


def array(row, names):
    return np.array([float(row[x]) for x in names], dtype=float)


def support(group):
    return (True, False, True) if group == "hard" else (True, True, True)


def design(mode, observed, measured):
    if mode == "B" or not measured:
        return C, array(observed, ("SiO2", "TiO2")), RB
    if mode == "R":
        return np.vstack([C, C]), array(observed, ("SiO2", "TiO2", "rep_SiO2", "rep_TiO2")), np.block([[RB, np.zeros((2,2))], [np.zeros((2,2)), RB]])
    if mode == "P":
        return np.vstack([C, np.eye(3)]), array(observed, ("SiO2", "TiO2", "QXRD_Q", "QXRD_A", "QXRD_R")), np.block([[RB, np.zeros((2,3))], [np.zeros((3,2)), RQ]])
    raise ValueError(mode)


def group_rows(rows, cell, cohort):
    return [x for x in rows if int(x["cell_id"]) == cell and int(x["cohort_id"]) == cohort]


def lookup(rows):
    return {x["sample_id"]: x for x in rows}

PRED_FIELDS = ("cell_id", "cohort_id", "support_group", "stage", "budget", "stratum", "arm",
               "sample_id", "status", "reason", "Q_hat", "A_hat", "R_hat", "bulk_SiO2_hat",
               "bulk_TiO2_hat", "bulk_resid_SiO2", "bulk_resid_TiO2", "assay_loss", "prior_loss",
               "objective", "constraint_violation", "kkt_residual", "measurement_design_sse",
               "measurement_design_n", "prior_id", "observation_design_id")
BOUND_FIELDS = ("cell_id", "cohort_id", "support_group", "stage", "budget", "stratum", "sample_id",
                "observation_design_id", "status", "reason", "lo_Q", "lo_A", "lo_R", "hi_Q", "hi_A",
                "hi_R", "width_Q", "width_A", "width_R")
METRIC_FIELDS = ("cell_id", "cohort_id", "support_group", "stage", "budget", "stratum", "arm",
                 "status", "reason", "n_planned", "n_success", "n_failure") + METRICS + (
                 "measurement_design_RMSE_observed", "n_envelope_nonempty", "n_envelope_infeasible",
                 "n_envelope_numerical_failure", "mean_width_Q", "mean_width_A", "mean_width_R")


def prediction_row(keys, est, prior_id, design_id, original_bulk):
    out = {**keys, "status": est.get("status", "numerical_failure"),
           "reason": est.get("reason", ""), "prior_id": prior_id,
           "observation_design_id": design_id}
    mass = est.get("mass")
    if out["status"] != "ok" or mass is None:
        return out
    mass = np.asarray(mass, dtype=float)
    if mass.shape != (3,) or not np.all(np.isfinite(mass)):
        out["status"], out["reason"] = "numerical_failure", "invalid estimator mass"
        return out
    bulk = C @ mass
    out.update(dict(zip(("Q_hat", "A_hat", "R_hat"), mass)))
    out.update({"bulk_SiO2_hat": bulk[0], "bulk_TiO2_hat": bulk[1],
                "bulk_resid_SiO2": bulk[0]-original_bulk[0],
                "bulk_resid_TiO2": bulk[1]-original_bulk[1]})
    for k in ("assay_loss", "prior_loss", "objective", "constraint_violation", "kkt_residual"):
        out[k] = est.get(k)
    residual = est.get("residual")
    if residual is not None:
        residual = np.asarray(residual, dtype=float)
        out["measurement_design_sse"] = float(residual @ residual)
        out["measurement_design_n"] = int(residual.size)
    return out


def run_predictions(output, cfg, deadline):
    # This function has no truth-file access. All fitting uses observation tables only.
    import phase_models as model
    obs = {kind: csv_rows(output / "data" / f"{kind}_observed.csv")
           for kind in ("calibration", "direct", "pool", "evaluation")}
    require([len(obs[x]) for x in obs] == [4000, 4000, 2000, 6000], "observed row counts")
    pf, pw = writer(output / "results" / "predictions.csv", PRED_FIELDS)
    bf, bw = writer(output / "results" / "bounds.csv", BOUND_FIELDS)
    prior_file = (output / "results" / "priors.jsonl").open("w")
    info_file = (output / "results" / "information.jsonl").open("w")
    af, aw = writer(output / "results" / "acquisition.csv",
                    ("cell_id", "cohort_id", "budget", "purchased_ids", "n_calibration_QXRD",
                     "n_pool_base_bulk", "n_replicate_bulk", "n_pool_QXRD", "coverage_pool"))
    counts = {"predictions": 0, "prediction_failures": 0, "bounds": 0,
              "infeasible_envelopes": 0, "numerical_bound_failures": 0,
              "prior_failures": 0, "negative_observations": 0}
    local_cache = {}
    emitted_information = set()
    try:
        for cell in cfg["cells"]:
            s = cell["cell_id"]
            for r in range(1, 51):
                require(time.monotonic() < deadline, "90 minute wall cap reached")
                cohort = {k: group_rows(v, s, r) for k, v in obs.items()}
                cal = np.stack([array(x, ("QXRD_Q", "QXRD_A", "QXRD_R")) for x in cohort["calibration"]])
                pool_ids = [x["sample_id"] for x in cohort["pool"]]
                for k in (0, 2, 5):
                    aw.writerow({"cell_id": s, "cohort_id": r, "budget": k,
                                 "purchased_ids": "|".join(pool_ids[:k]),
                                 "n_calibration_QXRD": 20, "n_pool_base_bulk": k,
                                 "n_replicate_bulk": k, "n_pool_QXRD": k,
                                 "coverage_pool": k/10})
                for group in cell["support_groups"]:
                    sup = support(group)
                    base = {"cell_id": s, "cohort_id": r, "support_group": group}
                    initial_id = f"ps1_s{s}_r{r}_{group}_initial"
                    initial = model.fit_initial(cal, RQ, center_support=tuple(cell["center_support"]))
                    jsonline(prior_file, {**base, "prior_id": initial_id, "update": "initial",
                                          "budget": 0, "excluded_unit": None, "prior": initial})
                    if initial["status"] != "ok":
                        counts["prior_failures"] += 1
                    priors = {("N", 0, None): (initial_id, initial)}
                    for k in (0, 2, 5):
                        for mode in ("B", "R", "P"):
                            units = []
                            for x in cohort["pool"][:k]:
                                D, y, cov = design(mode, x, True)
                                units.append({"unit_id": x["sample_id"], "D": D,
                                              "y": y, "cov": cov})
                            for excluded in [None] + [x["unit_id"] for x in units]:
                                pid = f"ps1_s{s}_r{r}_{group}_{mode}_k{k}_{excluded or 'all'}"
                                prior = model.update_mean(initial, units, support=sup, exclude_unit=excluded)
                                priors[(mode, k, excluded)] = (pid, prior)
                                jsonline(prior_file, {**base, "prior_id": pid, "update": mode,
                                                      "budget": k, "excluded_unit": excluded,
                                                      "prior": prior})
                                if prior["status"] != "ok":
                                    counts["prior_failures"] += 1
                    def emit_bounds(stage, budget, stratum, sample, kind, D, y, cov):
                        did = f"{kind}_{sample['sample_id']}_{group}"
                        key = (did, budget if stage == 2 else -1)
                        if key in local_cache:
                            return did
                        local_cache[key] = True
                        if did not in emitted_information:
                            inf = model.information(D, cov, support=sup)
                            jsonline(info_file, {**base, "observation_design_id": did,
                                                 "support": sup, "information": inf})
                            emitted_information.add(did)
                        b = model.abundance_bounds(y, D, cov, support=sup)
                        row = {**base, "stage": stage, "budget": budget,
                               "stratum": stratum, "sample_id": sample["sample_id"],
                               "observation_design_id": did, "status": b["status"],
                               "reason": b.get("reason", "")}
                        if b["status"] == "ok":
                            for prefix, value in (("lo", b["lower"]), ("hi", b["upper"]),
                                                  ("width", b["width"])):
                                row.update({f"{prefix}_{p}": v for p, v in zip("QAR", value)})
                        elif b["status"] == "infeasible_envelope":
                            counts["infeasible_envelopes"] += 1
                        else:
                            counts["numerical_bound_failures"] += 1
                        bw.writerow(row)
                        counts["bounds"] += 1
                        return did
                    def emit_prediction(stage, budget, stratum, arm, sample, D, y, cov, pid, prior, did):
                        method = arm[:2]
                        if method == "A0":
                            prior, pid = None, None
                        key = (method, sample["sample_id"], did, pid)
                        if key in fit_cache:
                            est = fit_cache[key]
                        else:
                            try:
                                est = model.estimate(None, None, None, prior, method="A1", support=sup) if method == "A1" else model.estimate(
                                    y, D, cov, prior, method=method, strength=cell["strength"], support=sup)
                            except (ValueError, np.linalg.LinAlgError, ArithmeticError) as ex:
                                est = {"status": "numerical_failure", "reason": str(ex), "mass": None}
                            fit_cache[key] = est
                        keys = {**base, "stage": stage, "budget": budget, "stratum": stratum,
                                "arm": arm, "sample_id": sample["sample_id"]}
                        bulk = array(sample, ("SiO2", "TiO2"))
                        row = prediction_row(keys, est, pid, did, bulk)
                        pw.writerow(row)
                        counts["predictions"] += 1
                        if row["status"] != "ok":
                            counts["prediction_failures"] += 1
                    fit_cache = {}
                    for x in cohort["direct"]:
                        for kind, indices in (("B", [0,1]), ("T", [0,1,2]), ("H", [0,1,3])):
                            D, y, cov = D4[indices], array(x, ("SiO2","TiO2","redundant_Ti","H_A"))[indices], R4[np.ix_(indices,indices)]
                            did = emit_bounds(1, -1, "direct", x, kind, D, y, cov)
                            for arm in ARMS1:
                                if arm == "A1" and kind == "B" or arm.endswith("_"+kind):
                                    emit_prediction(1, -1, "direct", arm, x, D, y, cov,
                                                    initial_id, initial, did)
                    for k in (0, 2, 5):
                        for stratum, rows in (("unmeasured", cohort["evaluation"]),
                                              ("measured", cohort["pool"][:k])):
                            for x in rows:
                                local = {}
                                for mode in (("B",) if stratum == "unmeasured" else ("B","R","P")):
                                    D, y, cov = design(mode, x, stratum == "measured")
                                    did = emit_bounds(2, k, stratum, x, mode, D, y, cov)
                                    local[mode] = (D, y, cov, did)
                                for arm in ARMS2:
                                    method = arm[:2]
                                    tag = arm.split("_", 1)[1]
                                    update = "N" if tag == "N" or tag.endswith("0") else tag
                                    local_mode = "B" if stratum == "unmeasured" or tag == "N" else tag[0]
                                    D, y, cov, did = local[local_mode]
                                    if update == "N":
                                        pid, prior = initial_id, initial
                                    else:
                                        exclude = x["sample_id"] if stratum == "measured" else None
                                        pid, prior = priors[(update, k, exclude)]
                                    emit_prediction(2, k, stratum, arm, x, D, y, cov, pid, prior, did)
                    fit_cache.clear()
                    local_cache.clear()
                pf.flush(); bf.flush(); prior_file.flush(); info_file.flush()
                if (output / "results" / "predictions.csv").stat().st_size > 500_000_000:
                    raise RuntimeError("500 MB output cap reached")
    finally:
        for file in (pf, bf, prior_file, info_file, af):
            file.close()
    require(counts["predictions"] == 350250, f"prediction count {counts['predictions']}")
    return counts


def score_block(rows, truths, observed, bounds):
    n = len(rows)
    good = [x for x in rows if x["status"] == "ok"]
    result = {"n_planned": n, "n_success": len(good), "n_failure": n-len(good)}
    if len(good) != n or not n:
        result["status"] = "empty_stratum" if n == 0 else "unavailable_matched_failure"
        result["reason"] = "no purchased specimens" if n == 0 else "one or more required predictions failed"
        return result
    hat = np.stack([array(x, ("Q_hat","A_hat","R_hat")) for x in good])
    truth = np.stack([array(truths[x["sample_id"]], ("Q","A","R")) for x in good])
    bulk_obs = np.stack([array(observed[x["sample_id"]], ("SiO2","TiO2")) for x in good])
    e = hat-truth
    result.update({"status": "ok", "reason": "",
                   "mass_RMSE": float(np.sqrt(np.mean(e*e))),
                   "totalTi_RMSE": float(np.sqrt(np.mean((e[:,1]+e[:,2])**2))),
                   "bulk_assay_RMSE_observed": float(np.sqrt(np.mean((hat@C.T-bulk_obs)**2))),
                   "bulk_assay_RMSE_latent": float(np.sqrt(np.mean((hat@C.T-truth@C.T)**2)))})
    for j, phase in enumerate("QAR"):
        result[f"phase_RMSE_{phase}"] = float(np.sqrt(np.mean(e[:,j]**2)))
        result[f"phase_bias_{phase}"] = float(np.mean(e[:,j]))
    # Own-design residual is incomparable across different observation designs.
    residual_rows = [x for x in good if x.get("measurement_design_n") not in (None, "")]
    result["measurement_design_RMSE_observed"] = (
        float(np.sqrt(sum(float(x["measurement_design_sse"]) for x in residual_rows) /
                      sum(int(x["measurement_design_n"]) for x in residual_rows)))
        if residual_rows else None)
    chosen = {x["observation_design_id"] for x in good}
    bounds = [b for b in bounds if b["observation_design_id"] in chosen]
    widths = [b for b in bounds if b["status"] == "ok"]
    result["n_envelope_nonempty"] = len(widths)
    result["n_envelope_infeasible"] = sum(b["status"] == "infeasible_envelope" for b in bounds)
    result["n_envelope_numerical_failure"] = sum(b["status"] not in ("ok", "infeasible_envelope") for b in bounds)
    for phase in "QAR":
        result[f"mean_width_{phase}"] = float(np.mean([float(b[f"width_{phase}"]) for b in widths])) if widths else None
    return result


def score(output, deadline):
    # Truth files are opened only after the complete prediction CSV is closed and hashed.
    truth = {}
    observed = {}
    for kind in ("direct", "pool", "evaluation"):
        truth.update(lookup(csv_rows(output / "data" / f"{kind}_truth.csv")))
        observed.update(lookup(csv_rows(output / "data" / f"{kind}_observed.csv")))
    bounds_by = {}
    for b in csv_rows(output / "results" / "bounds.csv"):
        k = (int(b["cell_id"]), int(b["cohort_id"]), b["support_group"],
             int(b["stage"]), int(b["budget"]), b["stratum"])
        bounds_by.setdefault(k, []).append(b)
    groups = {}
    for p in csv_rows(output / "results" / "predictions.csv"):
        k = (int(p["cell_id"]), int(p["cohort_id"]), p["support_group"],
             int(p["stage"]), int(p["budget"]), p["stratum"], p["arm"])
        groups.setdefault(k, []).append(p)
    mf, mw = writer(output / "results" / "cohort_metrics.csv", METRIC_FIELDS)
    records = {}
    blocked = 0
    try:
        for s in range(1,5):
            for r in range(1,51):
                require(time.monotonic() < deadline, "90 minute wall cap reached while scoring")
                supports = ("soft","hard") if s == 4 else ("full",)
                for sup in supports:
                    for stage, budget, stratum, arms, n in (
                        (1,-1,"direct",ARMS1,20),
                        *((2,k,st,ARMS2,(30 if st == "unmeasured" else k))
                          for k in (0,2,5) for st in ("unmeasured","measured"))):
                        block = (s,r,sup,stage,budget,stratum)
                        # A failure in any required arm invalidates the entire matched block.
                        bad = any(len(groups.get(block+(arm,), [])) != n or
                                  any(x["status"] != "ok" for x in groups.get(block+(arm,), []))
                                  for arm in arms) if n else False
                        if bad:
                            blocked += 1
                        b = bounds_by.get(block, [])
                        for arm in arms:
                            key = block+(arm,)
                            rows = groups.get(key, [])
                            result = score_block(rows, truth, observed, b)
                            result["n_planned"] = n
                            result["n_failure"] = n-result["n_success"]
                            if bad:
                                for metric in METRICS:
                                    result.pop(metric, None)
                                result["status"] = "unavailable_matched_failure"
                                result["reason"] = "matched arm failed or missing"
                            record = {"cell_id":s,"cohort_id":r,"support_group":sup,
                                      "stage":stage,"budget":budget,"stratum":stratum,"arm":arm,
                                      **result}
                            mw.writerow(record)
                            records[key] = record
    finally:
        mf.close()
    require(len(records) == 21250, "cohort metric row count")
    sf, sw = writer(output / "results" / "paired_summary.csv",
                    ("cell_id","support_group","stage","budget","stratum","contrast",
                     "metric","n_planned","n_complete","n_failed","mean_difference","descriptive_se"))
    try:
        for s in range(1,5):
            for sup in (("soft","hard") if s == 4 else ("full",)):
                for stage, budget, stratum, contrasts in (
                    (1,-1,"direct",CONTRASTS1),
                    *((2,k,st,CONTRASTS2) for k in (0,2,5)
                      for st in ("unmeasured","measured"))):
                    for contrast in contrasts:
                        a,b = contrast.split("-")
                        for metric in METRICS:
                            diffs = []
                            for r in range(1,51):
                                x = records[(s,r,sup,stage,budget,stratum,a)]
                                y = records[(s,r,sup,stage,budget,stratum,b)]
                                if x["status"] == y["status"] == "ok":
                                    diffs.append(float(x[metric])-float(y[metric]))
                            sw.writerow({"cell_id":s,"support_group":sup,"stage":stage,
                                "budget":budget,"stratum":stratum,"contrast":contrast,"metric":metric,
                                "n_planned":50,"n_complete":len(diffs),"n_failed":50-len(diffs),
                                "mean_difference": float(np.mean(diffs)) if len(diffs)==50 else None,
                                "descriptive_se": float(np.std(diffs,ddof=1)/np.sqrt(50)) if len(diffs)==50 else None})
    finally:
        sf.close()
    return {"cohort_metrics": len(records), "matched_failed_blocks": blocked,
            "paired_rows": 5*(len(CONTRASTS1)+6*len(CONTRASTS2))*len(METRICS)}


def figures(output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = csv_rows(output / "results" / "paired_summary.csv")
    selections = (
        ("direct.png", 1, -1, "direct", "A0_H-A0_T", "Direct: independent phase row minus redundant Ti"),
        ("transfer.png", 2, 5, "unmeasured", "A2_P-A2_R", "Independent transfer: phase versus replicated bulk"),
        ("occurrence.png", 2, 5, "unmeasured", "A2_P-A2_N", "Phase update versus no update"),
    )
    for filename, stage, budget, stratum, contrast, title in selections:
        chosen = [x for x in rows if int(x["stage"]) == stage and
                  int(x["budget"]) == budget and x["stratum"] == stratum and
                  x["contrast"] == contrast and x["metric"] == "mass_RMSE"]
        chosen.sort(key=lambda x: (int(x["cell_id"]), x["support_group"]))
        labels = [f"{x['cell_id']}:{x['support_group']}" for x in chosen]
        values = [float(x["mean_difference"]) if x["mean_difference"] else np.nan for x in chosen]
        errors = [float(x["descriptive_se"]) if x["descriptive_se"] else 0. for x in chosen]
        fig, ax = plt.subplots(figsize=(6.5, 3.7))
        ax.errorbar(labels, values, yerr=errors, fmt="o", capsize=3)
        ax.axhline(0, color="black", linewidth=.8)
        ax.set(title=title, ylabel="Paired cohort mass RMSE difference (descriptive SE)",
               xlabel="Cell : support group")
        fig.tight_layout()
        path = output / "figures" / filename
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=150)
        plt.close(fig)


def study_diagnostics(output, cfg, counts, prediction_hash):
    centers = {}
    eigen_t, eigen_m, diag_t, diag_m, correlations = [], [], [], [], []
    floor_count = 0
    priors_count = 0
    with (output / "results" / "priors.jsonl").open() as f:
        for line in f:
            x = json.loads(line)
            p = x["prior"]
            if p["status"] != "ok":
                continue
            priors_count += 1
            center = np.asarray(p["center"], float)
            if x["update"] == "initial":
                centers[(x["cell_id"],x["cohort_id"],x["support_group"])] = center
                d = p["diagnostics"]
                if "deconvolved_eigenvalues" in d and "floored_eigenvalues" in d:
                    floor_count += sum(float(a) < float(b) for a,b in zip(
                        d["deconvolved_eigenvalues"], d["floored_eigenvalues"]))
            T = np.asarray(p["tangent_cov"],float)
            M = np.asarray(p["mean_cov"],float)
            if T.shape == M.shape == (2,2):
                eigen_t.append(np.linalg.eigvalsh(T).tolist())
                eigen_m.append(np.linalg.eigvalsh(M).tolist())
                U = np.array([[2/np.sqrt(6),0],[-1/np.sqrt(6),1/np.sqrt(2)],
                              [-1/np.sqrt(6),-1/np.sqrt(2)]])
                physical = U@T@U.T
                diag_t.append(np.diag(physical).tolist())
                diag_m.append(np.diag(U@M@U.T).tolist())
                correlations.append(float(physical[1,2]/np.sqrt(physical[1,1]*physical[2,2])))
    cal = csv_rows(output / "data" / "calibration_truth.csv")
    ev = csv_rows(output / "data" / "evaluation_truth.csv")
    retrospective = []
    for (s,r,group), center in centers.items():
        for kind, rows in (("calibration",cal),("evaluation",ev)):
            latent = np.mean([array(x,("Q","A","R")) for x in group_rows(rows,s,r)],axis=0)
            retrospective.append({"cell_id":s,"cohort_id":r,"support_group":group,
                                  "population":kind,"center_minus_latent_mean":(center-latent).tolist()})
    by_group = {}
    for (s,r,g), center in centers.items():
        by_group.setdefault(f"{s}:{g}",[]).append(center)
    sd = {k: np.std(v,axis=0,ddof=1).tolist() for k,v in by_group.items()}
    negative = {}
    for kind in ("calibration","direct","pool","evaluation"):
        rows = csv_rows(output / "data" / f"{kind}_observed.csv")
        numeric = [key for key in rows[0] if key not in ("cell_id","cohort_id","sample_id")]
        negative[kind] = sum(float(x[key]) < 0 for x in rows for key in numeric if x[key] != "")
    return {"schema":"qar-phase-diagnostics-v1", "planned":planned_keys(cfg),
            "attempted":counts,"prediction_sha256_before_truth_access":prediction_hash,
            "negative_observations_by_stream":negative,"initial_prior_count":len(centers),
            "successful_prior_records":priors_count,"calibration_floor_eigenvalue_count":floor_count,
            "learned_center_sd_by_cell_support":sd,
            "mean_T_eigenvalues":np.mean(eigen_t,axis=0).tolist() if eigen_t else None,
            "mean_M_eigenvalues":np.mean(eigen_m,axis=0).tolist() if eigen_m else None,
            "mean_T_physical_diagonal":np.mean(diag_t,axis=0).tolist() if diag_t else None,
            "mean_M_physical_diagonal":np.mean(diag_m,axis=0).tolist() if diag_m else None,
            "mean_learned_A_R_correlation":float(np.mean(correlations)) if correlations else None,
            "retrospective_center_errors":retrospective,
            "check_evidence_links":["audit/checks.json","manifests/freeze.json"]}


def output_hashes(output):
    return {str(p.relative_to(output)): digest(p) for p in sorted(output.rglob("*"))
            if p.is_file() and p.name not in ("execution.json", "generation_log.json")
            and p != output / "audit" / "checks.json"}


def verify(output, cfg, q, f):
    frozen(output, cfg, q, f)
    execution = read(output / "manifests" / "execution.json")
    require(execution["schema"] == "qar-phase-execution-v1", "execution schema")
    require(execution["status"] == "complete", "execution not complete")
    require(execution["review_sha256"] == digest(HERE / "PREEXECUTION_REVIEW.json"), "execution review hash")
    require(execution["freeze_sha256"] == digest(output / "manifests" / "freeze.json"), "execution freeze hash")
    pre_checks = output / "manifests" / "pre_run_checks.json"
    require(execution["certification_sha256"] == digest(pre_checks), "pre-run certification hash")
    current = {c.get("check_id"): c for c in read(output / "audit" / "checks.json")["cases"]}
    for case in read(pre_checks)["cases"]:
        if case.get("check_id") in PRE_CASES:
            require(current.get(case["check_id"]) == case, f"pre-run certification changed {case['check_id']}")
    certification_gate(output, final=True)
    review_time = datetime.fromisoformat(q["reviewed_utc"])
    freeze_time = datetime.fromisoformat(read(output / "manifests" / "freeze.json")["frozen_utc"])
    draw_time = datetime.fromisoformat(execution["first_draw_utc"])
    require(review_time <= freeze_time < draw_time, "review/freeze/draw order")
    current_hashes = output_hashes(output)
    require(all(current_hashes.get(path) == sha for path, sha in execution["output_hashes"].items()),
            "output hash mismatch")
    require(all(path in execution["output_hashes"] or path.startswith("audit/")
                for path in current_hashes), "unexpected post-run output")
    return {"verdict": "PASS", "verified_utc": utc(), "counts": execution["counts"]}


def execute(output, cfg, q, f):
    freeze_record = frozen(output, cfg, q, f)
    certification_sha = certification_gate(output)
    if "attempts" in output.parts:
        correction = read(output / "manifests" / "correction.json")
        require(correction["schema"] == "qar-phase-correction-v1" and
                correction["design_digest"] == f["design_digest"] and
                correction["seeds_sha256"] == digest(output / "manifests" / "seeds.json"),
                "correction lineage invalid")
    require(not (output / "data").exists() and not (output / "results").exists(),
            "existing data/results; use fresh correction attempt")
    require(not (output / "manifests" / "execution.json").exists(), "existing execution ledger")
    pre_checks = output / "manifests" / "pre_run_checks.json"
    require(not pre_checks.exists(), "existing pre-run certification snapshot")
    pre_checks.write_bytes((output / "audit" / "checks.json").read_bytes())
    require(digest(pre_checks) == certification_sha, "pre-run certification snapshot mismatch")
    deadline = time.monotonic() + 90*60
    manifests = output / "manifests"
    ledger = {"schema":"qar-phase-execution-v1", "status":"started",
              "review_sha256": digest(HERE / "PREEXECUTION_REVIEW.json"),
              "preregistration_sha256": digest(HERE / "preregistration.json"),
              "design_digest": f["design_digest"],
              "freeze_sha256": digest(manifests / "freeze.json"),
              "certification_sha256": certification_sha,
              "reviewed_utc": q["reviewed_utc"], "frozen_utc": freeze_record["frozen_utc"],
              "first_draw_utc": utc(), "planned": planned_keys(cfg)}
    require(datetime.fromisoformat(ledger["frozen_utc"]) < datetime.fromisoformat(ledger["first_draw_utc"]),
            "freeze must precede first draw")
    dump(manifests / "execution.json", ledger)  # durable start record before R process
    try:
        result = subprocess.run(["Rscript", str(HERE / "generation.R"),
                                 str(HERE / "preregistration.json"), str(output),
                                 cfg["dgp"]["library"]], capture_output=True, text=True,
                                timeout=max(1, int(deadline-time.monotonic())))
        dump(manifests / "generation_log.json", {"returncode": result.returncode,
             "stdout": result.stdout, "stderr": result.stderr, "ended_utc": utc()})
        require(result.returncode == 0, f"generation failed: {result.stderr[-500:]}")
        require(time.monotonic() < deadline, "90 minute wall cap reached after generation")
        data_hashes = {str(p.relative_to(output)): digest(p) for p in sorted((output/"data").glob("*.csv"))}
        require(len(data_hashes) == 8, "eight data tables required")
        counts = run_predictions(output, cfg, deadline)
        prediction_hash = digest(output / "results" / "predictions.csv")
        dump(manifests / "access.json", {"schema": "qar-phase-access-v1",
             "fit_input_paths": [f"data/{k}_observed.csv" for k in
                 ("calibration","direct","pool","evaluation")],
             "truth_opened_only_after_prediction_sha256": prediction_hash,
             "data_sha256": data_hashes, "recorded_utc": utc()})
        counts.update(score(output, deadline))
        require(counts["prediction_failures"] == counts["prior_failures"] ==
                counts["numerical_bound_failures"] == counts["matched_failed_blocks"] == 0,
                "numerical/fit/LP failure retained; run BLOCK")
        figures(output)
        dump(output / "results" / "diagnostics.json",
             study_diagnostics(output, cfg, counts, prediction_hash))
        dump(manifests / "deviations.json", {"schema": "qar-phase-deviations-v1",
            "scientific_design_changes": [], "automatic_resampling": False,
            "status": "none", "recorded_utc": utc()})
        dump(manifests / "preservation.json", {"schema": "qar-phase-preservation-v1",
             "old_study_untouched": True, "accepted_parent": cfg["accepted_parent"]})
        require(sum(p.stat().st_size for p in output.rglob("*") if p.is_file()) <= 500_000_000,
                "500 MB output cap exceeded")
        ledger.update({"status":"complete", "ended_utc":utc(), "counts":counts,
                       "output_hashes":output_hashes(output)})
        dump(manifests / "execution.json", ledger)
    except Exception as exc:
        ledger.update({"status":"blocked", "ended_utc":utc(), "reason":str(exc),
                       "output_hashes":output_hashes(output)})
        dump(manifests / "execution.json", ledger)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--review", type=Path, required=True)
    parser.add_argument("--attempt", help="correction-NNN; preserves base evidence")
    parser.add_argument("--correction-reason", help="required with --attempt --freeze")
    mode = parser.add_mutually_exclusive_group(required=True)
    for name in ("check", "freeze", "run", "verify"):
        mode.add_argument("--" + name, action="store_true")
    args = parser.parse_args()
    try:
        cfg, q, f = design_gate(args.config, args.review)
        base = args.output.resolve()
        require(base == HERE / "run", "output must be declared study run directory")
        output = base
        if args.attempt:
            require(re.fullmatch(r"correction-[0-9]{3}", args.attempt), "attempt name")
            output = base / "attempts" / args.attempt
            require((base / "manifests" / "execution.json").exists(), "base attempt must be preserved")
            require(read(base / "manifests" / "execution.json")["status"] == "blocked",
                    "correction requires blocked base attempt")
        else:
            require(args.correction_reason is None, "correction reason requires an attempt")
        if args.check:
            result = check(cfg, q, f)
        elif args.freeze:
            freeze(output, cfg, q, f, args.correction_reason, base if args.attempt else None)
            result = {"verdict":"FROZEN", "output":str(output)}
        elif args.run:
            execute(output, cfg, q, f)
            result = {"verdict":"RUN_COMPLETE_PENDING_POSTCHECKS", "output":str(output)}
        else:
            result = verify(output, cfg, q, f)
        print(json.dumps(result, sort_keys=True, allow_nan=False))
    except Exception as exc:
        print(f"BLOCK: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
