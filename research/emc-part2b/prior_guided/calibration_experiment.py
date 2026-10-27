#!/usr/bin/env python3
"""Frozen finite-calibration QAR pilot. No random draws occur in --check."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import importlib.util
import json
import math
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import tarfile
import tempfile
import time
from datetime import datetime, timezone

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PLAN = HERE / "FINITE_CALIBRATION_PLAN.md"
MODEL = HERE / "calibration_models.py"
GENERATOR = HERE / "finite_calibration" / "generation.R"
AMENDMENT = HERE / "finite_calibration" / "manifests" / "correction_amendment.json"
PIN_REPO = Path("/home/matt/Documents/GitHub/simCoDaR")
PHASES = ("Quartz", "Anatase", "Rutile")
CHANNELS = ("SiO2", "TiO2")
STREAMS = ("calibration_truth", "evaluation_truth", "calibration_noise", "assay_noise")
CONTRAST_METRICS = ("mass_RMSE", "assay_RMSE_observed", "assay_RMSE_truth")
PIN = "d9a0b27676a0cb3a4acb300a3f6b2ebc845d205d"
FIRST_SOURCE_COMMIT = "ccf2c23"
CORRECTION_ATTEMPT = "correction-001"


class ContractError(RuntimeError):
    """Frozen contract, environment, or output violation."""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def strict_json(path: Path) -> dict:
    def reject(value: str) -> None:
        raise ContractError(f"nonfinite JSON token in {path}: {value}")

    value = json.loads(path.read_text(), parse_constant=reject)
    if not isinstance(value, dict):
        raise ContractError(f"expected JSON object: {path}")
    return value


def dump_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(REPO.resolve()))


def validate_config(config: dict) -> None:
    """Refuse amendments to the ten-cell frozen plan."""
    if config.get("schema") != "qar-finite-calibration-v1":
        raise ContractError("wrong config schema")
    if config.get("accepted_parent") != "8ce5da876c23503819b4fbcdfed2d8a6393cfbe0":
        raise ContractError("wrong accepted parent")
    if config["dependency"]["commit"] != PIN or config["dependency"]["library"] != "/tmp/emc-simcodar-d9/lib":
        raise ContractError("wrong public dependency pin")
    if config["phase_order"] != list(PHASES) or config["channels"] != list(CHANNELS):
        raise ContractError("wrong phase/channel order")
    if config["C"] != [[1, 0, 0], [0, 1, 1]]:
        raise ContractError("wrong chemistry matrix")
    if config["rng"]["engine"] != ["Mersenne-Twister", "Inversion", "Rejection"] or config["rng"]["master_seed"] != 20261027:
        raise ContractError("wrong RNG")
    if config["replicates"] != 50 or config["evaluation_n"] != 40 or config["methods"] != ["A0", "A1", "A2"]:
        raise ContractError("wrong pilot bounds")
    cells = config["cells"]
    expected = [
        (1, 10, .02, 1, .6, (.6, .1, .3), (True, True, True)),
        (2, 30, .02, 1, .6, (.6, .1, .3), (True, True, True)),
        (3, 100, .02, 1, .6, (.6, .1, .3), (True, True, True)),
        (4, 10, .02, 1, None, (.6, .1, .3), (True, True, True)),
        (5, 30, .02, 1, None, (.6, .1, .3), (True, True, True)),
        (6, 100, .02, 1, None, (.6, .1, .3), (True, True, True)),
        (7, 30, .02, .01, None, (.6, .1, .3), (True, True, True)),
        (8, 30, .02, 1, None, (.6, .3, .1), (True, True, True)),
        (9, 30, .08, 1, None, (.6, .1, .3), (True, True, True)),
        (10, 30, .02, 1, None, (.6, .1, .3), (True, False, True)),
    ]
    if len(cells) != 10:
        raise ContractError("expected ten cells")
    for cell, (ident, n, sd, strength, q, center, support) in zip(cells, expected):
        if (cell["cell_id"], cell["calibration_n"], cell["qxrd_sd"],
            cell["strength"], cell["fixed_quartz"], tuple(cell["calibration_center"]),
            tuple(cell["support"])) != (ident, n, sd, strength, q, center, support):
            raise ContractError(f"cell {ident} differs from frozen plan")
    if config["budget"] != {"cells": 10, "cohorts": 500, "fits": 500,
                            "evaluation_specimens": 20000, "arm_predictions": 60000,
                            "public_generator_calls": 1000, "streams": 2000,
                            "max_wall_seconds": 3600, "max_workers": 2,
                            "max_output_bytes": 250000000, "paid_cost": 0}:
        raise ContractError("wrong budget")


def prereg_check(config_path: Path, output: Path) -> dict:
    config = strict_json(config_path)
    validate_config(config)
    if output.resolve() != (HERE / "finite_calibration").resolve() or config_path.resolve() != (output / "preregistration.json").resolve():
        raise ContractError("config and output must be the frozen in-repository paths")
    freeze = strict_json(output / "manifests" / "preregistration_freeze.json")
    if freeze.get("status") != "PREREG_PASS" or freeze.get("generation_performed") is not False:
        raise ContractError("preliminary PREREG PASS absent")
    for file_name, record in freeze["files"].items():
        path = REPO / file_name
        if not path.is_file() or sha256(path) != record["sha256"]:
            raise ContractError(f"preregistration hash changed: {file_name}")
    return config


def seed_rows(config: dict) -> list[dict]:
    rows = []
    seen = set()
    for cell in config["cells"]:
        s = cell["cell_id"]
        for r in range(1, 51):
            for k, name in enumerate(STREAMS, 1):
                seed = 20261027 + 100000 * s + 100 * r + k
                if seed in seen:
                    raise ContractError("duplicate RNG stream seed")
                seen.add(seed)
                rows.append({"cell_id": s, "cohort_id": r, "stream_id": k,
                             "stream": name, "seed": seed})
    if len(rows) != 2000:
        raise ContractError("wrong seed count")
    return rows


def source_fingerprint(config_path: Path) -> dict:
    paths = (PLAN, config_path, Path(__file__), GENERATOR, MODEL,
             HERE / "test_calibration_models.py", HERE / "test_calibration_experiment.py",
             AMENDMENT)
    if any(not p.is_file() for p in paths):
        raise ContractError("merged estimator source is not present yet")
    return {relative(p): sha256(p) for p in paths}


def run_command(args: list[str], **kwargs) -> subprocess.CompletedProcess:
    p = subprocess.run(args, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                       check=False, **kwargs)
    if p.returncode:
        raise ContractError(f"command failed ({p.returncode}): {' '.join(args[:3])}\n{p.stderr[-3000:]}")
    return p


def pin_environment() -> dict:
    if not PIN_REPO.is_dir():
        raise ContractError("pinned public source checkout unavailable")
    commit = run_command(["git", "-C", str(PIN_REPO), "rev-parse", PIN]).stdout.strip()
    if commit != PIN:
        raise ContractError("public source commit mismatch")
    tree = run_command(["git", "-C", str(PIN_REPO), "rev-parse", f"{PIN}^{{tree}}"]).stdout.strip()
    r = run_command(["Rscript", "-e", "cat(as.character(getRversion()),'|',as.character(packageVersion('jsonlite')))"]).stdout.strip()
    if not (Path("/tmp/emc-simcodar-d9/lib/simCoDaR").is_dir() and
            Path("/tmp/emc-simcodar-d9/src").is_dir()):
        raise ContractError("declared isolated package/source absent")
    # Verify the exported source against the git object, rather than package metadata.
    for name in ("DESCRIPTION", "NAMESPACE", "R/scenario.R", "R/mineral_scenario.R",
                 "R/mineral_conversion.R", "R/utils-ilr.R"):
        expected = subprocess.run(["git", "-C", str(PIN_REPO), "show", f"{PIN}:{name}"],
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        source = Path("/tmp/emc-simcodar-d9/src") / name
        if expected.returncode or not source.is_file() or expected.stdout != source.read_bytes():
            raise ContractError(f"isolated source differs from pinned export: {name}")
    return {"python": platform.python_version(), "python_implementation": platform.python_implementation(),
            "r_jsonlite": r, "platform": platform.platform(), "pin": PIN,
            "pin_tree": tree, "RNGkind": ["Mersenne-Twister", "Inversion", "Rejection"],
            "package_version": "0.1.0.9000"}


def isolated_install(temp: Path) -> Path:
    """Reinstall exact git export into a private temporary R library."""
    archive = subprocess.run(["git", "-C", str(PIN_REPO), "archive", "--format=tar", PIN],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if archive.returncode:
        raise ContractError("cannot export pinned public source")
    src = temp / "src"
    lib = temp / "lib"
    src.mkdir()
    lib.mkdir()
    with tarfile.open(fileobj=io.BytesIO(archive.stdout), mode="r:") as tar:
        for member in tar.getmembers():
            dest = (src / member.name).resolve()
            if not dest.is_relative_to(src.resolve()) or member.issym() or member.islnk():
                raise ContractError("unsafe pinned source archive")
        tar.extractall(src, filter="data")
    run_command(["R", "CMD", "INSTALL", "--library=" + str(lib), str(src)], timeout=180)
    check = run_command(["Rscript", "-e", "cat(as.character(packageVersion('simCoDaR', lib.loc=commandArgs(TRUE)[1])))", str(lib)])
    if check.stdout.strip() != "0.1.0.9000":
        raise ContractError("pinned package installation version mismatch")
    return lib


def freeze_check(config_path: Path, output: Path, *, create: bool) -> dict:
    config = prereg_check(config_path, output)
    sources = source_fingerprint(config_path)
    environment = pin_environment()
    manifests = output / "manifests"
    source_path = manifests / "source.json"
    environment_path = manifests / "environment.json"
    freeze_path = manifests / "freeze.json"
    if all(p.exists() for p in (source_path, environment_path, freeze_path)):
        frozen = strict_json(freeze_path)
        if strict_json(source_path).get("files") != sources or strict_json(environment_path).get("environment") != environment:
            raise ContractError("source or environment changed after freeze")
        if frozen.get("status") != "SOURCE_ENVIRONMENT_FROZEN" or frozen.get("source_sha256") != sha256(source_path) or frozen.get("environment_sha256") != sha256(environment_path):
            raise ContractError("source/environment freeze hash mismatch")
    elif create:
        if any(p.exists() for p in (source_path, environment_path, freeze_path)):
            raise ContractError("partial source/environment freeze")
        if (output / "data").exists() or (output / "results").exists():
            raise ContractError("cannot freeze after output generation")
        dump_json(source_path, {"schema": "qar-source-v1", "files": sources})
        dump_json(environment_path, {"schema": "qar-environment-v1", "environment": environment})
        dump_json(freeze_path, {"schema": "qar-source-environment-freeze-v1",
                                "status": "SOURCE_ENVIRONMENT_FROZEN", "timestamp_utc": utc(),
                                "preregistration_sha256": sha256(output / "manifests" / "preregistration_freeze.json"),
                                "source_sha256": sha256(source_path),
                                "environment_sha256": sha256(environment_path)})
    else:
        raise ContractError("source/environment freeze missing; run --check after merge")
    return config


def first_attempt_snapshot(output: Path) -> dict:
    """Verify the failed first source against git and fingerprint its retained files."""
    first_source = strict_json(output / "manifests" / "source.json")
    first_environment = output / "manifests" / "environment.json"
    first_freeze = strict_json(output / "manifests" / "freeze.json")
    first_execution = strict_json(output / "manifests" / "execution.json")
    if first_execution.get("status") != "INCOMPLETE" or "is.finite(x)" not in str(first_execution.get("failure")):
        raise ContractError("expected retained first draw failure is absent")
    if first_freeze.get("source_sha256") != sha256(output / "manifests" / "source.json") or first_freeze.get("environment_sha256") != sha256(first_environment):
        raise ContractError("initial source/environment freeze changed")
    if first_execution.get("source_freeze_sha256") != sha256(output / "manifests" / "freeze.json"):
        raise ContractError("initial execution references a different source freeze")
    for name, digest in first_source["files"].items():
        old = subprocess.run(["git", "-C", str(REPO), "show", f"{FIRST_SOURCE_COMMIT}:{name}"],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        if old.returncode or hashlib.sha256(old.stdout).hexdigest() != digest:
            raise ContractError(f"first source cannot be resolved at {FIRST_SOURCE_COMMIT}: {name}")
    files = {}
    for section in ("data", "results", "figures", "manifests"):
        directory = output / section
        if directory.is_dir():
            for path in sorted(directory.rglob("*")):
                if path.is_file() and path.resolve() != AMENDMENT.resolve():
                    files[str(path.relative_to(output))] = sha256(path)
    return {"first_source_commit": FIRST_SOURCE_COMMIT,
            "first_source_freeze_sha256": sha256(output / "manifests" / "freeze.json"),
            "first_execution_sha256": sha256(output / "manifests" / "execution.json"),
            "retained_files_sha256": files}


def correction_freeze_check(config_path: Path, output: Path, *, create: bool) -> tuple[dict, Path]:
    config = prereg_check(config_path, output)
    snapshot = first_attempt_snapshot(output)
    attempt = output / "attempts" / CORRECTION_ATTEMPT
    manifests = attempt / "manifests"
    source_path = manifests / "source.json"
    environment_path = manifests / "environment.json"
    freeze_path = manifests / "freeze.json"
    correction_path = manifests / "correction.json"
    sources = source_fingerprint(config_path)
    environment = pin_environment()
    if all(p.exists() for p in (source_path, environment_path, freeze_path, correction_path)):
        if strict_json(source_path).get("files") != sources or strict_json(environment_path).get("environment") != environment:
            raise ContractError("correction source or environment changed after freeze")
        correction = strict_json(correction_path)
        if correction.get("initial_attempt") != snapshot or correction.get("seed_schedule") != config["rng"]["formula"]:
            raise ContractError("initial attempt or seed schedule changed after correction freeze")
        frozen = strict_json(freeze_path)
        if frozen.get("status") != "SOURCE_ENVIRONMENT_FROZEN" or frozen.get("source_sha256") != sha256(source_path) or frozen.get("environment_sha256") != sha256(environment_path) or frozen.get("correction_sha256") != sha256(correction_path):
            raise ContractError("correction source/environment freeze hash mismatch")
    elif create:
        if manifests.exists() or any((attempt / p).exists() for p in ("data", "results", "figures")):
            raise ContractError("partial correction attempt already exists")
        dump_json(correction_path, {"schema": "qar-correction-v1", "attempt": CORRECTION_ATTEMPT,
                                    "reason": "public mineral_fraction data-frame requires numeric matrix conversion",
                                    "initial_attempt": snapshot, "seed_schedule": config["rng"]["formula"],
                                    "plan_sha256": sha256(PLAN), "config_sha256": sha256(config_path)})
        dump_json(source_path, {"schema": "qar-source-v1", "files": sources})
        dump_json(environment_path, {"schema": "qar-environment-v1", "environment": environment})
        dump_json(freeze_path, {"schema": "qar-source-environment-freeze-v1",
                                "status": "SOURCE_ENVIRONMENT_FROZEN", "timestamp_utc": utc(),
                                "preregistration_sha256": sha256(output / "manifests" / "preregistration_freeze.json"),
                                "source_sha256": sha256(source_path),
                                "environment_sha256": sha256(environment_path),
                                "correction_sha256": sha256(correction_path)})
    else:
        raise ContractError("correction freeze missing; run --check --attempt correction-001")
    return config, attempt


def check_mode(config_path: Path, output: Path, attempt: str | None = None) -> None:
    config = prereg_check(config_path, output)
    seed_rows(config)
    environment = pin_environment()
    with tempfile.TemporaryDirectory(prefix="qar-pinned-check-") as temp_name:
        isolated_install(Path(temp_name))
    if attempt is None:
        freeze_check(config_path, output, create=True)
    else:
        correction_freeze_check(config_path, output, create=True)
    print(json.dumps({"status": "READY", "cells": 10, "cohorts": 500,
                      "streams": 2000, "pin_tree": environment["pin_tree"],
                      "generation": "none", "attempt": attempt}, allow_nan=False))


def read_csv(path: Path, cols: tuple[str, ...]) -> list[dict]:
    with path.open(newline="") as file:
        reader = csv.DictReader(file)
        if reader.fieldnames != ["cell_id", "cohort_id", "sample_id", *cols]:
            raise ContractError(f"wrong CSV schema: {path}")
        result = []
        for row in reader:
            vals = {key: float(row[key]) for key in cols}
            if not all(math.isfinite(x) for x in vals.values()):
                raise ContractError(f"nonfinite generated value: {path}")
            result.append({"cell_id": int(row["cell_id"]), "cohort_id": int(row["cohort_id"]),
                           "sample_id": row["sample_id"], **vals})
        return result


def write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise ContractError(f"refusing to overwrite: {path}")
    with path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def group_rows(rows: list[dict]) -> dict[tuple[int, int], list[dict]]:
    groups: dict[tuple[int, int], list[dict]] = {}
    seen = set()
    for row in rows:
        key = (row["cell_id"], row["cohort_id"])
        if row["sample_id"] in seen:
            raise ContractError("duplicate sample ID")
        seen.add(row["sample_id"])
        groups.setdefault(key, []).append(row)
    return groups


def model_api():
    spec = importlib.util.spec_from_file_location("calibration_models", MODEL)
    if spec is None or spec.loader is None:
        raise ContractError("cannot load merged estimator")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def finite(value):
    if value is None or isinstance(value, (bool, str, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ContractError("estimator returned nonfinite numeric value")
        return value
    if isinstance(value, (tuple, list)):
        return [finite(x) for x in value]
    if isinstance(value, dict):
        return {k: finite(v) for k, v in value.items()}
    raise ContractError(f"non-JSON result type: {type(value)}")


def rmse(errors: list[float]) -> float:
    return math.sqrt(sum(e * e for e in errors) / len(errors))


def metrics(truth: list[dict], obs: list[dict], predictions: list[dict]) -> dict:
    if len(truth) != 40 or len(obs) != 40 or len(predictions) != 40:
        raise ContractError("cohort size mismatch")
    errors = [[p[f"mass_{phase}"] - t[phase] for phase in PHASES]
              for p, t in zip(predictions, truth)]
    all_errors = [e for row in errors for e in row]
    observed_errors = [p[f"reconstructed_{ch}"] - a[ch]
                       for p, a in zip(predictions, obs) for ch in CHANNELS]
    truth_errors = [p["reconstructed_SiO2"] - t["Quartz"] for p, t in zip(predictions, truth)]
    truth_errors += [p["reconstructed_TiO2"] - t["Anatase"] - t["Rutile"] for p, t in zip(predictions, truth)]
    result = {"mass_RMSE": rmse(all_errors),
              "totalTi_RMSE": rmse([e[1] + e[2] for e in errors]),
              "assay_RMSE_observed": rmse(observed_errors),
              "assay_RMSE_truth": rmse(truth_errors),
              "max_constraint_violation": max(p["constraint_violation"] for p in predictions),
              "max_kkt_residual": max(p["kkt_residual"] for p in predictions)}
    for j, phase in enumerate(PHASES):
        result[f"phase_RMSE_{phase}"] = rmse([e[j] for e in errors])
        result[f"phase_bias_{phase}"] = statistics.mean(e[j] for e in errors)
    return result


def matrix_eigenvalues_2x2(a: list[list[float]]) -> list[float]:
    trace = a[0][0] + a[1][1]
    delta = math.sqrt((a[0][0] - a[1][1]) ** 2 + 4 * a[0][1] ** 2)
    return [(trace - delta) / 2, (trace + delta) / 2]


def physical_diagonal(tangent: list[list[float]]) -> list[float]:
    u = ((2 / math.sqrt(6), 0), (-1 / math.sqrt(6), 1 / math.sqrt(2)),
         (-1 / math.sqrt(6), -1 / math.sqrt(2)))
    return [sum(u[j][a] * tangent[a][b] * u[j][b] for a in range(2) for b in range(2))
            for j in range(3)]


def prediction_row(cell_id: int, cohort_id: int, sample_id: str, method: str, result: dict) -> dict:
    mass = result.get("mass") or [None] * 3
    reconstructed = result.get("reconstructed") or [None] * 2
    residual = result.get("residual") or [None] * 2
    return {"cell_id": cell_id, "cohort_id": cohort_id, "sample_id": sample_id,
            "method": method, "status": result["status"],
            **{f"mass_{p}": v for p, v in zip(PHASES, mass)},
            **{f"reconstructed_{c}": v for c, v in zip(CHANNELS, reconstructed)},
            **{f"residual_{c}": v for c, v in zip(CHANNELS, residual)},
            **{k: result.get(k) for k in ("assay_loss", "prior_loss", "objective",
                                           "constraint_violation", "kkt_residual")},
            "active_support": json.dumps(result.get("active_support", []))}


def analyse(config: dict, output: Path, deadline: float) -> dict:
    api = model_api()
    data = output / "data"
    cal_obs = group_rows(read_csv(data / "calibration_observed.csv", PHASES))
    cal_truth = group_rows(read_csv(data / "calibration_truth.csv", PHASES))
    eval_obs = group_rows(read_csv(data / "evaluation_observed.csv", CHANNELS))
    eval_truth = group_rows(read_csv(data / "evaluation_truth.csv", PHASES))
    expected_keys = {(s, r) for s in range(1, 11) for r in range(1, 51)}
    if any(set(d) != expected_keys for d in (cal_obs, cal_truth, eval_obs, eval_truth)):
        raise ContractError("missing or extra cohort")
    predictions, priors, cohort_rows = [], [], []
    fit_failures = 0
    arm_failures = {m: 0 for m in config["methods"]}
    negative_measurements = 0
    support_contradictions = 0
    max_violation = max_kkt = max_fixed_difference = 0.0
    floor_count = 0
    cell_counts = {s: {"fits_failed": 0, "predictions_failed": {m: 0 for m in config["methods"]},
                       "negative_calibration_measurements": 0, "support_contradictions": 0,
                       "calibration_support_contradictions": 0,
                       "evaluation_support_contradictions": 0, "eigenvalue_floor_cohorts": 0}
                   for s in range(1, 11)}
    for cell in config["cells"]:
        s = cell["cell_id"]
        support = tuple(cell["support"])
        fixed_q = cell["fixed_quartz"]
        for r in range(1, 51):
            if time.monotonic() > deadline:
                raise TimeoutError("analysis wall-time cap")
            key = (s, r)
            y_rows, c_rows, a_rows, t_rows = (d[key] for d in (cal_obs, cal_truth, eval_obs, eval_truth))
            if len(y_rows) != cell["calibration_n"] or len(c_rows) != len(y_rows) or len(a_rows) != 40 or len(t_rows) != 40:
                raise ContractError("wrong cohort specimen count")
            if [x["sample_id"] for x in y_rows] != [x["sample_id"] for x in c_rows] or [x["sample_id"] for x in a_rows] != [x["sample_id"] for x in t_rows]:
                raise ContractError("observation/truth key mismatch")
            negatives = sum(x[p] < 0 for x in y_rows for p in PHASES)
            cal_contradictions = sum(t[p] > 1e-12 for t in c_rows for p, ok in zip(PHASES, support) if not ok)
            eval_contradictions = sum(t[p] > 1e-12 for t in t_rows for p, ok in zip(PHASES, support) if not ok)
            contradictions = cal_contradictions + eval_contradictions
            negative_measurements += negatives
            support_contradictions += contradictions
            cell_counts[s]["negative_calibration_measurements"] += negatives
            cell_counts[s]["support_contradictions"] += contradictions
            cell_counts[s]["calibration_support_contradictions"] += cal_contradictions
            cell_counts[s]["evaluation_support_contradictions"] += eval_contradictions
            R = [[cell["qxrd_sd"] ** 2 if i == j else 0.0 for j in range(3)] for i in range(3)]
            # Only signed calibration observations and declared laboratory covariance enter fit.
            prior = finite(api.fit_prior([[x[p] for p in PHASES] for x in y_rows], R,
                                         support=support, fixed_quartz=fixed_q,
                                         shrinkage=config["prior"]["shrinkage"],
                                         variance_floor=config["prior"]["variance_floor"]))
            priors.append({"cell_id": s, "cohort_id": r, "support": support,
                           "strength": cell["strength"], "fit": prior})
            fit_ok = prior["status"] == "ok"
            if not fit_ok:
                fit_failures += 1
                cell_counts[s]["fits_failed"] += 1
            else:
                raw_eigen = prior["diagnostics"].get("deconvolved_eigenvalues")
                if raw_eigen is not None and min(raw_eigen) < config["prior"]["variance_floor"]:
                    floor_count += 1
                    cell_counts[s]["eigenvalue_floor_cohorts"] += 1
            by_method = {m: [] for m in config["methods"]}
            for a_row in a_rows:
                a = [a_row[ch] for ch in CHANNELS]
                for method in config["methods"]:
                    if method != "A0" and not fit_ok:
                        result = {"status": "numerical_failure"}
                    else:
                        result = finite(api.estimate(None if method == "A1" else a,
                          None if method == "A1" else [[.0001, 0.0], [0.0, .0001]],
                          None if method == "A0" else prior, method=method,
                          support=support, fixed_quartz=fixed_q, strength=cell["strength"]))
                    row = prediction_row(s, r, a_row["sample_id"], method, result)
                    if result["status"] != "ok":
                        arm_failures[method] += 1
                        cell_counts[s]["predictions_failed"][method] += 1
                    else:
                        max_violation = max(max_violation, result["constraint_violation"])
                        max_kkt = max(max_kkt, result["kkt_residual"])
                    predictions.append(row)
                    by_method[method].append(row)
            cohort_ok = fit_ok and all(all(row["status"] == "ok" for row in by_method[m]) for m in by_method)
            if cohort_ok and fixed_q is not None and cell["strength"] > 0:
                max_fixed_difference = max(max_fixed_difference, *(abs(x[f"mass_{p}"] - y[f"mass_{p}"])
                    for x, y in zip(by_method["A1"], by_method["A2"]) for p in PHASES))
            for method, rows in by_method.items():
                values = metrics(t_rows, a_rows, rows) if cohort_ok else {}
                cohort_rows.append({"cell_id": s, "cohort_id": r, "method": method,
                                    "status": "ok" if cohort_ok else "unavailable",
                                    "reason": None if cohort_ok else ("fit_failure" if not fit_ok else "solver_failure"),
                                    **{k: values.get(k) for k in METRIC_FIELDS}})
    diagnostics = {"status": "COMPLETE" if not fit_failures and not any(arm_failures.values()) else "NUMERICAL_BLOCK",
                   "fit_failures": fit_failures, "fit_denominator": 500,
                   "prediction_failures": arm_failures, "prediction_denominator_per_arm": 20000,
                   "negative_calibration_measurements": negative_measurements,
                   "support_contradictions": support_contradictions,
                   "eigenvalue_floor_cohorts": floor_count,
                   "max_constraint_violation": max_violation, "max_kkt_residual": max_kkt,
                   "max_fixed_quartz_A1_A2_difference": max_fixed_difference,
                   "cell_counts": [{"cell_id": s, "fit_denominator": 50,
                                    "prediction_denominator_per_arm": 2000, **cell_counts[s]}
                                   for s in range(1, 11)]}
    write_csv(output / "results" / "predictions.csv", PREDICTION_FIELDS, predictions)
    prior_path = output / "results" / "priors.jsonl"
    if prior_path.exists():
        raise ContractError("refusing prior overwrite")
    prior_path.write_text("".join(json.dumps(x, allow_nan=False) + "\n" for x in priors))
    write_csv(output / "results" / "cohort_metrics.csv", COHORT_FIELDS, cohort_rows)
    summaries = paired_summaries(config, cohort_rows)
    write_csv(output / "results" / "paired_summary.csv", SUMMARY_FIELDS, summaries)
    diagnostics["calibration"] = calibration_diagnostics(config, priors, cal_truth, eval_truth)
    dump_json(output / "results" / "diagnostics.json", diagnostics)
    return diagnostics


METRIC_FIELDS = ("mass_RMSE", "phase_RMSE_Quartz", "phase_RMSE_Anatase",
                 "phase_RMSE_Rutile", "phase_bias_Quartz", "phase_bias_Anatase",
                 "phase_bias_Rutile", "totalTi_RMSE", "assay_RMSE_observed",
                 "assay_RMSE_truth", "max_constraint_violation", "max_kkt_residual")
PREDICTION_FIELDS = ["cell_id", "cohort_id", "sample_id", "method", "status",
    *(f"mass_{p}" for p in PHASES), *(f"reconstructed_{c}" for c in CHANNELS),
    *(f"residual_{c}" for c in CHANNELS), "assay_loss", "prior_loss", "objective",
    "constraint_violation", "kkt_residual", "active_support"]
COHORT_FIELDS = ["cell_id", "cohort_id", "method", "status", "reason", *METRIC_FIELDS]
SUMMARY_FIELDS = ["cell_id", "contrast", "metric", "n_planned", "n_complete",
                  "mean_difference", "descriptive_se"]


def paired_summaries(config: dict, rows: list[dict]) -> list[dict]:
    lookup = {(x["cell_id"], x["cohort_id"], x["method"]): x for x in rows}
    result = []
    for cell in config["cells"]:
        s = cell["cell_id"]
        for control in ("A0", "A1"):
            for metric in CONTRAST_METRICS:
                differences = []
                for r in range(1, 51):
                    a = lookup[s, r, "A2"]
                    b = lookup[s, r, control]
                    if a["status"] == b["status"] == "ok":
                        differences.append(a[metric] - b[metric])
                complete = len(differences) == 50
                result.append({"cell_id": s, "contrast": f"A2-A{control[-1]}",
                               "metric": metric, "n_planned": 50,
                               "n_complete": len(differences),
                               "mean_difference": statistics.mean(differences) if complete else None,
                               "descriptive_se": statistics.stdev(differences) / math.sqrt(50) if complete else None})
    return result


def calibration_diagnostics(config: dict, priors: list[dict], cal_truth: dict,
                            eval_truth: dict) -> list[dict]:
    by_cell: dict[int, list[dict]] = {}
    for p in priors:
        by_cell.setdefault(p["cell_id"], []).append(p)
    result = []
    for cell in config["cells"]:
        s = cell["cell_id"]
        ok = [p for p in by_cell[s] if p["fit"]["status"] == "ok"]
        centers = [p["fit"]["center"] for p in ok]
        cal_deltas = [[p["fit"]["center"][j] - statistics.mean(x[PHASES[j]] for x in cal_truth[s, p["cohort_id"]]) for j in range(3)] for p in ok]
        eval_deltas = [[p["fit"]["center"][j] - statistics.mean(x[PHASES[j]] for x in eval_truth[s, p["cohort_id"]]) for j in range(3)] for p in ok]
        result.append({"cell_id": s, "n_fits": len(ok), "n_planned": 50,
            "center_across_cohort_sd": [statistics.stdev(x[j] for x in centers) for j in range(3)] if len(ok) > 1 else None,
            "mean_T_physical_diagonal": [statistics.mean(physical_diagonal(p["fit"]["tangent_cov"])[j] for p in ok) for j in range(3)] if ok else None,
            "mean_M_physical_diagonal": [statistics.mean(physical_diagonal(p["fit"]["mean_cov"])[j] for p in ok) for j in range(3)] if ok else None,
            "center_minus_calibration_latent_mean": cal_deltas,
            "center_minus_heldout_cohort_truth_mean": eval_deltas,
            "mean_center_minus_calibration_latent_mean": [statistics.mean(x[j] for x in cal_deltas) for j in range(3)] if ok else None,
            "mean_center_minus_heldout_cohort_truth_mean": [statistics.mean(x[j] for x in eval_deltas) for j in range(3)] if ok else None,
            "T_eigenvalues": [matrix_eigenvalues_2x2(p["fit"]["tangent_cov"]) for p in ok],
            "M_eigenvalues": [matrix_eigenvalues_2x2(p["fit"]["mean_cov"]) for p in ok],
            "B_eigenvalues": [matrix_eigenvalues_2x2([[p["fit"]["tangent_cov"][i][j] + p["fit"]["mean_cov"][i][j] for j in range(2)] for i in range(2)]) for p in ok],
            "precision": [p["fit"]["precision"] for p in ok],
            "projection_distance": [p["fit"]["diagnostics"]["projection_distance"] for p in ok]})
    return result


def output_bytes(output: Path) -> int:
    return sum(p.stat().st_size for p in output.rglob("*") if p.is_file())


def generate(config_path: Path, output: Path, library: Path, deadline: float,
             cap_bytes: int) -> dict:
    command = ["Rscript", str(GENERATOR), str(config_path), str(output), str(library)]
    p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    reason = None
    try:
        while p.poll() is None:
            if time.monotonic() > deadline:
                p.kill()
                reason = "generation wall-time cap"
                break
            if output_bytes(output) > cap_bytes:
                p.kill()
                reason = "generation output-size cap"
                break
            time.sleep(.5)
    finally:
        if p.poll() is None:
            p.kill()
        stdout, stderr = p.communicate()
        dump_json(output / "manifests" / "generation_log.json",
                  {"command": command, "returncode": p.returncode,
                   "stdout": stdout, "stderr": stderr, "stop_reason": reason})
    if reason:
        raise ContractError(reason)
    if p.returncode:
        raise ContractError(f"R generation failed ({p.returncode}): {stderr[-3000:]}")
    if "generation_complete cells=10 cohorts=500" not in stdout:
        raise ContractError("R generation completion marker absent")
    return {"command": command, "returncode": p.returncode, "stdout": stdout,
            "stderr": stderr}


def draw_figures(output: Path) -> None:
    """Base R figures use cohort summaries; no new synthetic draws."""
    code = r'''
args <- commandArgs(TRUE)
o <- args[[1]]
d <- read.csv(file.path(o,"results","cohort_metrics.csv"))
p <- read.csv(file.path(o,"results","paired_summary.csv"))
dir.create(file.path(o,"figures"), showWarnings=FALSE)
means <- aggregate(mass_RMSE ~ cell_id + method, d, mean)
png(file.path(o,"figures","mass_rmse.png"), width=1100, height=650)
plot(NA, xlim=c(.5,10.5), ylim=range(means$mass_RMSE), xlab="Frozen cell", ylab="Mean cohort mass RMSE", main="All ten cells; independent cohorts")
for (m in c("A0","A1","A2")) { z <- means[means$method==m,]; lines(z$cell_id,z$mass_RMSE,type="b",pch=match(m,c("A0","A1","A2")),col=match(m,c("A0","A1","A2"))) }
legend("topright",legend=c("A0","A1","A2"),col=1:3,pch=1:3,lty=1)
dev.off()
z <- p[p$metric=="mass_RMSE" & p$contrast=="A2-A0",]
png(file.path(o,"figures","calibration_size.png"), width=1100, height=650)
plot(z$cell_id,z$mean_difference,type="b",pch=19,xlim=c(.5,10.5),xlab="Frozen cell (1-3 fixed; 4-6 variable)",ylab="A2-A0 mass RMSE; cohort SE",main="Paired cohort differences; descriptive SE")
arrows(z$cell_id,z$mean_difference-z$descriptive_se,z$cell_id,z$mean_difference+z$descriptive_se,angle=90,code=3,length=.04)
abline(h=0,lty=2)
dev.off()
png(file.path(o,"figures","stress.png"), width=1100, height=650)
plot(z$cell_id,z$mean_difference,type="b",pch=19,xlim=c(6.5,10.5),xlab="Stress cell (7 weak; 8 mismatch; 9 noisy; 10 wrong occurrence)",ylab="A2-A0 mass RMSE; cohort SE",main="Descriptive stress contrasts")
arrows(z$cell_id,z$mean_difference-z$descriptive_se,z$cell_id,z$mean_difference+z$descriptive_se,angle=90,code=3,length=.04)
abline(h=0,lty=2)
dev.off()
'''
    run_command(["Rscript", "-e", code, str(output)], timeout=90)


def execute(config_path: Path, output: Path, attempt: str | None = None) -> int:
    start = time.monotonic()
    if attempt is None:
        config = freeze_check(config_path, output, create=False)
    else:
        config, output = correction_freeze_check(config_path, output, create=False)
    cap_time = config["budget"]["max_wall_seconds"]
    cap_bytes = config["budget"]["max_output_bytes"]
    for part in ("data", "results", "figures"):
        if (output / part).exists():
            raise ContractError(f"execution output already exists: {part}; use a reviewed correction attempt")
    execution_path = output / "manifests" / "execution.json"
    seeds_path = output / "manifests" / "seeds.json"
    deviations_path = output / "manifests" / "deviations.json"
    if any(p.exists() for p in (execution_path, seeds_path, deviations_path)):
        raise ContractError("execution manifest already exists; refusing overwrite")
    rows = seed_rows(config)
    manifest = {"schema": "qar-execution-v1", "status": "INCOMPLETE",
                "start_utc": utc(), "end_utc": None,
                "seed_formula": config["rng"]["formula"],
                "seed_count": len(rows), "r_command": None, "r_stdout": None,
                "r_stderr": None, "failure": None, "elapsed_seconds": None,
                "output_bytes": None, "source_freeze_sha256": sha256(output / "manifests" / "freeze.json")}
    dump_json(seeds_path, {"schema": "qar-seeds-v1", "streams": rows})
    dump_json(deviations_path, {"schema": "qar-deviations-v1", "deviations": []})
    dump_json(execution_path, manifest)
    deadline = start + cap_time
    try:
        with tempfile.TemporaryDirectory(prefix="qar-pinned-run-") as temp_name:
            library = isolated_install(Path(temp_name))
            if time.monotonic() > deadline:
                raise TimeoutError("package installation exceeded wall-time cap")
            r_log = generate(config_path, output, library, deadline, cap_bytes)
            manifest["r_command"] = r_log["command"]
            manifest["r_stdout"] = r_log["stdout"]
            manifest["r_stderr"] = r_log["stderr"]
        if output_bytes(output) > cap_bytes:
            raise ContractError("output-size cap after generation")
        diagnostics = analyse(config, output, deadline)
        if output_bytes(output) > cap_bytes:
            raise ContractError("output-size cap after analysis")
        if diagnostics["status"] == "COMPLETE":
            draw_figures(output)
        if output_bytes(output) > cap_bytes or time.monotonic() > deadline:
            raise ContractError("final output/time cap")
        if diagnostics["max_constraint_violation"] > 1e-9 or diagnostics["max_kkt_residual"] > 1e-8 or diagnostics["max_fixed_quartz_A1_A2_difference"] > 1e-8:
            manifest["status"] = "NUMERICAL_BLOCK"
            manifest["failure"] = "declared numerical acceptance tolerance exceeded"
        elif diagnostics["status"] != "COMPLETE":
            manifest["status"] = "NUMERICAL_BLOCK"
            manifest["failure"] = "fit or solver failure"
        else:
            manifest["status"] = "COMPLETE"
    except BaseException as exc:
        manifest["status"] = "INCOMPLETE"
        manifest["failure"] = f"{type(exc).__name__}: {exc}"
        log_path = output / "manifests" / "generation_log.json"
        if log_path.is_file():
            log = strict_json(log_path)
            manifest["r_command"] = log["command"]
            manifest["r_stdout"] = log["stdout"]
            manifest["r_stderr"] = log["stderr"]
        raise
    finally:
        manifest["end_utc"] = utc()
        manifest["elapsed_seconds"] = time.monotonic() - start
        manifest["output_bytes"] = output_bytes(output)
        manifest["output_files_sha256"] = {
            str(p.relative_to(output)): sha256(p)
            for section in ("data", "results", "figures")
            for p in sorted((output / section).rglob("*")) if p.is_file()
        }
        dump_json(execution_path, manifest)
    return 0 if manifest["status"] == "COMPLETE" else 3


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--attempt", choices=[CORRECTION_ATTEMPT],
                        help="reviewed correction attempt in a separate output directory")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true", help="verify and freeze merged source; no draws")
    mode.add_argument("--run", action="store_true", help="execute the sole frozen pilot")
    args = parser.parse_args()
    try:
        if args.check:
            check_mode(args.config, args.output, args.attempt)
            return 0
        return execute(args.config, args.output, args.attempt)
    except (ContractError, TimeoutError, ValueError, OSError, subprocess.TimeoutExpired) as exc:
        print(f"BLOCK: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
