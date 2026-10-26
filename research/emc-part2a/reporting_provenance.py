#!/usr/bin/env python3
"""Add row provenance to immutable EMC result CSVs without recomputing fits.

Run from any directory. Part B is processed only when its metrics file exists;
the source CSVs and their timing values are never rewritten.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DEPENDENCY = "d9a0b27676a0cb3a4acb300a3f6b2ebc845d205d"
UNITS = "whole-sample mass fraction"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def enrich(source: Path, target: Path, sidecar: Path, freeze_hash: str) -> int:
    with source.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        if not fields or ("units" in fields) != ("dependency_hash" in fields):
            raise ValueError(f"Unexpected source columns: {source}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"Empty source CSV: {source}")
    target.parent.mkdir(parents=True, exist_ok=True)
    already_enriched = "units" in fields
    if already_enriched and any(row["units"] != UNITS or row["dependency_hash"] != DEPENDENCY for row in rows):
        raise ValueError(f"Conflicting row provenance: {source}")
    with target.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields if already_enriched else fields + ["units", "dependency_hash"])
        writer.writeheader()
        for row in rows:
            writer.writerow(row if already_enriched else {**row, "units": UNITS, "dependency_hash": DEPENDENCY})
    metadata = {
        "schema": "emc-reporting-provenance-1",
        "source": str(source.relative_to(ROOT.parent)),
        "derived": str(target.relative_to(ROOT.parent)),
        "source_sha256": sha256(source),
        "derived_sha256": sha256(target),
        "baseline_freeze_sha256": freeze_hash,
        "dependency_hash": DEPENDENCY,
        "rows": len(rows),
        "units": UNITS,
        "column_units": {
            "solver_objective": "squared standardized assay mass fraction for weighted methods; squared assay mass fraction otherwise",
            "assay_rmse": "whole-sample mass fraction",
            "standardized_assay_rmse": "dimensionless",
            "mineral_rmse": "whole-sample mass fraction",
            "closure_error": "whole-sample mass fraction",
            "unallocated": "whole-sample mass fraction",
            "negative_mass": "whole-sample mass fraction",
            "exact_residual": "whole-sample mass fraction",
            "fit_wall_seconds": "seconds",
            "fit_cpu_seconds": "seconds",
            "wall_seconds": "seconds",
            "cpu_seconds": "seconds",
            "row_loss": "proposed auxiliary score; dimensionless after scaling",
            "mean_loss": "proposed auxiliary score; dimensionless after scaling",
            "covariance_loss": "proposed auxiliary score; dimensionless after scaling",
            "js_loss": "dimensionless",
            "total_loss": "proposed auxiliary score; dimensionless after scaling",
            "train_score": "proposed auxiliary score; dimensionless after scaling",
            "objective_gap": "proposed auxiliary score; dimensionless after scaling",
            "n_total": "samples",
            "n_supported": "samples",
            "unsupported_count": "samples",
            "failure_count": "samples",
            "evaluations": "candidate evaluations",
            "proposals": "candidate proposals",
        },
        "derivation": "CSV text fields copied unchanged; units and dependency_hash appended only if absent",
    }
    sidecar.write_text(json.dumps(metadata, indent=2) + "\n")
    return len(rows)


def main() -> None:
    freeze = ROOT / "baseline_freeze.json"
    frozen = json.loads(freeze.read_text())
    if frozen["dependency_commit"] != DEPENDENCY:
        raise ValueError("Pinned dependency differs from current A freeze")
    for relative, expected in frozen["hashes"].items():
        actual = sha256(ROOT / relative)
        if actual != expected:
            raise ValueError(f"Frozen A hash mismatch: {relative}")
    freeze_hash = sha256(freeze)
    a_source = ROOT / "results/baselines.csv"
    count = enrich(a_source, ROOT / "results/baselines_provenance.csv",
                   ROOT / "results/baselines_provenance.json", freeze_hash)
    print(f"A: {count} provenance rows; freeze references: {len(frozen['hashes'])}")
    b_source = ROOT.parent / "emc-part2b/results/metrics.csv"
    if b_source.exists():
        count = enrich(b_source, b_source.with_name("metrics_provenance.csv"),
                       b_source.with_name("metrics_provenance.json"), freeze_hash)
        print(f"B: {count} provenance rows")


if __name__ == "__main__":
    main()
