"""Deterministic harness checks; no generator or scientific draw is invoked."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

import measurement_experiment as study


def args():
    return (study.HERE / "preregistration.json", study.HERE / "PREEXECUTION_REVIEW.json")


def test_design_gate_exact_bytes_and_review(tmp_path):
    cfg, review, freeze = study.design_gate(*args())
    assert cfg["replicates"] == 50
    assert freeze["design_digest"] == review["design_digest"]
    assert study.planned_keys(cfg)["predictions"] == 350250
    assert study.planned_keys(cfg)["cohort_metrics"] == 21250
    assert review["label_diff_verdict"] == "PASS"


def test_seed_inventory_and_disjointness():
    manifest = study.seed_manifest()
    assert len(manifest["streams"]) == 2000
    assert manifest["collision_count"] == 0
    assert any("finite_calibration/manifests/seeds.json" in x["path"]
               for x in manifest["historical_manifests"])
    assert 20261027 in manifest["historical_seed_set"]


def test_covariance_channels_and_full_stack():
    assert np.linalg.matrix_rank(study.C) == 2
    assert np.linalg.matrix_rank(np.vstack([study.C, study.H])) == 3
    assert np.linalg.eigvalsh(study.R4).min() > 0
    observed = {"SiO2":".4", "TiO2":".6", "rep_SiO2":".39",
                "rep_TiO2":".61", "QXRD_Q":".4", "QXRD_A":".2", "QXRD_R":".4"}
    D, y, cov = study.design("P", observed, True)
    assert D.shape == (5,3) and y.shape == (5,) and cov.shape == (5,5)
    assert np.allclose(cov[:2,2:], 0)
    D, y, cov = study.design("R", observed, True)
    assert D.shape == (4,3) and np.allclose(cov[:2,2:], 0)


def test_no_truth_in_fit_path_static():
    import inspect
    src = inspect.getsource(study.run_predictions)
    assert "_truth.csv" not in src
    assert 'f"{kind}_observed.csv"' in src
    score_src = inspect.getsource(study.score)
    assert "_truth.csv" in score_src


def test_r_truth_matrix_conversion_precedes_finite_check():
    source = (study.HERE / "generation.R").read_text()
    assert source.index("x <- as.matrix(raw)") < source.index("storage.mode(x) <- \"double\"")
    assert source.index("storage.mode(x) <- \"double\"") < source.index("any(!is.finite(x))")


def test_failed_predictions_block_matched_metric():
    sample = {"sample_id":"u1", "status":"numerical_failure"}
    result = study.score_block([sample], {}, {}, [])
    assert result["status"] == "unavailable_matched_failure"
    assert result["n_planned"] == 1 and result["n_failure"] == 1
    assert study.score_block([], {}, {}, [])["status"] == "empty_stratum"


def test_metric_recalculation():
    p = [{"sample_id":"u1","status":"ok","Q_hat":".5","A_hat":".2","R_hat":".3","assay_loss":".01","observation_design_id":"B_u1"},
         {"sample_id":"u2","status":"ok","Q_hat":".6","A_hat":".1","R_hat":".3","assay_loss":".02","observation_design_id":"B_u2"}]
    t = {"u1":{"Q":".4","A":".3","R":".3"},
         "u2":{"Q":".6","A":".2","R":".2"}}
    o = {"u1":{"SiO2":".41","TiO2":".59"},
         "u2":{"SiO2":".59","TiO2":".41"}}
    b = [{"status":"ok","width_Q":".1","width_A":".2","width_R":".3","observation_design_id":"B_u1"},
         {"status":"infeasible_envelope","observation_design_id":"B_u2"}]
    x = study.score_block(p,t,o,b)
    assert x["status"] == "ok"
    assert x["mass_RMSE"] == pytest.approx(np.sqrt(.04/6))
    assert x["phase_bias_A"] == pytest.approx(-.1)
    assert x["n_envelope_nonempty"] == 1 and x["n_envelope_infeasible"] == 1
    assert x["mean_width_A"] == pytest.approx(.2)


def test_prediction_failure_preserves_keys():
    keys = {"cell_id":1,"cohort_id":1,"support_group":"full","stage":2,
            "budget":5,"stratum":"unmeasured","arm":"A2_P","sample_id":"u1"}
    x = study.prediction_row(keys,{"status":"numerical_failure","reason":"oracle"},
                             "prior1","B_u1",np.array([.4,.6]))
    assert x["status"] == "numerical_failure"
    assert x["sample_id"] == "u1" and "Q_hat" not in x


def test_gate_rejects_review_changes_without_sampling(tmp_path, monkeypatch):
    cfg, review, freeze = study.design_gate(*args())
    bad = dict(review, verdict="PASS WITH NOTE")
    review_path = tmp_path / "review.json"
    review_path.write_text(json.dumps(bad))
    with pytest.raises(RuntimeError):
        study.design_gate(args()[0], review_path)


def test_run_requires_all_tester_certificates(tmp_path):
    path = tmp_path / "audit" / "checks.json"
    path.parent.mkdir()
    path.write_text(json.dumps({"cases":[{"check_id":"D01","command":"deterministic",
        "expected":"PASS","actual":"PASS","verdict":"PASS","evidence_path":"audit/evidence.txt"}]}))
    (path.parent / "evidence.txt").write_text("checked")
    with pytest.raises(RuntimeError, match="tester certification incomplete"):
        study.certification_gate(tmp_path)


def test_certification_phase_boundary_without_draws(tmp_path):
    audit = tmp_path / "audit"
    audit.mkdir()
    (audit / "evidence.txt").write_text("deterministic evidence")
    cases = [{"check_id": check_id, "command": "deterministic", "expected": "PASS",
              "actual": "PASS", "verdict": "PASS", "evidence_path": "audit/evidence.txt"}
             for check_id in sorted(study.PRE_CASES)]
    cases += [{"check_id": check_id, "command": "after run", "expected": "PASS",
               "actual": None, "verdict": "PENDING", "evidence_path": None}
              for check_id in sorted(study.POST_CASES)]
    path = audit / "checks.json"
    path.write_text(json.dumps({"cases": cases}))
    study.certification_gate(tmp_path)
    with pytest.raises(RuntimeError, match="unexecuted/failed certification D13"):
        study.certification_gate(tmp_path, final=True)
    for case in cases:
        if case["check_id"] in study.POST_CASES:
            case.update(actual="PASS", verdict="PASS", evidence_path="audit/evidence.txt")
    path.write_text(json.dumps({"cases": cases}))
    study.certification_gate(tmp_path, final=True)
    with pytest.raises(RuntimeError, match="post-run certification executed before run"):
        study.certification_gate(tmp_path)
