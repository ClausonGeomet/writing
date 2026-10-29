#!/usr/bin/env python3
"""Reproduce compact EMC calculations with frozen Part 2A/B research functions.

This generates no additional scientific experiment. Run generate_data.R first.
Import the two unchanged frozen Python research modules alongside this file.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from baselines import solve
from emc_prototypes import allocate, fit_order, predict_order, score_assays

PHASES = ("Quartz", "Hematite", "Kaolinite", "Forsterite")
ANALYTES = ("SiO2", "Fe2O3_total", "Al2O3", "MgO")

def run(directory: Path) -> None:
    frame = pd.read_csv(directory / "emc_data.csv").sort_values("sample_id").reset_index(drop=True)
    response = pd.read_csv(directory / "emc_C.csv").set_index("analyte")
    C = response.loc[list(ANALYTES), list(PHASES)].to_numpy(dtype=float)
    M = frame[[f"m_{p}" for p in PHASES]].to_numpy(dtype=float)
    A = frame[[f"a_{a}" for a in ANALYTES]].to_numpy(dtype=float)
    if len(frame) != 60:
        raise ValueError("Expected 60 specimens from the R generator")
    np.testing.assert_allclose(M.sum(axis=1), 1.0, atol=1e-10)
    np.testing.assert_allclose(A, M @ C.T, atol=1e-10)

    # Theoretical forward moment identities, checked using simulated truth.
    mu_error = np.max(np.abs(A.mean(axis=0) - C @ M.mean(axis=0)))
    cov_error = np.max(np.abs(np.cov(A, rowvar=False) - C @ np.cov(M, rowvar=False) @ C.T))
    print(f"Forward moments: max mean error {mu_error:.3e}; covariance error {cov_error:.3e}")

    # Fixed greedy order can change the allocation on an identical assay.
    demo_m = np.array([.45, .15, .20, .20])
    demo_a = C @ demo_m
    quartz_first = allocate(C, demo_a, order=(0, 1, 2, 3))
    silicates_first = allocate(C, demo_a, order=(2, 3, 0, 1))
    print("Quartz-first order:", np.round(quartz_first["m"], 4))
    print("Silicates-first order:", np.round(silicates_first["m"], 4))
    np.testing.assert_allclose(silicates_first["m"], demo_m, atol=1e-8)

    # One 30/10/20 split, consistent in spirit with the frozen main design.
    train, validation, test = A[:30], A[30:40], A[40:]
    truth = M[40:]
    baseline = [solve(C, a, method="cls_eq") for a in test]
    if any(row["status"] != "ok" for row in baseline):
        raise RuntimeError("Constrained least-squares baseline failed")
    M_cls = np.vstack([row["m"] for row in baseline])
    exhaustive = fit_order(C, train, method="exhaustive", objective="distribution")
    ga = fit_order(C, train, method="ga", objective="distribution",
                   seed=101, budget=24, beam_width=2)
    pred_exh = predict_order(C, test, exhaustive["order"])
    pred_ga = predict_order(C, test, ga["order"])
    for label, estimate in [
        ("Constrained least squares", M_cls),
        ("Exhaustive-order allocation", pred_exh["M"]),
        ("GA-order allocation", pred_ga["M"]),
    ]:
        rmse = np.sqrt(np.mean((estimate - truth)**2))
        print(f"{label:30s} held-out mineral RMSE: {rmse:.6f}")
    print("Exhaustive order/score:", exhaustive["order"], exhaustive["train_score"])
    print("GA order/score:", ga["order"], ga["train_score"])
    print("GA unique evaluations / proposals:", ga["evaluations"], ga["proposals"])
    assert exhaustive["train_score"] <= ga["train_score"] + 1e-10

    # The actual fixed experimental scorer decomposes observable assay loss.
    parts = score_assays(test, pred_ga["A_hat"], objective="distribution",
                         scales=ga["scales"])
    print("Distribution-score components:",
          {k: round(float(parts[k]), 6) for k in
           ("row", "mean", "covariance", "js", "total")})

    # Stress-test different failure modes on the same independent evaluation rows.
    # These are teaching perturbations, NOT the frozen six-scenario benchmark.
    rng = np.random.default_rng(20261029)
    wrong_C = C.copy()
    wrong_C[2, 2] *= 0.8   # Fitted kaolinite Al2O3 differs from the generator.
    conditions = [
        ("ideal", C, test, "cls_eq"),
        ("noisy assay", C, test + rng.normal(size=test.shape)
         * np.array([.010, .008, .003, .006]), "cls_eq"),
        ("missing forsterite", C[:, :3], test, "cls_le"),
        ("wrong mineral chemistry", wrong_C, test, "cls_eq"),
        ("infeasible assay x1.2", C, 1.2 * test, "cls_eq"),
    ]
    print("One-cohort exploratory stress checks (not the historical study):")
    for label, fit_C, observed, method in conditions:
        fitted = [solve(fit_C, a, method=method) for a in observed]
        if any(row["status"] != "ok" for row in fitted):
            raise RuntimeError(f"Fit did not converge: {label}")
        estimate = np.zeros_like(truth)
        estimate[:, :fit_C.shape[1]] = np.vstack([row["m"] for row in fitted])
        fitted_assays = estimate[:, :fit_C.shape[1]] @ fit_C.T
        mineral_rmse = np.sqrt(np.mean((estimate - truth) ** 2))
        assay_rmse = np.sqrt(np.mean((fitted_assays - observed) ** 2))
        print(f"  {label:25s} mineral RMSE={mineral_rmse:.6f}; assay RMSE={assay_rmse:.6f}")

    # Two distinct mineral populations, exactly the same assay distribution.
    C_qar = np.array([[1., 0., 0.], [0., 1., 1.]])
    M_qar = np.array([[.6, .25, .15], [.5, .1, .4],
                      [.7, .15, .15], [.4, .35, .25]])
    assay_qar = M_qar @ C_qar.T
    assay_swap = M_qar[:, [0, 2, 1]] @ C_qar.T
    np.testing.assert_allclose(assay_qar, assay_swap, atol=1e-12)
    identical = score_assays(assay_qar, assay_swap, objective="distribution")
    assert identical["total"] == 0.0
    print("QAR phase swap: identical assay observations and zero JS divergence")
    print("Teaching example completed; for all six scenarios consult the frozen study.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", nargs="?", type=Path, default=Path("."))
    run(parser.parse_args().directory)
