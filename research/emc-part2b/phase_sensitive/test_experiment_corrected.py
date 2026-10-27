"""Focused draw-free tests for the additive post-execution harness."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

import measurement_experiment as old
import measurement_experiment_corrected as corrected


class CorrectedHarnessTests(unittest.TestCase):
    def test_elapsed_preflight_rejects_before_callback(self):
        with tempfile.TemporaryDirectory() as directory:
            now = [0.]
            invoked = []
            def preflight(gate):
                now[0] = 10.
            with self.assertRaisesRegex(RuntimeError, "wall-time cap"):
                corrected.non_generating_boundary(Path(directory), lambda: invoked.append(1),
                                                   preflight=preflight, clock=lambda: now[0], seconds=10)
            self.assertEqual(invoked, [])

    def test_existing_bytes_and_intervening_write_reject(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            (output / "old-ledger").write_bytes(b"123456")
            calls = []
            with self.assertRaisesRegex(RuntimeError, "byte cap"):
                corrected.non_generating_boundary(output, lambda: calls.append(1),
                                                   preflight=lambda gate: None, byte_cap=5)
            self.assertEqual(calls, [])
            def preflight(gate):
                (output / "new-evidence").write_bytes(b"12345")
            with self.assertRaisesRegex(RuntimeError, "byte cap"):
                corrected.non_generating_boundary(output, lambda: calls.append(1),
                                                   preflight=preflight, byte_cap=10)
            self.assertEqual(calls, [])

    def test_exact_deadline_rejects(self):
        with tempfile.TemporaryDirectory() as directory:
            now = [0.]
            called = []
            def preflight(gate):
                now[0] = 2.
            with self.assertRaisesRegex(RuntimeError, "wall-time"):
                corrected.non_generating_boundary(directory, lambda: called.append(1),
                                                   preflight=preflight, clock=lambda: now[0], seconds=2)
            self.assertFalse(called)

    def test_mask_retains_only_accounting(self):
        keys = {"cell_id": 1, "arm": "A2_P"}
        rows = [{"status": "ok", "mass_RMSE": .2}, {"status": "numerical_failure"}]
        masked = corrected.masked_record(keys, rows, 2)
        self.assertEqual((masked["n_planned"], masked["n_success"], masked["n_failure"]), (2, 1, 1))
        self.assertEqual(set(masked), set(keys) | {"status", "reason", "n_planned", "n_success", "n_failure"})

    def test_independent_score_fixture(self):
        row = {"status": "ok", "sample_id": "unit-1", "Q_hat": ".6", "A_hat": ".2",
               "R_hat": ".2", "measurement_design_sse": ".02",
               "measurement_design_n": "2", "observation_design_id": "B-unit-1"}
        truth = {"unit-1": {"Q": ".5", "A": ".2", "R": ".3"}}
        observed = {"unit-1": {"SiO2": ".52", "TiO2": ".48"}}
        bounds = [{"observation_design_id": "B-unit-1", "status": "ok",
                   "width_Q": ".1", "width_A": ".4", "width_R": ".3"}]
        result = corrected.score_block_corrected([row], truth, observed, bounds)
        self.assertAlmostEqual(result["mass_RMSE"], np.sqrt(.02 / 3))
        self.assertAlmostEqual(result["totalTi_RMSE"], .1)
        self.assertAlmostEqual(result["bulk_assay_RMSE_observed"], .08)
        self.assertAlmostEqual(result["bulk_assay_RMSE_latent"], .1)
        self.assertAlmostEqual(result["measurement_design_RMSE_observed"], .1)
        self.assertEqual(result["n_envelope_nonempty"], 1)
        self.assertAlmostEqual(result["mean_width_A"], .4)

    def test_fake_reported_kkt_does_not_authorize_wrong_mass(self):
        D, R = old.C, old.RB
        y = np.array([.62, .39])
        m = np.array([.99, .005, .005])
        inverse = np.linalg.inv(R)
        H = D.T @ inverse @ D
        f = D.T @ inverse @ y
        with self.assertRaisesRegex(RuntimeError, "KKT"):
            corrected.kkt_certificate(m, H, f, (True, True, True))

    def test_hard_zero_is_equality(self):
        H = np.eye(3)
        m = np.array([.6, 0., .4])
        f = m.copy()
        f[1] = 100.  # A gradient is irrelevant on the hard-zero coordinate.
        corrected.kkt_certificate(m, H, f, (True, False, True))

    def test_output_overlap_rejected(self):
        with self.assertRaisesRegex(RuntimeError, "overlap"):
            corrected.safe_paths(
                "/home/matt/Documents/GitHub/writing/research/emc-part2b/phase_sensitive/run/attempts/correction-002",
                "/home/matt/Documents/GitHub/writing/research/emc-part2b/phase_sensitive/run/attempts/correction-002/evidence")


if __name__ == "__main__":
    unittest.main()
