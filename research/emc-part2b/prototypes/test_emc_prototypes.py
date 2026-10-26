"""Focused contract checks for the proposed standalone EMC prototypes."""

import importlib.util
from pathlib import Path
import unittest

import numpy as np


MODULE = Path(__file__).with_name("emc_prototypes.py")
SPEC = importlib.util.spec_from_file_location("emc_prototypes", MODULE)
emc = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(emc)


class PrototypeChecks(unittest.TestCase):
    def test_allocation_and_original_residual(self):
        C = np.array([[1.0, 0.0], [0.0, 1.0]])
        result = emc.allocate(C, np.array([0.2, 0.3]), (0, 1))
        np.testing.assert_allclose(result["m"], [0.2, 0.3])
        np.testing.assert_allclose(result["residual"], [0.0, 0.0])
        self.assertAlmostEqual(result["unallocated"], 0.5)
        self.assertEqual(result["status"], "ok")
        tiny = emc.allocate(C, np.array([-5e-13, 0.2]), (0, 1))
        self.assertAlmostEqual(tiny["correction_max"], 5e-13)
        self.assertAlmostEqual(tiny["residual"][0], -5e-13)

    def test_signed_assay_is_unsupported_not_invalid(self):
        C = np.array([[1.0]])
        result = emc.allocate(C, np.array([-1e-4]), (0,))
        self.assertEqual(result["status"], "unsupported_negative_assay")
        self.assertTrue(np.isnan(result["m"]).all())
        fit = emc.fit_order(C, np.array([[0.2], [-1e-4]]))
        self.assertIsNone(fit["order"])
        self.assertEqual(fit["proposals"], 0)
        with self.assertRaises(ValueError):
            emc.allocate(C, [np.nan], (0,))

    def test_score_components_and_single_channel_covariance(self):
        A = np.array([[0.0], [1.0]])
        same = emc.score_assays(A, A, objective="distribution")
        for name in ("row", "mean", "covariance", "js", "total"):
            self.assertAlmostEqual(same[name], 0.0)
        shifted = emc.score_assays(A, np.array([[0.5], [0.5]]), objective="distribution")
        self.assertAlmostEqual(shifted["row"], 0.25)
        self.assertGreater(shifted["covariance"], 0.0)
        self.assertGreater(shifted["js"], 0.0)
        with self.assertRaises(ValueError):
            emc.score_assays(A, A, scales=[0.0])

    def test_search_methods_and_cache_accounting(self):
        C = np.array([[1.0, 0.5], [0.0, 0.5]])
        A = np.array([[0.3, 0.2], [0.4, 0.1], [0.2, 0.3]])
        exhaustive = emc.fit_order(C, A, method="exhaustive", budget=1)
        self.assertEqual(exhaustive["proposals"], 2)
        self.assertEqual(exhaustive["evaluations"], 2)
        for method in ("fixed", "greedy", "beam", "random", "ga"):
            result = emc.fit_order(C, A, method=method, budget=4, seed=7)
            self.assertEqual(result["status"], "ok")
            self.assertLessEqual(result["proposals"], 4)
            self.assertLessEqual(result["evaluations"], result["proposals"])
            self.assertEqual(set(result["order"]), {0, 1})
            self.assertEqual(len(result["trace"]), result["proposals"])
            self.assertAlmostEqual(result["train_score"], result["components"]["total"])
        prediction = emc.predict_order(C, A, exhaustive["order"])
        self.assertEqual(prediction["A_hat"].shape, A.shape)
        self.assertEqual(prediction["M"].shape, (3, 2))
        with self.assertRaises(ValueError):
            emc.fit_order(C, A, budget=1.5)
        with self.assertRaises(ValueError):
            emc.predict_order(C, A, (0, 0))


if __name__ == "__main__":
    unittest.main()
