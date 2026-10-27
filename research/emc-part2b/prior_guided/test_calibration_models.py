"""Focused contract and analytic QAR tests for finite calibration models."""

import importlib.util
from pathlib import Path
import unittest

import numpy as np

MODULE = Path(__file__).with_name('calibration_models.py')
SPEC = importlib.util.spec_from_file_location('calibration_models', MODULE)
models = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(models)


class CalibrationModelTests(unittest.TestCase):
    def setUp(self):
        self.y = [[0.2, 0.5, 0.3], [0.4, 0.2, 0.4],
                  [0.3, 0.4, 0.3], [0.1, 0.6, 0.3]]
        self.R = np.diag([0.001, 0.002, 0.003])
        self.prior = models.fit_prior(self.y, self.R)

    def test_prior_formula_and_unclosed_signed_input(self):
        y = np.array(self.y)
        self.assertEqual(self.prior['status'], 'ok')
        z = y @ models.U
        S = np.cov(z, rowvar=False, ddof=1)
        d, V = np.linalg.eigh(S - models.U.T @ self.R @ models.U)
        G = (V * np.maximum(d, 0.0001)) @ V.T
        T = 0.75 * G + 0.25 * np.trace(G) * np.eye(2) / 2
        M = (T + models.U.T @ self.R @ models.U) / len(y)
        np.testing.assert_allclose(self.prior['tangent_cov'], T, atol=1e-12)
        np.testing.assert_allclose(self.prior['mean_cov'], M, atol=1e-12)
        np.testing.assert_allclose(self.prior['precision'],
                                   models.U @ np.linalg.inv(T + M) @ models.U.T,
                                   atol=1e-10)
        self.assertEqual(self.prior['diagnostics']['n'], 4)
        self.assertEqual(set(self.prior), {'status', 'center', 'tangent_cov',
                                           'mean_cov', 'precision', 'diagnostics'})
        signed = models.fit_prior([[-0.2, 0.9, 0.5], [0.2, 0.7, 0.4]],
                                  np.zeros((3, 3)))
        self.assertEqual(signed['status'], 'ok')
        self.assertLess(signed['diagnostics']['raw_mean'][0], 0.01)
        self.assertGreaterEqual(min(signed['center']), -1e-9)
        self.assertAlmostEqual(sum(signed['center']), 1)

    def test_project_with_fixed_quartz_and_support(self):
        prior = models.fit_prior(self.y, self.R,
                                 support=(True, True, False), fixed_quartz=0.25)
        self.assertEqual(prior['status'], 'ok')
        np.testing.assert_allclose(prior['center'], [0.25, 0.75, 0], atol=1e-10)
        self.assertEqual(prior['diagnostics']['projection_distance'],
                         float(np.linalg.norm(np.mean(self.y, axis=0) - prior['center'])))

    def test_a0_minnorm_tie_and_signed_assay(self):
        result = models.estimate([0.2, 0.8], np.eye(2), None, method='A0')
        np.testing.assert_allclose(result['mass'], [0.2, 0.4, 0.4], atol=1e-10)
        self.assertLess(result['kkt_residual'], 1e-8)
        fixed = models.estimate([-2, 3], np.eye(2), None, method='A0', fixed_quartz=0.3)
        np.testing.assert_allclose(fixed['mass'], [0.3, 0.35, 0.35], atol=1e-10)
        one = models.estimate([0.2, 0.8], np.eye(2), None, method='A0',
                              support=(True, False, True))
        np.testing.assert_allclose(one['mass'], [0.2, 0, 0.8], atol=1e-10)

    def test_a1_no_assay_and_a2_global_faces(self):
        a1 = models.estimate(None, None, self.prior, method='A1')
        self.assertIsNone(a1['residual'])
        self.assertIsNone(a1['assay_loss'])
        np.testing.assert_allclose(a1['mass'], self.prior['center'])
        a2 = models.estimate([0.8, 0.2], np.eye(2), self.prior, method='A2')
        self.assertEqual(a2['status'], 'ok')
        self.assertLess(a2['constraint_violation'], 1e-9)
        self.assertLess(a2['kkt_residual'], 1e-8)
        zero = models.estimate([0.8, 0.2], np.eye(2), self.prior,
                               method='A2', strength=0)
        a0 = models.estimate([0.8, 0.2], np.eye(2), None, method='A0')
        np.testing.assert_allclose(zero['mass'], a0['mass'])
        fixed_prior = models.fit_prior(self.y, self.R, fixed_quartz=0.2)
        fixed_a1 = models.estimate(None, None, fixed_prior, method='A1',
                                   fixed_quartz=0.2)
        fixed_a2 = models.estimate([0.9, 0.1], np.eye(2), fixed_prior,
                                   method='A2', fixed_quartz=0.2, strength=3)
        np.testing.assert_allclose(fixed_a2['mass'], fixed_a1['mass'], atol=1e-8)

    def test_boundary_and_degenerate_feasible_set(self):
        prior = models.fit_prior(self.y, self.R,
                                 support=(True, False, True), fixed_quartz=0)
        for method in ('A0', 'A2'):
            result = models.estimate([1.7, -1.4], np.diag([0.02, 0.03]),
                                     None if method == 'A0' else prior,
                                     method=method, support=(True, False, True),
                                     fixed_quartz=0)
            self.assertEqual(result['status'], 'ok')
            np.testing.assert_allclose(result['mass'], [0, 0, 1], atol=1e-8)
            self.assertLess(result['kkt_residual'], 1e-8)
        only_q = models.estimate([-1, 2], np.eye(2), None, method='A0',
                                 support=(True, False, False))
        np.testing.assert_allclose(only_q['mass'], [1, 0, 0], atol=1e-10)

    def test_invalid_inputs(self):
        bad = [
            lambda: models.fit_prior([[1, 0, 0]], self.R),
            lambda: models.fit_prior([[True, 0, 0], [0, 1, 0]], self.R),
            lambda: models.fit_prior(self.y, np.diag([-0.1, 1, 1])),
            lambda: models.fit_prior(self.y, self.R, support=(True, 1, True)),
            lambda: models.fit_prior(self.y, self.R, fixed_quartz=1.1),
            lambda: models.fit_prior(self.y, self.R, variance_floor=0),
            lambda: models.estimate([0, 1], np.eye(2), None, method='A2'),
            lambda: models.estimate(None, np.eye(2), None, method='A0'),
            lambda: models.estimate([0, 1], np.zeros((2, 2)), None, method='A0'),
            lambda: models.estimate([0, 1], np.eye(2), self.prior, method='A0'),
            lambda: models.estimate([0, 1], np.eye(2), self.prior, method='A2',
                                    fixed_quartz=0.1),
            lambda: models.estimate([0, 1], np.eye(2), self.prior, method='A2',
                                    strength=-1),
        ]
        for call in bad:
            with self.subTest(call=call), self.assertRaises(ValueError):
                call()


if __name__ == '__main__':
    unittest.main()
