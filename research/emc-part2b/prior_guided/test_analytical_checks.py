"""Focused exact-arithmetic checks for the mathematical witness module."""

import importlib.util
from fractions import Fraction
from pathlib import Path
import unittest


MODULE = Path(__file__).with_name("analytical_checks.py")
SPEC = importlib.util.spec_from_file_location("analytical_checks", MODULE)
checks = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(checks)


class AnalyticalChecksTests(unittest.TestCase):
    def test_encode_reaches_same_target_in_each_order(self):
        chemistry = [[1, 1], [0, 1]]
        assay = [Fraction(3, 4), Fraction(1, 4)]
        target = [Fraction(1, 4), Fraction(1, 4)]
        for order in ([0, 1], [1, 0]):
            with self.subTest(order=order):
                fractions = checks.encode(chemistry, assay, order, target)
                result = checks.decode(chemistry, assay, order, fractions)
                self.assertEqual(result, {
                    "status": "ok",
                    "m": target,
                    "reconstructed": [Fraction(1, 2), Fraction(1, 4)],
                    "residual": [Fraction(1, 4), Fraction(0)],
                    "remaining": Fraction(1, 2),
                })

    def test_zero_capacity_and_zero_budget(self):
        chemistry = [[1, 1]]
        self.assertEqual(checks.encode(chemistry, [1], [0, 1], [1, 0]),
                         [Fraction(1), Fraction(0)])
        self.assertEqual(checks.decode(chemistry, [1], [0, 1], [1, 1])["m"],
                         [Fraction(1), Fraction(0)])
        self.assertEqual(checks.encode(chemistry, [1], [1, 0], [0, 0], budget=0),
                         [Fraction(0), Fraction(0)])

    def test_negative_assay_status_and_validation_precedence(self):
        result = checks.decode([[1]], [-Fraction(1, 100)], [0], [0])
        self.assertEqual(result, {
            "status": "unsupported_negative_assay", "m": None,
            "reconstructed": None, "residual": None, "remaining": None,
        })
        with self.assertRaises(ValueError):
            checks.decode([[1]], [-1], [0], [2])
        with self.assertRaises(ValueError):
            checks.encode([[1]], [-1], [0], [0])
        with self.assertRaises(ValueError):
            checks.encode([[1]], [0], [0], [1])

    def test_exact_score_accepts_signed_values(self):
        score = checks.quadratic_score(
            [[1, 0], [0, 1]], [-1, 1], [-2, 2], [2, 0], [0, 0], [1, 2],
            strength=Fraction(1, 2),
        )
        self.assertEqual(score, {
            "assay": Fraction(1), "prior": Fraction(3), "total": Fraction(4),
        })
        self.assertEqual(checks.quadratic_score(
            [[1]], [0], [2], [0], [5], [0], strength=0,
        )["total"], Fraction(0))

    def test_strict_types_shapes_and_nonnegative_precisions(self):
        invalid = [
            lambda: checks.decode([[0]], [0], [0], [0]),
            lambda: checks.decode([[-1]], [0], [0], [0]),
            lambda: checks.decode([[1]], [0], [True], [0]),
            lambda: checks.decode([[1]], [0], [0], [False]),
            lambda: checks.decode([[1]], [float("nan")], [0], [0]),
            lambda: checks.decode([[1]], [0], [0], [0], budget=-1),
            lambda: checks.decode([[1], [1, 1]], [0, 0], [0], [0]),
            lambda: checks.quadratic_score([[1]], [0], [0], [-1], [0], [1]),
            lambda: checks.quadratic_score([[1]], [0], [0], [1], [0], [1], -1),
        ]
        for call in invalid:
            with self.subTest(call=call), self.assertRaises(ValueError):
                call()


if __name__ == "__main__":
    unittest.main()
