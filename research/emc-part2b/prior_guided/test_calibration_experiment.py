"""Draw-free checks of the frozen finite-calibration harness contract."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("calibration_experiment", ROOT / "calibration_experiment.py")
assert SPEC is not None and SPEC.loader is not None
harness = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(harness)


class CheckOnlyTests(unittest.TestCase):
    def setUp(self) -> None:
        self.output = ROOT / "finite_calibration"
        self.config_path = self.output / "preregistration.json"

    def test_preregistration_hashes_and_unique_streams(self) -> None:
        config = harness.prereg_check(self.config_path, self.output)
        streams = harness.seed_rows(config)
        self.assertEqual(len(streams), 2000)
        self.assertEqual(len({x["seed"] for x in streams}), 2000)
        self.assertEqual(streams[0]["seed"], 20361128)
        self.assertEqual(streams[-1]["seed"], 21266031)

    def test_frozen_grid_refuses_extra_cell(self) -> None:
        config = harness.strict_json(self.config_path)
        config["cells"].append(config["cells"][-1].copy())
        with self.assertRaises(harness.ContractError):
            harness.validate_config(config)

    def test_metric_formulas_use_cohort_denominator(self) -> None:
        truth = [{"Quartz": .6, "Anatase": .1, "Rutile": .3} for _ in range(40)]
        obs = [{"SiO2": .6, "TiO2": .4} for _ in range(40)]
        predictions = [{"mass_Quartz": .6, "mass_Anatase": .2,
                        "mass_Rutile": .2, "reconstructed_SiO2": .6,
                        "reconstructed_TiO2": .4, "constraint_violation": 0,
                        "kkt_residual": 0} for _ in range(40)]
        result = harness.metrics(truth, obs, predictions)
        self.assertAlmostEqual(result["mass_RMSE"], (2 * .1**2 / 3) ** .5)
        self.assertAlmostEqual(result["totalTi_RMSE"], 0)
        self.assertAlmostEqual(result["assay_RMSE_truth"], 0)


if __name__ == "__main__":
    unittest.main()
