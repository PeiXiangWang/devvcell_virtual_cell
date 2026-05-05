from __future__ import annotations

import sys
import unittest
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from devvcell.response_recovery import classify_from_latent_tables, summarize_response_recovery  # noqa: E402


class ResponseRecoveryClassificationTest(unittest.TestCase):
    def test_delay_and_fate_classes_are_detected_from_latent_geometry(self) -> None:
        centroids = pd.DataFrame(
            [
                {"stage": "Theiler stage 14", "stage_num": 14, "cell_type": "A", "latent_01": 0.0, "latent_02": 0.0},
                {"stage": "Theiler stage 15", "stage_num": 15, "cell_type": "A", "latent_01": 1.0, "latent_02": 0.0},
                {"stage": "Theiler stage 16", "stage_num": 16, "cell_type": "A", "latent_01": 2.0, "latent_02": 0.0},
                {"stage": "Theiler stage 15", "stage_num": 15, "cell_type": "B", "latent_01": 1.0, "latent_02": 1.0},
                {"stage": "Theiler stage 16", "stage_num": 16, "cell_type": "B", "latent_01": 2.0, "latent_02": 1.0},
            ]
        )
        transferred = pd.DataFrame(
            [
                {
                    "perturbation": "delay",
                    "stage": "Theiler stage 15",
                    "cell_type": "A",
                    "response_latent_01": -0.9,
                    "response_latent_02": 0.0,
                },
                {
                    "perturbation": "fate",
                    "stage": "Theiler stage 15",
                    "cell_type": "A",
                    "response_latent_01": 0.0,
                    "response_latent_02": 0.9,
                },
            ]
        )
        config = {
            "classification_thresholds": {
                "delay_margin": 0.05,
                "fate_margin": 0.05,
                "off_manifold_quantile": 0.95,
                "off_manifold_multiplier": 5.0,
            }
        }
        classes = classify_from_latent_tables(transferred, centroids, config)
        labels = dict(zip(classes["perturbation"], classes["response_recovery_class"]))
        self.assertEqual(labels["delay"], "developmental_delay")
        self.assertEqual(labels["fate"], "fate_deflection")

        summary = summarize_response_recovery(classes)
        self.assertEqual(int(summary["n_cases"].sum()), 2)


if __name__ == "__main__":
    unittest.main()
