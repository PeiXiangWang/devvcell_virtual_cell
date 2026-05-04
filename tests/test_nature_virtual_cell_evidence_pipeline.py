from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts import run_nature_virtual_cell_evidence_pipeline as pipeline  # noqa: E402


def _row_for_table(table_name: str, model: str | None = None) -> dict[str, object]:
    spec = pipeline.REQUIRED_COLUMNS[table_name]
    numeric = set(spec["numeric"])
    row: dict[str, object] = {}
    for column in spec["columns"]:
        if column in numeric:
            row[column] = 1.0
        else:
            row[column] = f"{column}_value"

    if "model" in row and model is not None:
        row["model"] = model
    if "split" in row:
        row["split"] = "test"
    if "system" in row:
        row["system"] = "hindbrain"
    if "condition" in row:
        row["condition"] = "guide_A"
    if "stage" in row:
        row["stage"] = "TS14"
    if "stage_num" in row:
        row["stage_num"] = 14
    if "src_stage" in row:
        row["src_stage"] = 14
    if "tgt_stage" in row:
        row["tgt_stage"] = 15
    if "tf_name" in row:
        row["tf_name"] = "TfA"
    if "best_rescuer_tf" in row:
        row["best_rescuer_tf"] = "TfB"
    return row


def _write_required_tables(root: Path, mutate: callable | None = None) -> dict[str, str]:
    input_tables: dict[str, str] = {}
    for table_name in pipeline.REQUIRED_COLUMNS:
        rows = [_row_for_table(table_name), _row_for_table(table_name)]
        if table_name == "transition_metrics":
            rows = [
                _row_for_table(table_name, "context_residual_mlp"),
                _row_for_table(table_name, "identity"),
                _row_for_table(table_name, "mean_shift"),
                _row_for_table(table_name, "ridge"),
                _row_for_table(table_name, "mlp"),
            ]
        if table_name == "external_condition_metrics":
            rows = [
                _row_for_table(table_name, "gene_context_ridge_residual"),
                _row_for_table(table_name, "identity"),
                _row_for_table(table_name, "global_mean_shift"),
                _row_for_table(table_name, "gene_context_shift"),
                _row_for_table(table_name, "gene_context_mlp_residual"),
            ]
            rows[1]["effect_cosine"] = None
        if mutate is not None:
            mutate(table_name, rows)
        path = root / f"{table_name}.csv"
        pd.DataFrame(rows).to_csv(path, index=False)
        input_tables[table_name] = str(path)
    return input_tables


def _config(input_tables: dict[str, str], output_root: Path) -> dict[str, object]:
    return {
        "analysis": "unit_test_nature_virtual_cell_evidence",
        "seed": 1,
        "bootstrap_reps": 10,
        "output_dir": str(output_root),
        "input_tables": input_tables,
        "primary_models": {
            "transition": "context_residual_mlp",
            "external_perturbation": "gene_context_ridge_residual",
        },
        "baseline_models": {
            "transition": ["identity", "mean_shift", "ridge", "mlp"],
            "external_perturbation": [
                "identity",
                "global_mean_shift",
                "gene_context_shift",
                "gene_context_mlp_residual",
            ],
        },
        "competence_window": {"label": "TS14-TS19", "stage_min": 14, "stage_max": 19},
        "schema": {
            "minimum_rows": {name: 1 for name in pipeline.REQUIRED_COLUMNS},
            "nullable_numeric_columns": {"external_condition_metrics": ["effect_cosine"]},
        },
        "claim_thresholds": {},
        "failure_policy": {
            "fail_on_missing_inputs": True,
            "fail_on_schema_error": True,
            "fail_on_missing_required_models": True,
            "fail_on_nonfinite_required_metrics": True,
            "fail_on_claim_threshold_failure": False,
        },
    }


def _paths(root: Path) -> pipeline.PipelinePaths:
    return pipeline.PipelinePaths(
        root=root,
        tables=root / "tables",
        figures=root / "figures",
        methods=root / "paper_methods_evidence_pipeline.md",
        manifest=root / "evidence_manifest.json",
    )


class NatureVirtualCellEvidencePipelineTest(unittest.TestCase):
    def test_schema_accepts_explicit_nullable_numeric_column(self) -> None:
        with tempfile.TemporaryDirectory(dir=PROJECT_ROOT / "results") as tmp_dir:
            root = Path(tmp_dir)
            paths = _paths(root / "outputs")
            paths.tables.mkdir(parents=True)
            paths.figures.mkdir(parents=True)

            config = _config(_write_required_tables(root), paths.root)
            validation, tables = pipeline.validate_inputs(config, paths)
            model_validation = pipeline.validate_required_models(config, tables, paths)

            self.assertTrue(validation["schema_pass"].all())
            self.assertTrue(model_validation["present_in_main_metrics"].fillna(True).all())

    def test_schema_rejects_undeclared_missing_required_numeric_value(self) -> None:
        def mutate(table_name: str, rows: list[dict[str, object]]) -> None:
            if table_name == "stimulus_response":
                rows[0]["stimulus_response_norm"] = None

        with tempfile.TemporaryDirectory(dir=PROJECT_ROOT / "results") as tmp_dir:
            root = Path(tmp_dir)
            paths = _paths(root / "outputs")
            paths.tables.mkdir(parents=True)
            paths.figures.mkdir(parents=True)

            config = _config(_write_required_tables(root, mutate=mutate), paths.root)
            with self.assertRaises(pipeline.EvidencePipelineError):
                pipeline.validate_inputs(config, paths)


if __name__ == "__main__":
    unittest.main()
