from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.discovery.common import cell_feature_frame, configure_paths, linear_effects, load_teacher, output_dirs, write_report


def run(config: str = "configs/train.yaml", quick_fixture: bool = False):
    train_cfg, model_cfg = configure_paths(config, quick_fixture)
    table_dir, report_dir, fig_dir = output_dirs(train_cfg)
    adata = load_teacher(model_cfg)
    frame = cell_feature_frame(adata, model_cfg)
    predictors = ["local_density", "ot_growth", "cell_cycle_score", "fate_entropy", "cci_signal"]
    stats = linear_effects(frame, "net_growth_hazard", predictors)
    events_path = Path(model_cfg.get("event_log_path", "tables/birth_death_event_log.csv"))
    events = pd.read_csv(events_path) if events_path.exists() else pd.DataFrame()
    density_bins = pd.qcut(frame["local_density"].rank(method="first"), q=min(8, max(2, frame.shape[0] // 20)), duplicates="drop")
    curve = frame.groupby(density_bins, observed=False).agg(
        local_density=("local_density", "mean"),
        birth_hazard=("birth_hazard", "mean"),
        death_hazard=("death_hazard", "mean"),
        net_growth_hazard=("net_growth_hazard", "mean"),
    )
    curve["event_birth_count"] = int((events.get("event", pd.Series(dtype=str)) == "birth").sum()) if not events.empty else 0
    curve["event_death_count"] = int((events.get("event", pd.Series(dtype=str)) == "death").sum()) if not events.empty else 0
    curve.reset_index(drop=True).to_csv(table_dir / "birth_death_law.csv", index=False)
    fig, ax = plt.subplots(figsize=(6, 4), dpi=160)
    ax.plot(curve["local_density"], curve["birth_hazard"], marker="o", label="birth")
    ax.plot(curve["local_density"], curve["death_hazard"], marker="o", label="death")
    ax.plot(curve["local_density"], curve["net_growth_hazard"], marker="o", label="net")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("local density")
    ax.set_ylabel("hazard")
    ax.set_title("Birth/death density phase curve")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(fig_dir / "birth_density_phase_curve.png")
    plt.close(fig)
    density_coef = stats.loc[stats["predictor"] == "local_density", "coef"]
    threshold_like = bool((curve["net_growth_hazard"].min() < 0) and (curve["net_growth_hazard"].max() > 0))
    gate = bool((not density_coef.empty and abs(float(density_coef.iloc[0])) > 1e-4) or threshold_like)
    write_report(
        report_dir / "discovery_birth_death_law.md",
        "Discovery Birth Death Law",
        [
            "Net growth hazard is analysed against local density, OT mass expansion, cell-cycle score, fate entropy and CCI signal.",
            "",
            "## Regression",
            "",
            stats.to_markdown(index=False) if not stats.empty else "No stable regression fit.",
            "",
            f"- carrying_capacity_like_threshold: {threshold_like}",
            f"- birth_death_law_gate: {gate}",
            f"- event counts: births={int((events.get('event', pd.Series(dtype=str)) == 'birth').sum()) if not events.empty else 0}, deaths={int((events.get('event', pd.Series(dtype=str)) == 'death').sum()) if not events.empty else 0}",
        ],
    )
    return {"law": "birth_death", "gate": gate, "table": str(table_dir / "birth_death_law.csv")}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/train.yaml")
    parser.add_argument("--quick-fixture", action="store_true")
    args = parser.parse_args()
    print(run(args.config, args.quick_fixture))


if __name__ == "__main__":
    main()

