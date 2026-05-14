from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st


ROOT = Path(__file__).resolve().parents[2]


def read_csv(path: str) -> pd.DataFrame:
    file_path = ROOT / path
    if not file_path.exists():
        return pd.DataFrame()
    return pd.read_csv(file_path)


def read_text(path: str) -> str:
    file_path = ROOT / path
    if not file_path.exists():
        return f"Missing: `{path}`"
    return file_path.read_text(encoding="utf-8", errors="ignore")


def metric_from(frame: pd.DataFrame, column: str, default: str = "missing") -> str:
    if frame.empty or column not in frame:
        return default
    value = frame[column].iloc[0]
    return str(value)


def show_table(title: str, frame: pd.DataFrame, height: int = 280) -> None:
    st.subheader(title)
    if frame.empty:
        st.info("Table not found or empty.")
        return
    st.dataframe(frame, width="stretch", height=height)


def show_png(path: str, caption: str) -> None:
    file_path = ROOT / path
    if file_path.exists():
        st.image(str(file_path), caption=caption, width="stretch")
    else:
        st.caption(f"Not generated yet: `{path}`")


def main() -> None:
    st.set_page_config(page_title="SwarmLineage-OT v0.6", layout="wide")
    st.title("SwarmLineage-OT v0.6")
    st.caption(
        "Native moscot teacher sensitivity, evidence-selected primary agent, "
        "and branch-nucleation discovery audit."
    )

    teacher = read_csv("tables/teacher_backend_status.csv")
    mech = read_csv("tables/mechanistic_usefulness_summary.csv")
    laws = read_csv("tables/emergent_law_gate_summary.csv")
    primary = read_csv("tables/primary_agent_selection.csv")
    sensitivity = read_csv("tables/native_teacher_sensitivity.csv")
    external = read_csv("tables/external_dataset_registry.csv")

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Teacher backend", metric_from(teacher, "teacher_backend"))
    c2.metric("Teacher fidelity", metric_from(mech, "teacher_fidelity_tier"))
    c3.metric("Emergent laws", metric_from(mech, "emergent_law_tier"))
    c4.metric("Mechanistic usefulness", metric_from(mech, "mechanistic_usefulness_tier"))

    if not primary.empty:
        selected = primary[primary["recommendation"] == "primary_mechanistic_model"]
        if not selected.empty:
            st.success(f"Primary mechanistic model: `{selected['model'].iloc[0]}`")

    st.warning(
        "This dashboard summarizes computational evidence only. It does not claim "
        "experimental validation, causal proof, or high-impact submission readiness."
    )

    tabs = st.tabs(
        [
            "Overview",
            "Native Teacher",
            "Branch Nucleation",
            "Primary Agent",
            "External Registry",
            "Reports",
        ]
    )

    with tabs[0]:
        show_table("Emergent-law gate summary", laws)
        show_png("figures/main/figure1_framework.png", "Framework draft")

    with tabs[1]:
        show_table("Native teacher sensitivity", sensitivity, height=360)
        left, right = st.columns(2)
        with left:
            show_png("figures/discovery/native_teacher_sensitivity.png", "Native teacher sensitivity")
        with right:
            show_png(
                "figures/discovery/native_teacher_barycentric_stability.png",
                "Barycentric stability",
            )

    with tabs[2]:
        show_table("Branch-nucleation model comparison", read_csv("tables/branch_nucleation_model_comparison.csv"))
        show_table("Branch-nucleation event windows", read_csv("tables/branch_nucleation_event_windows.csv"))
        left, right = st.columns(2)
        with left:
            show_png("figures/discovery/branch_nucleation_event_window.png", "Event-window signature")
        with right:
            show_png("figures/main/figure4_branch_nucleation.png", "Model comparison draft")

    with tabs[3]:
        show_table("Primary agent selection", primary)
        st.markdown(read_text("reports/primary_agent_selection.md"))

    with tabs[4]:
        show_table("External dataset registry", external, height=260)
        show_table(
            "External branch-nucleation validation status",
            read_csv("tables/external_branch_nucleation_validation.csv"),
            height=220,
        )
        show_png("figures/discovery/external_branch_nucleation.png", "External validation status")

    with tabs[5]:
        report_choice = st.selectbox(
            "Report",
            [
                "reports/native_teacher_sensitivity.md",
                "reports/branch_nucleation_mechanism_summary.md",
                "reports/discovery_branch_nucleation.md",
                "reports/scientific_gap_audit.md",
                "reports/editorial_assessment.md",
                "manuscript/final_retained_results_and_methods.md",
            ],
        )
        st.markdown(read_text(report_choice))


if __name__ == "__main__":
    main()
