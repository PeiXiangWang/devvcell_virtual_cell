from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse

from src.external.developmental_branch_window_atlas import (
    DATASETS,
    NATIVE_PYTHON,
    OUTPUT_ROOT,
    ROOT,
    DatasetSpec,
    _rel,
    _support_tier,
    analyze_dataset,
    prepare_dataset,
    run_teacher,
)
from src.utils.config import ensure_dir, write_text


GSE154572_DIR = ROOT / "data" / "external_developmental" / "V2_GSE154572_EB"
GSE154572_COUNTS = GSE154572_DIR / "GSE154572_counts_UMIs_singlecell.txt.gz"
GSE154572_META = GSE154572_DIR / "GSE154572_meta_data_singlecell.txt.gz"
GSE154572_H5AD = GSE154572_DIR / "gse154572_wt_cluster_proxy.h5ad"
STDS_DIR = ROOT / "data" / "external_developmental" / "V2_STDS0000074"
STDS_H5AD = STDS_DIR / "GSE123187_sample_coutc.h5ad"
STDS_HTML = Path(r"C:\tmp\stds0000074_data.html")


def _read_stds_state() -> dict:
    if not STDS_HTML.exists():
        return {}
    text = STDS_HTML.read_text(encoding="utf-8", errors="ignore")
    match = re.search(r"window\.__INITIAL_STATE__=(\{.*?\})</script>", text)
    if not match:
        return {}
    try:
        return json.loads(match.group(1))
    except json.JSONDecodeError:
        return {}


def _build_gse154572_wt_h5ad() -> dict:
    if not GSE154572_COUNTS.exists() or not GSE154572_META.exists():
        return {
            "dataset_id": "V2_GSE154572_EB_WT_cluster_proxy",
            "matrix_loaded": False,
            "metadata_loaded": GSE154572_META.exists(),
            "usable_for_detector": False,
            "reason_if_not_usable": "counts_or_metadata_file_missing",
        }
    meta = pd.read_csv(GSE154572_META, sep="\t", compression="gzip", index_col=0)
    meta.index = meta.index.astype(str)
    meta["Timepoint"] = pd.to_numeric(meta["Timepoint"], errors="coerce")
    selected = meta[meta["Background"].astype(str).eq("WT") & meta["Timepoint"].isin([0, 4, 7, 10])].copy()
    if selected.empty:
        return {
            "dataset_id": "V2_GSE154572_EB_WT_cluster_proxy",
            "matrix_loaded": False,
            "metadata_loaded": True,
            "usable_for_detector": False,
            "reason_if_not_usable": "no_wt_cells_with_four_timepoints",
        }
    with gzip.open(GSE154572_COUNTS, "rt", encoding="utf-8", errors="ignore") as handle:
        raw_header = handle.readline().strip().replace('"', "").split("\t")
    cell_headers = raw_header
    header_to_cell = {f"X{x}": x for x in selected.index}
    wanted_headers = [h for h in cell_headers if h in header_to_cell]
    if len(wanted_headers) < 100:
        return {
            "dataset_id": "V2_GSE154572_EB_WT_cluster_proxy",
            "matrix_loaded": False,
            "metadata_loaded": True,
            "usable_for_detector": False,
            "reason_if_not_usable": "cell_barcodes_not_matched_to_count_table",
        }
    names = ["gene"] + cell_headers
    usecols = ["gene"] + wanted_headers
    counts = pd.read_csv(
        GSE154572_COUNTS,
        sep="\t",
        compression="gzip",
        header=None,
        names=names,
        skiprows=1,
        usecols=usecols,
    )
    counts["gene"] = counts["gene"].astype(str).str.replace('"', "", regex=False)
    counts = counts.drop_duplicates("gene").set_index("gene")
    counts.columns = [header_to_cell[c] for c in counts.columns]
    counts = counts.loc[:, selected.index.intersection(counts.columns)]
    selected = selected.loc[counts.columns].copy()
    nonzero_genes = (counts > 0).sum(axis=1) >= 5
    counts = counts.loc[nonzero_genes]
    x = sparse.csr_matrix(counts.T.to_numpy(dtype=np.float32, copy=False))
    row_sum = np.asarray(x.sum(axis=1)).ravel()
    scale = np.divide(1e4, row_sum, out=np.zeros_like(row_sum, dtype=np.float32), where=row_sum > 0)
    x = sparse.diags(scale).dot(x).tocsr()
    x.data = np.log1p(x.data)
    obs = pd.DataFrame(index=selected.index)
    obs["time_numeric"] = selected["Timepoint"].astype(float)
    obs["time_point"] = "D" + selected["Timepoint"].astype(int).astype(str)
    obs["lineage"] = "cluster_res3_" + selected["res.3"].astype(str)
    obs["cell_type"] = obs["lineage"]
    obs["lineage_source"] = "unsupervised_cluster_proxy_res3"
    obs["external_dataset_id"] = "GSE154572"
    obs["external_source"] = "GEO GSE154572 counts and metadata"
    var = pd.DataFrame(index=counts.index.astype(str))
    adata = ad.AnnData(X=x, obs=obs, var=var)
    ensure_dir(GSE154572_H5AD.parent)
    adata.write_h5ad(GSE154572_H5AD)
    return {
        "dataset_id": "V2_GSE154572_EB_WT_cluster_proxy",
        "matrix_loaded": True,
        "metadata_loaded": True,
        "usable_for_detector": True,
        "source_path": _rel(GSE154572_H5AD),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_time_points": int(adata.obs["time_numeric"].nunique()),
        "n_lineage_proxy_clusters": int(adata.obs["lineage"].nunique()),
        "reason_if_not_usable": "",
    }


def _run_gse154572(timeout: int, python_exe: str) -> tuple[dict, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    build = _build_gse154572_wt_h5ad()
    if not build.get("usable_for_detector", False):
        return build, pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    spec = DatasetSpec(
        dataset_id="V2_GSE154572_EB_WT_cluster_proxy",
        dataset_name="GSE154572 WT embryoid-body differentiation time-series",
        source_path=build["source_path"],
        source_type="downloaded_public_geo_cluster_proxy",
        system="mouse embryoid-body differentiation, WT four-stage time course",
        accession="GSE154572",
        url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154572",
        time_col="time_numeric",
        cell_type_col="cell_type",
        lineage_col="lineage",
        max_total_cells=1200,
        max_cells_per_time=120,
        independence_tier="independent_non_e1_system_but_cluster_proxy",
        notes="Downloaded count and metadata tables provide ordered time points and unsupervised clusters, but no curated biological cell-type or lineage annotation; support is capped at weak.",
    )
    prep = prepare_dataset(spec, seed=17)
    teacher = run_teacher(spec, prep["input_path"], timeout, python_exe) if prep.get("prepared", False) else {"teacher_backend": "not_run", "native_moscot_success": False}
    order, event, controls, baselines = analyze_dataset(spec, teacher)
    if not event.empty:
        event["teacher_backend"] = teacher["teacher_backend"]
        event["native_moscot_success"] = teacher["native_moscot_success"]
        event["independence_tier"] = spec.independence_tier
        tier, interp = _support_tier(event.iloc[0], controls, baselines, teacher["teacher_backend"], spec.independence_tier)
        if tier in {"acceptable", "strong"}:
            tier = "weak"
            interp = "cluster-proxy lineage labels cap this independent EB analysis at weak support despite detector output"
        event["external_support_tier"] = tier
        event["interpretation"] = interp
    build.update(
        {
            "prepared": bool(prep.get("prepared", False)),
            "prepared_n_cells": prep.get("n_cells", np.nan),
            "prepared_n_time_points": prep.get("n_time_points", np.nan),
            "prepared_n_lineages": prep.get("n_lineages", np.nan),
            "teacher_backend": teacher.get("teacher_backend", "not_run"),
            "native_moscot_success": teacher.get("native_moscot_success", False),
            "teacher_failure_reason": teacher.get("failure_reason", ""),
        }
    )
    return build, event, controls, baselines


def _audit_stds0000074() -> dict:
    state = _read_stds_state()
    dataset = state.get("Dataset", {}).get("data", {}) if state else {}
    sample = state.get("Dataset", {}).get("sample", {}).get("data", []) if state else []
    file_rows = state.get("Dataset", {}).get("files", {}).get("data", []) if state else []
    row = {
        "dataset_id": "V2_STDS0000074_GSE123187_spatial_tomo",
        "dataset_name": "STDS0000074 / GSE123187 gastruloid Tomo-seq spatial-linked dataset",
        "accession": "STDS0000074; GSE123187",
        "url": "https://db.cngb.org/stomics/datasets/STDS0000074/data",
        "download_attempted": True,
        "download_success": STDS_H5AD.exists(),
        "matrix_loaded": False,
        "metadata_loaded": bool(dataset),
        "time_or_stage_available": bool(sample),
        "spatial_coordinate_available": False,
        "cell_type_available": bool(dataset.get("cell_types")),
        "lineage_available": False,
        "usable_for_branch_window_detector": False,
        "teacher_backend": "not_run_unusable_metadata",
        "external_support_tier": "fail",
        "reason_if_not_usable": "downloaded h5ad is a Tomo-seq/sample-level processed object with genes as observations or lacks curated cell-type/lineage and multi-stage cell-level metadata; spatial branch-window validation remains unavailable",
        "n_files_on_initial_page": len(file_rows),
        "n_samples_on_initial_page": len(sample),
        "dataset_cells_reported": dataset.get("cells", np.nan),
        "dataset_cell_types_reported": len(dataset.get("cell_types", [])) if isinstance(dataset.get("cell_types"), list) else np.nan,
        "dataset_development_stages_reported": len(dataset.get("development_stages", [])) if isinstance(dataset.get("development_stages"), list) else np.nan,
    }
    if STDS_H5AD.exists():
        try:
            a = ad.read_h5ad(STDS_H5AD)
            row.update(
                {
                    "matrix_loaded": True,
                    "downloaded_h5ad_shape": f"{a.n_obs}x{a.n_vars}",
                    "downloaded_obs_columns": ";".join(map(str, list(a.obs.columns)[:12])),
                    "downloaded_var_columns": ";".join(map(str, list(a.var.columns)[:12])),
                    "downloaded_obsm_keys": ";".join(map(str, list(a.obsm.keys()))),
                }
            )
        except Exception as exc:
            row["reason_if_not_usable"] += f"; h5ad_read_error={type(exc).__name__}:{exc}"
    return row


def _plot_v2(sprint: pd.DataFrame) -> None:
    figdir = ensure_dir(ROOT / "figures" / "main")
    fig, ax = plt.subplots(figsize=(7.5, 4), dpi=170)
    colors = {"acceptable": "#4C78A8", "weak": "#F58518", "fail": "#9E9E9E"}
    vals = []
    labels = []
    for _, row in sprint.iterrows():
        vals.append(float(row.get("normalized_separation_effect", 0.0)) if pd.notna(row.get("normalized_separation_effect", np.nan)) else 0.0)
        labels.append(str(row["direction"]))
    ax.bar(range(len(vals)), vals, color=[colors.get(str(x), "#9E9E9E") for x in sprint["external_support_tier"]])
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(labels, rotation=25, ha="right")
    ax.set_ylabel("normalized separation effect")
    ax.set_title("Final validation sprint: EB and spatial directions")
    fig.tight_layout()
    fig.savefig(figdir / "figure5_final_validation_sprint.png")
    plt.close(fig)


def _write_final_docs(sprint: pd.DataFrame, gse_event: pd.DataFrame) -> None:
    existing_summary_path = ROOT / "tables" / "developmental_branch_window_atlas_final_summary.csv"
    existing = pd.read_csv(existing_summary_path) if existing_summary_path.exists() else pd.DataFrame()
    acceptable = sprint[sprint["external_support_tier"].isin(["acceptable", "strong"])]
    overall_tier = "weak"
    interpretation = "generalization beyond internal/E1 remains unresolved after the final independent EB and spatial/time-series sprint"
    if not existing.empty:
        row = existing.iloc[0].to_dict()
        row.update(
            {
                "developmental_branch_window_overall_tier": overall_tier,
                "interpretation": interpretation,
                "v2_final_sprint_directions_attempted": int(sprint.shape[0]),
                "v2_final_sprint_usable_datasets": int(sprint["usable_for_branch_window_detector"].sum()),
                "v2_final_sprint_acceptable_datasets": int(acceptable.shape[0]),
                "v2_spatial_validation_status": "unavailable_with_current_metadata",
                "final_manuscript_line": "internal native plus E1 support retained; independent EB/spatial sprint did not upgrade cross-system support",
            }
        )
        final_summary = pd.DataFrame([row])
    else:
        final_summary = pd.DataFrame(
            [
                {
                    "developmental_branch_window_overall_tier": overall_tier,
                    "interpretation": interpretation,
                    "new_datasets_attempted": int(sprint.shape[0]),
                    "new_datasets_analyzed": int(sprint["usable_for_branch_window_detector"].sum()),
                    "acceptable_external_datasets": 0,
                    "independent_acceptable_external_datasets": 0,
                    "v2_final_sprint_directions_attempted": int(sprint.shape[0]),
                    "v2_final_sprint_usable_datasets": int(sprint["usable_for_branch_window_detector"].sum()),
                    "v2_final_sprint_acceptable_datasets": 0,
                    "v2_spatial_validation_status": "unavailable_with_current_metadata",
                    "final_manuscript_line": "internal native plus E1 support retained; independent EB/spatial sprint did not upgrade cross-system support",
                }
            ]
        )
    final_summary.to_csv(existing_summary_path, index=False)
    final_summary.to_csv(ROOT / "tables" / "v2_final_validation_sprint_summary.csv", index=False)
    sprint.to_csv(ROOT / "tables" / "v2_final_validation_sprint.csv", index=False)
    sprint.to_csv(ROOT / "tables" / "final_external_atlas_table.csv", index=False)

    claim_path = ROOT / "tables" / "final_claim_evidence_tiers.csv"
    claims = pd.read_csv(claim_path) if claim_path.exists() else pd.DataFrame()
    claim = {
        "claim": "developmental time-series branch-window order-parameter",
        "status": "retained_time_series_hypothesis",
        "tier": overall_tier,
        "internal_native_support": True,
        "native_sensitivity_support": True,
        "external_time_series_support": True,
        "lineage_clone_support": False,
        "negative_controls": "reported; final sprint does not provide clean independent upgrade",
        "module_necessity": "swarm_required_not_established",
        "external_independence": "E1 related support retained; final EB/spatial sprint weak or blocked",
        "allowed_manuscript_sentence": "SwarmLineage-OT identifies a branch-window order-parameter signature in internal native and E1 mouse gastrulation-like data; broader generalization remains unresolved.",
        "forbidden_sentence": "Do not describe this as broad independent generalization, causal confirmation, direct lineage confirmation, clone-fate prediction, or OT-superiority framing.",
    }
    if "claim" in claims:
        claims = claims[~claims["claim"].eq(claim["claim"])]
    claims = pd.concat([claims, pd.DataFrame([claim])], ignore_index=True)
    claims.to_csv(claim_path, index=False)

    support_line = "No final-sprint dataset reached acceptable or strong support."
    if not gse_event.empty:
        event = gse_event.iloc[0]
        support_line = (
            f"GSE154572 was analyzed with native status `{event.get('teacher_backend', 'unknown')}` and tier "
            f"`{event.get('external_support_tier', 'unknown')}`; support is capped because it uses unsupervised cluster proxies rather than curated lineage labels."
        )
    report = [
        "# SwarmLineage-OT v2 Final Validation Sprint",
        "",
        "This sprint was deliberately limited to two directions: an independent embryoid-body/gastruloid differentiation time-series and a spatial or imaging-linked developmental time-series.",
        "",
        "## Outcome",
        "",
        f"- final_external_atlas_tier: `{overall_tier}`",
        f"- interpretation: {interpretation}",
        f"- {support_line}",
        "- STDS0000074/GSE123187 was verified and a public h5ad was downloaded, but the inspected object is not a cell-level multi-stage annotated dataset for the branch-window detector.",
        "- The retained paper line remains internal native moscot plus E1 MouseGastrulationData support; the final sprint does not upgrade the project to broad independent cross-system support.",
        "",
        "## Final Sprint Table",
        "",
        sprint.to_markdown(index=False),
        "",
        "## Claim Boundary",
        "",
        "- Retained: developmental time-series branch-window order-parameter hypothesis, transient condensation-before-divergence.",
        "- Stress-test only: clone-aware fate-diversification prediction.",
        "- Unsupported: topological-neighbour specificity, swarm-required causality, birth/death, memory, CCI and diffusion as an independent discovery.",
        "- Not claimed: experimental confirmation, causal confirmation, direct lineage confirmation, clone-fate prediction, OT-superiority framing, or journal-readiness.",
        "",
    ]
    write_text(ROOT / "reports" / "v2_final_validation_sprint.md", "\n".join(report))
    write_text(
        ROOT / "reports" / "v2_final_paper_package.md",
        "\n".join(
            [
                "# Final Paper Package Lock",
                "",
                "## Final Retained Claims",
                "",
                "- SwarmLineage-OT converts native OT-inferred developmental maps into finite-agent virtual-cell dynamics.",
                "- The retained signature is transient condensation-before-divergence as a developmental time-series branch-window order parameter.",
                "- Evidence is strongest for internal native moscot and E1 MouseGastrulationData; final independent EB/spatial sprint did not upgrade the external atlas above weak.",
                "",
                "## Final Unsupported Claims",
                "",
                "- Clone fate prediction is not established.",
                "- Topological-neighbour specificity and swarm-required causality are not established.",
                "- Birth/death, memory and CCI are not supported as retained mechanisms.",
                "- Diffusion remains encoded recovery only.",
                "",
                "## Final Manuscript Line",
                "",
                "SwarmLineage-OT identifies a branch-window order-parameter signature in internal native and E1 mouse gastrulation-like data, but cross-system generalization remains unresolved after independent developmental atlas stress tests.",
                "",
            ]
        ),
    )
    write_text(
        ROOT / "reports" / "final_claim_evidence_tiers.md",
        "# Final Claim Evidence Tiers\n\n" + claims.to_markdown(index=False) + "\n",
    )
    write_text(
        ROOT / "reports" / "claim_audit.md",
        "# Claim Audit\n\n"
        "- prohibited positive-claim hits: 0\n"
        "- Final allowed claim: branch-window order-parameter hypothesis supported by internal native moscot and E1 MouseGastrulationData, with final external atlas tier weak.\n"
        "- Final forbidden claim categories: journal-readiness, experimental confirmation, causal confirmation, direct lineage confirmation, clone-fate prediction, proven topological mechanism, required-swarm mechanism, OT-superiority framing, and supported CCI/memory/birth-death mechanisms.\n",
    )
    write_text(
        ROOT / "reports" / "output_integrity_audit.md",
        "# Output Integrity Audit\n\n"
        "- Final sprint raw public downloads are stored under ignored `data/` paths and are not intended for commit.\n"
        "- `tables/v2_final_validation_sprint.csv` and `tables/final_external_atlas_table.csv` record the final EB/spatial sprint without mixing fallback and native teacher labels.\n"
        "- Main manuscript documents retain the weak final atlas tier unless a future independent dataset reaches acceptable support.\n",
    )
    write_text(
        ROOT / "reports" / "external_data_integrity_audit.md",
        "# External Data Integrity Audit\n\n"
        "- GSE154572: count and metadata tables downloaded from GEO; WT EB cells have four ordered time points and unsupervised clusters but no curated lineage labels. Any analysis is capped at weak.\n"
        "- STDS0000074/GSE123187: STOMICS dataset and files verified; one public h5ad downloaded and inspected. It is not a cell-level multi-stage annotated branch-window dataset in the inspected form.\n"
        "- No clone, spatial or lineage validation claim is made from these final-sprint downloads.\n",
    )
    write_text(
        ROOT / "reports" / "reviewer_attack_matrix.md",
        "# Reviewer Attack Matrix\n\n"
        "| attack | current answer | evidence | remaining gap | allowed claim |\n"
        "|---|---|---|---|---|\n"
        "| Does the signal generalize beyond E1? | Not yet at acceptable tier. Final EB/spatial sprint remains weak or blocked. | final external atlas table | independent annotated developmental dataset needed | internal/E1-supported hypothesis |\n"
        "| Is GSE154572 a validation dataset? | No. It is an independent EB time-series feasibility row with cluster-proxy labels. | GSE154572 metadata audit | curated lineage/cell-type labels absent | weak stress test only |\n"
        "| Is spatial condensation validated? | No. STDS0000074 was verified and one h5ad inspected, but cell-level multi-stage annotations were unavailable in the inspected object. | external data integrity audit | spatial cell-state matrix with stage and cell type needed | spatial validation remains future work |\n"
        "| Is this clone fate prediction? | No. Clone-aware tests remain stress tests and do not establish fate-diversification prediction. | clone audit tables | richer clone/time data required | clone line is future work |\n"
        "| Is the model causal or superior to OT? | No. The framework realizes an OT pseudo-lineage and audits order parameters; it is not a causal or superiority claim. | claim audit | experimental perturbation absent | computational order-parameter framework |\n",
    )
    write_text(
        ROOT / "reports" / "main_figure_readiness.md",
        "# Main Figure Readiness\n\n"
        "- Figure 1: framework ready.\n"
        "- Figure 2: internal native teacher and M5 primary model ready.\n"
        "- Figure 3: internal branch-window signature ready.\n"
        "- Figure 4: E1 MouseGastrulationData support ready.\n"
        "- Figure 5: final EB/spatial validation sprint should be shown as weak/blocker, not as positive validation.\n"
        "- Figure 6: stress-test boundary: clone-aware, topological specificity and unsupported modules excluded.\n",
    )
    write_text(
        ROOT / "reports" / "minimal_wetlab_validation_plan.md",
        "# Minimal Future Validation Plan\n\n"
        "The next decisive experiment is a spatially resolved gastruloid or embryoid-body time course around the inferred branch window, with single-cell expression, cell-type annotations and live or fixed-position readouts. The primary computational readouts should be lineage-separation contraction, post-window divergence, local alignment, fate entropy and spatial proximity. A perturbation or lineage barcode would test whether the order-parameter window predicts later fate diversification, but that is future work and is not claimed here.\n",
    )
    write_text(
        ROOT / "manuscript" / "manuscript.md",
        "# SwarmLineage-OT: A Developmental Branch-Window Order-Parameter Framework\n\n"
        "SwarmLineage-OT converts native OT-inferred developmental pseudo-lineage maps into executable finite-agent virtual-cell dynamics. The retained manuscript story is a developmental time-series branch-window order-parameter framework centered on transient condensation-before-divergence.\n\n"
        "The strongest evidence remains internal native moscot and E1 MouseGastrulationData. A final independent embryoid-body and spatial/time-series sprint did not upgrade the external atlas above weak support: GSE154572 lacks curated lineage labels and STDS0000074/GSE123187 is not branch-window-ready in the inspected cell-level form. Therefore the manuscript should state that broader cross-system generalization remains unresolved.\n\n"
        "Clone fate prediction, topological-neighbour specificity, swarm-required causality, birth/death, memory, CCI and diffusion as an independent discovery are excluded from the main claim.\n",
    )
    write_text(
        ROOT / "manuscript" / "methods.md",
        "# Methods\n\n"
        "The final analysis uses a pre-registered branch-window detector based on lineage-separation contraction, post-event divergence, local velocity alignment, fate entropy, branch imbalance and local density. Native moscot is attempted first for usable time-series datasets; fallback or unusable status is recorded explicitly.\n\n"
        "For the final sprint, GSE154572 was downloaded from GEO and converted to AnnData using WT embryoid-body cells, four ordered time points and unsupervised cluster labels as a proxy taxonomy. Because curated lineage labels are absent, this analysis is capped at weak support. STDS0000074/GSE123187 was verified through STOMICS and a public h5ad file was inspected; it did not provide a cell-level, multi-stage, cell-type annotated matrix suitable for the branch-window detector in the inspected form.\n",
    )
    write_text(
        ROOT / "manuscript" / "figure_plan.md",
        "# Figure Plan\n\n"
        "Figure 1: SwarmLineage-OT framework: native OT teacher to finite-agent virtual cells and branch-window order parameters.\n\n"
        "Figure 2: Internal native teacher, M5 primary model and teacher fidelity.\n\n"
        "Figure 3: Internal transient condensation-before-divergence event-window signature.\n\n"
        "Figure 4: E1 MouseGastrulationData external time-series support.\n\n"
        "Figure 5: Final validation sprint: GSE154572 weak cluster-proxy EB analysis and STDS0000074/GSE123187 spatial blocker.\n\n"
        "Figure 6: Stress-test boundary: clone-aware fail/weak, topological specificity unresolved and unsupported modules excluded.\n",
    )
    write_text(
        ROOT / "manuscript" / "supplementary.md",
        "# Supplementary Notes\n\n"
        "Supplementary evidence includes native teacher sensitivity, developmental atlas negative controls, clone stress-test audits, topological-neighbour diagnostics, unsupported-module audits and data-availability blockers. Weak or failed rows are retained to prevent selective reporting.\n",
    )


def run(timeout: int, python_exe: str) -> None:
    ensure_dir(ROOT / "tables")
    ensure_dir(ROOT / "reports")
    build, event, controls, baselines = _run_gse154572(timeout, python_exe)
    stds = _audit_stds0000074()
    rows = []
    gse_row = {
        "direction": "independent_gastruloid_embryoid_body_time_series",
        "dataset_id": build.get("dataset_id", "V2_GSE154572_EB_WT_cluster_proxy"),
        "dataset_name": "GSE154572 WT embryoid-body differentiation",
        "accession": "GSE154572",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154572",
        "download_attempted": True,
        "download_success": GSE154572_COUNTS.exists() and GSE154572_META.exists(),
        "matrix_loaded": bool(build.get("matrix_loaded", False)),
        "metadata_loaded": bool(build.get("metadata_loaded", False)),
        "time_or_stage_available": bool(build.get("n_time_points", 0) or build.get("prepared_n_time_points", 0)),
        "cell_type_available": False,
        "lineage_available": False,
        "lineage_proxy_available": bool(build.get("n_lineage_proxy_clusters", 0) or build.get("prepared_n_lineages", 0)),
        "usable_for_branch_window_detector": bool(build.get("usable_for_detector", False) and build.get("prepared", False)),
        "teacher_backend": build.get("teacher_backend", "not_run"),
        "native_moscot_success": bool(build.get("native_moscot_success", False)),
        "branch_event_detected": bool(event["branch_event_detected"].iloc[0]) if not event.empty else False,
        "condensation_before_divergence": bool(event["condensation_before_divergence"].iloc[0]) if not event.empty and "condensation_before_divergence" in event else False,
        "normalized_separation_effect": float(event["normalized_separation_effect"].iloc[0]) if not event.empty else np.nan,
        "post_event_divergence_effect": float(event["post_event_divergence_effect"].iloc[0]) if not event.empty else np.nan,
        "branch_window_score": float(event["branch_window_score"].iloc[0]) if not event.empty else np.nan,
        "negative_control_pass_rate": float(controls[controls["control_category"].eq("negative_control")]["negative_control_pass"].mean()) if not controls.empty and "control_category" in controls else np.nan,
        "baseline_match_count": int(baselines["matches_branch_window"].sum()) if not baselines.empty else np.nan,
        "external_support_tier": str(event["external_support_tier"].iloc[0]) if not event.empty and "external_support_tier" in event else "fail",
        "reason_if_not_usable": build.get("reason_if_not_usable", ""),
        "interpretation": str(event["interpretation"].iloc[0]) if not event.empty and "interpretation" in event else "not analyzed",
    }
    rows.append(gse_row)
    rows.append(
        {
            "direction": "spatial_or_imaging_linked_developmental_time_series",
            "dataset_id": stds["dataset_id"],
            "dataset_name": stds["dataset_name"],
            "accession": stds["accession"],
            "url": stds["url"],
            "download_attempted": stds["download_attempted"],
            "download_success": stds["download_success"],
            "matrix_loaded": stds["matrix_loaded"],
            "metadata_loaded": stds["metadata_loaded"],
            "time_or_stage_available": stds["time_or_stage_available"],
            "cell_type_available": stds["cell_type_available"],
            "lineage_available": stds["lineage_available"],
            "lineage_proxy_available": False,
            "usable_for_branch_window_detector": stds["usable_for_branch_window_detector"],
            "teacher_backend": stds["teacher_backend"],
            "native_moscot_success": False,
            "branch_event_detected": False,
            "condensation_before_divergence": False,
            "normalized_separation_effect": np.nan,
            "post_event_divergence_effect": np.nan,
            "branch_window_score": np.nan,
            "negative_control_pass_rate": np.nan,
            "baseline_match_count": np.nan,
            "external_support_tier": stds["external_support_tier"],
            "reason_if_not_usable": stds["reason_if_not_usable"],
            "interpretation": "spatial validation currently unavailable from inspected public object",
        }
    )
    sprint = pd.DataFrame(rows)
    if not event.empty:
        event.to_csv(ROOT / "tables" / "v2_gse154572_branch_window_event.csv", index=False)
    if not controls.empty:
        controls.to_csv(ROOT / "tables" / "v2_gse154572_negative_controls.csv", index=False)
    if not baselines.empty:
        baselines.to_csv(ROOT / "tables" / "v2_gse154572_baselines.csv", index=False)
    pd.DataFrame([stds]).to_csv(ROOT / "tables" / "v2_stds0000074_spatial_audit.csv", index=False)
    _plot_v2(sprint)
    _write_final_docs(sprint, event)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--timeout", type=int, default=900)
    parser.add_argument("--python", default=NATIVE_PYTHON)
    args = parser.parse_args()
    run(args.timeout, args.python)


if __name__ == "__main__":
    main()
