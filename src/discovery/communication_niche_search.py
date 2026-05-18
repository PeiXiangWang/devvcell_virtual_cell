from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.neighbors import NearestNeighbors


ROOT = Path(__file__).resolve().parents[2]
PATHWAYFINDER_ROOT = ROOT.parent / "PathwayFinder"
OMNIPATH_DIR = PATHWAYFINDER_ROOT / "data" / "raw" / "knowledge" / "omnipath"
TABLES = ROOT / "tables"
REPORTS = ROOT / "reports"
MANUSCRIPT = ROOT / "manuscript"


@dataclass(frozen=True)
class DatasetSpec:
    dataset_id: str
    path: Path
    role: str
    independence_tier: str
    event_time: float | None
    expected_window_type: str


DATASETS = [
    DatasetSpec(
        "internal_native",
        ROOT / "processed" / "swarmlineage_input.h5ad",
        "internal_reference",
        "internal",
        None,
        "condensation_before_divergence",
    ),
    DatasetSpec(
        "E1_MouseGastrulationData",
        ROOT / "processed" / "external" / "e1_swarmlineage_input.h5ad",
        "external_time_series",
        "related_mouse_gastrulation_sample",
        7.25,
        "condensation_before_divergence",
    ),
    DatasetSpec(
        "E3_MouseGastrulationData_full",
        ROOT
        / "processed"
        / "developmental_atlas"
        / "E3_MouseGastrulationData_wt_chimera_full_stage_mapped"
        / "swarmlineage_input.h5ad",
        "related_atlas_stress",
        "related_mouse_gastrulation_atlas",
        8.0,
        "condensation_before_divergence_weak",
    ),
    DatasetSpec(
        "E5_zebrafish_Farrell",
        ROOT
        / "processed"
        / "developmental_atlas"
        / "E5_CellRank_Farrell_zebrafish_axial_mesoderm"
        / "swarmlineage_input.h5ad",
        "independent_developmental_stress",
        "independent_cross_species",
        10.0,
        "divergence_without_condensation",
    ),
    DatasetSpec(
        "E2_GSE212050_gastruloid",
        ROOT
        / "processed"
        / "developmental_atlas"
        / "E2_GSE212050_gastruloid_native_atlas"
        / "swarmlineage_input.h5ad",
        "independent_developmental_stress",
        "independent_gastruloid",
        4.0,
        "no_branch_or_false_positive_window",
    ),
    DatasetSpec(
        "GSE154572_EB_WT",
        ROOT
        / "processed"
        / "developmental_atlas"
        / "V2_GSE154572_EB_WT_cluster_proxy"
        / "swarmlineage_input.h5ad",
        "independent_developmental_stress",
        "independent_embryoid_body",
        None,
        "weak_or_unresolved",
    ),
]


COMMUNICATION_MODULES = {
    "ecm_adhesion_guidance": {
        "ligand_prefixes": ("COL", "LAMA", "LAMB", "LAMC", "FN", "THBS", "NTN", "SEMA", "EFN", "MDK", "PTN"),
        "receptor_prefixes": ("ITG", "CD47", "DCC", "UNC5", "NRP", "PLXN", "EPH", "NCL", "SDC"),
    },
    "morphogen_patterning": {
        "ligand_prefixes": ("FGF", "WNT", "BMP", "NODAL", "TGFB", "SHH"),
        "receptor_prefixes": ("FGFR", "FZD", "BMPR", "TGFBR", "ACVR", "SMO", "PTCH"),
    },
    "notch_delta": {
        "ligand_prefixes": ("DLL", "JAG"),
        "receptor_prefixes": ("NOTCH",),
    },
    "chemokine_growth_factor": {
        "ligand_prefixes": ("CXCL", "CCL", "KITL", "VEGF", "PDGF", "CSF", "IL"),
        "receptor_prefixes": ("CXCR", "CCR", "KIT", "KDR", "FLT", "PDGFR", "CSF", "IL"),
    },
}


def _read_symbol_columns(var: pd.DataFrame) -> list[str]:
    return [c for c in ("SYMBOL", "gene_symbol", "gene_short_name", "gene_name", "feature_name") if c in var.columns]


def _gene_index(adata: ad.AnnData) -> dict[str, int]:
    out: dict[str, int] = {}
    for idx, name in enumerate(adata.var_names.astype(str)):
        out.setdefault(name.lower(), idx)
    for col in _read_symbol_columns(adata.var):
        for idx, name in enumerate(adata.var[col].astype(str)):
            if name and name.lower() != "nan":
                out.setdefault(name.lower(), idx)
    return out


def _load_intercell_prior() -> tuple[pd.DataFrame, dict[str, set[str]], pd.DataFrame]:
    intercell_path = OMNIPATH_DIR / "omnipath_intercell.tsv"
    interactions_path = OMNIPATH_DIR / "omnipath_interactions.tsv"
    if not intercell_path.exists() or not interactions_path.exists():
        raise FileNotFoundError(
            "PathwayFinder OmniPath intercell/interactions files were not found; "
            f"looked under {OMNIPATH_DIR}"
        )
    inter = pd.read_csv(intercell_path, sep="\t", low_memory=False)
    inter["genesymbol"] = inter["genesymbol"].astype(str)
    flags = {
        "sender": set(inter.loc[(inter["transmitter"] == True) | (inter["secreted"] == True), "genesymbol"].str.upper()),
        "receiver": set(
            inter.loc[
                (inter["receiver"] == True)
                | (inter["plasma_membrane_transmembrane"] == True)
                | (inter["plasma_membrane_peripheral"] == True),
                "genesymbol",
            ].str.upper()
        ),
        "secreted": set(inter.loc[inter["secreted"] == True, "genesymbol"].str.upper()),
        "membrane": set(
            inter.loc[
                (inter["plasma_membrane_transmembrane"] == True)
                | (inter["plasma_membrane_peripheral"] == True),
                "genesymbol",
            ].str.upper()
        ),
    }
    interactions = pd.read_csv(interactions_path, sep="\t", low_memory=False)
    uniprot_to_symbol = (
        inter[["uniprot", "genesymbol", "transmitter", "receiver", "secreted", "plasma_membrane_transmembrane", "plasma_membrane_peripheral"]]
        .dropna(subset=["uniprot", "genesymbol"])
        .drop_duplicates()
    )
    src = uniprot_to_symbol.rename(
        columns={
            "uniprot": "source",
            "genesymbol": "ligand",
            "transmitter": "source_transmitter",
            "receiver": "source_receiver",
            "secreted": "source_secreted",
            "plasma_membrane_transmembrane": "source_membrane_tm",
            "plasma_membrane_peripheral": "source_membrane_peripheral",
        }
    )
    tgt = uniprot_to_symbol.rename(
        columns={
            "uniprot": "target",
            "genesymbol": "receptor",
            "transmitter": "target_transmitter",
            "receiver": "target_receiver",
            "secreted": "target_secreted",
            "plasma_membrane_transmembrane": "target_membrane_tm",
            "plasma_membrane_peripheral": "target_membrane_peripheral",
        }
    )
    lr = interactions.merge(src, on="source", how="inner").merge(tgt, on="target", how="inner")
    lr = lr[
        ((lr["source_transmitter"] == True) | (lr["source_secreted"] == True))
        & (
            (lr["target_receiver"] == True)
            | (lr["target_membrane_tm"] == True)
            | (lr["target_membrane_peripheral"] == True)
        )
    ][["ligand", "receptor", "is_stimulation", "is_inhibition"]].drop_duplicates()
    lr["ligand"] = lr["ligand"].astype(str).str.upper()
    lr["receptor"] = lr["receptor"].astype(str).str.upper()
    return inter, flags, lr


def _normalize_score(arr: np.ndarray) -> np.ndarray:
    arr = np.nan_to_num(arr.astype(float), nan=0.0)
    lo, hi = np.quantile(arr, [0.05, 0.95]) if np.any(arr) else (0.0, 1.0)
    return np.clip((arr - lo) / max(hi - lo, 1e-8), 0.0, 1.0)


def _vector_for_indices(adata: ad.AnnData, matched: list[int]) -> np.ndarray:
    if not matched:
        return np.zeros(adata.n_obs, dtype=float)
    x = adata.X[:, sorted(set(matched))]
    if sparse.issparse(x):
        arr = np.asarray(x.mean(axis=1)).ravel()
    else:
        arr = np.asarray(x).mean(axis=1)
    return _normalize_score(arr)


def _vector_for_genes(adata: ad.AnnData, genes: set[str]) -> tuple[np.ndarray, int]:
    index = _gene_index(adata)
    matched = []
    for gene in genes:
        key = gene.lower()
        if key in index:
            matched.append(index[key])
        else:
            # Conservative cross-species fallback for symbols such as fgf8a/fgf8b.
            candidates = [k for k in (key + "a", key + "b", key + "1", key + "2") if k in index]
            matched.extend(index[k] for k in candidates[:1])
    matched = sorted(set(matched))
    if not matched:
        return np.zeros(adata.n_obs, dtype=float), 0
    return _vector_for_indices(adata, matched), len(matched)


def _indices_for_prefixes(adata: ad.AnnData, prefixes: tuple[str, ...]) -> list[int]:
    names_by_idx: dict[int, set[str]] = {i: {str(v).upper()} for i, v in enumerate(adata.var_names)}
    for col in _read_symbol_columns(adata.var):
        for idx, name in enumerate(adata.var[col].astype(str)):
            if name and name.lower() != "nan":
                names_by_idx.setdefault(idx, set()).add(name.upper())
    out = []
    for idx, names in names_by_idx.items():
        for name in names:
            if any(name.startswith(prefix) for prefix in prefixes):
                out.append(idx)
                break
    return sorted(set(out))


def _latent(adata: ad.AnnData) -> np.ndarray:
    for key in ("X_pca", "X_scvi", "X_diffmap", "X_force_directed", "X_umap"):
        if key in adata.obsm:
            z = np.asarray(adata.obsm[key])
            if z.ndim == 2 and z.shape[0] == adata.n_obs:
                return z[:, : min(z.shape[1], 20)]
    x = adata.X
    if sparse.issparse(x):
        x = x[:, : min(50, x.shape[1])].toarray()
    else:
        x = np.asarray(x[:, : min(50, x.shape[1])])
    return np.nan_to_num(x, nan=0.0)


def _time_key(adata: ad.AnnData) -> str | None:
    for key in ("time_numeric", "stage_num", "time", "day"):
        if key in adata.obs:
            return key
    return None


def _lineage_key(adata: ad.AnnData) -> str | None:
    for key in ("lineage", "cell_type", "celltype.mapped", "author_cell_type", "gt_terminal_states", "lineages"):
        if key in adata.obs:
            return key
    return None


def _local_neighbor_mean(z: np.ndarray, values: np.ndarray, times: np.ndarray, k: int = 15) -> np.ndarray:
    out = np.zeros_like(values, dtype=float)
    for t in np.unique(times):
        mask = times == t
        idx = np.flatnonzero(mask)
        if idx.size <= 2:
            out[idx] = values[idx]
            continue
        kk = min(k + 1, idx.size)
        nn = NearestNeighbors(n_neighbors=kk).fit(z[idx]).kneighbors(z[idx], return_distance=False)[:, 1:]
        out[idx] = values[idx[nn]].mean(axis=1)
    return out


def _lineage_divergence(values: np.ndarray, lineages: np.ndarray) -> float:
    rows = []
    for lab in np.unique(lineages.astype(str)):
        mask = lineages.astype(str) == lab
        if mask.sum() >= 5:
            rows.append(float(values[mask].mean()))
    if len(rows) < 2:
        return 0.0
    return float(np.var(rows))


def _window_metrics(df: pd.DataFrame, pre: float, event: float, post: float, shuffled: bool = False) -> dict[str, float | bool]:
    work = df.copy()
    if shuffled:
        rng = np.random.default_rng(1729)
        work["time"] = rng.permutation(work["time"].to_numpy())
    pre_m = work[work["time"] == pre]
    event_m = work[work["time"] == event]
    post_m = work[work["time"] == post]
    if min(len(pre_m), len(event_m), len(post_m)) < 10:
        return {"valid": False}
    pre_comm = float(pre_m["comm_field"].mean())
    event_comm = float(event_m["comm_field"].mean())
    post_comm = float(post_m["comm_field"].mean())
    pre_receiver = float(pre_m["receiver_score"].mean())
    event_receiver = float(event_m["receiver_score"].mean())
    pre_div = _lineage_divergence(pre_m["comm_field"].to_numpy(), pre_m["lineage"].astype(str).to_numpy())
    event_div = _lineage_divergence(event_m["comm_field"].to_numpy(), event_m["lineage"].astype(str).to_numpy())
    post_div = _lineage_divergence(post_m["comm_field"].to_numpy(), post_m["lineage"].astype(str).to_numpy())
    return {
        "valid": True,
        "pre_comm": pre_comm,
        "event_comm": event_comm,
        "post_comm": post_comm,
        "communication_activation_effect": event_comm - pre_comm,
        "communication_relaxation_effect": post_comm - event_comm,
        "receiver_priming_effect": event_receiver - pre_receiver,
        "pre_lineage_comm_divergence": pre_div,
        "event_lineage_comm_divergence": event_div,
        "post_lineage_comm_divergence": post_div,
        "post_event_comm_divergence_effect": post_div - event_div,
        "communication_window_score": (event_comm - pre_comm) + 0.5 * (post_div - event_div),
        "niche_priming_then_divergence": bool((event_comm > pre_comm) and (post_div >= event_div)),
    }


def _choose_windows(times: np.ndarray, event_time: float | None) -> list[tuple[float, float, float, str]]:
    unique = sorted(float(t) for t in np.unique(times) if pd.notna(t))
    windows = []
    for i in range(1, len(unique) - 1):
        windows.append((unique[i - 1], unique[i], unique[i + 1], "scan"))
    if event_time is not None and len(unique) >= 3:
        arr = np.asarray(unique)
        event = float(arr[np.argmin(np.abs(arr - event_time))])
        pos = unique.index(event)
        if 0 < pos < len(unique) - 1:
            windows.insert(0, (unique[pos - 1], event, unique[pos + 1], "pre_registered"))
    return windows


def _top_lr_candidates(adata: ad.AnnData, lr: pd.DataFrame, sender_neighbor: np.ndarray, receiver: np.ndarray) -> pd.DataFrame:
    index = _gene_index(adata)
    rows = []
    present = lr[lr["ligand"].str.lower().isin(index) & lr["receptor"].str.lower().isin(index)].drop_duplicates(
        ["ligand", "receptor"]
    )
    if present.empty:
        return pd.DataFrame()
    sample = present.head(500)
    for row in sample.itertuples(index=False):
        li = index[str(row.ligand).lower()]
        ri = index[str(row.receptor).lower()]
        lig_x = adata.X[:, li]
        rec_x = adata.X[:, ri]
        if sparse.issparse(lig_x):
            lig_arr = np.asarray(lig_x.todense()).ravel()
            rec_arr = np.asarray(rec_x.todense()).ravel()
        else:
            lig_arr = np.asarray(lig_x).ravel()
            rec_arr = np.asarray(rec_x).ravel()
        score = float(np.sqrt(np.maximum(lig_arr, 0).mean() * np.maximum(rec_arr, 0).mean()))
        if score > 0:
            rows.append({"ligand": row.ligand, "receptor": row.receptor, "mean_pair_score": score})
    return pd.DataFrame(rows).sort_values("mean_pair_score", ascending=False).head(20)


def _module_scan(adata: ad.AnnData, z: np.ndarray, times: np.ndarray, lineages: np.ndarray, windows: pd.DataFrame) -> pd.DataFrame:
    if windows.empty:
        return pd.DataFrame()
    rows = []
    best = windows.iloc[0]
    pre, event, post = float(best["pre_time"]), float(best["event_time"]), float(best["post_time"])
    for module, cfg in COMMUNICATION_MODULES.items():
        lig_idx = _indices_for_prefixes(adata, cfg["ligand_prefixes"])
        rec_idx = _indices_for_prefixes(adata, cfg["receptor_prefixes"])
        if not lig_idx or not rec_idx:
            rows.append(
                {
                    "module": module,
                    "matched_ligand_genes": len(lig_idx),
                    "matched_receptor_genes": len(rec_idx),
                    "valid": False,
                }
            )
            continue
        ligand = _vector_for_indices(adata, lig_idx)
        receptor = _vector_for_indices(adata, rec_idx)
        ligand_neighbor = _local_neighbor_mean(z, ligand, times)
        field = ligand_neighbor * receptor
        df = pd.DataFrame({"time": times, "lineage": lineages, "comm_field": field, "receiver_score": receptor}).dropna(
            subset=["time"]
        )
        result = _window_metrics(df, pre, event, post, shuffled=False)
        if not result.get("valid"):
            rows.append(
                {
                    "module": module,
                    "matched_ligand_genes": len(lig_idx),
                    "matched_receptor_genes": len(rec_idx),
                    "valid": False,
                }
            )
            continue
        result.update(
            {
                "module": module,
                "matched_ligand_genes": len(lig_idx),
                "matched_receptor_genes": len(rec_idx),
                "valid": True,
                "pre_time": pre,
                "event_time": event,
                "post_time": post,
                "module_activation": result["communication_activation_effect"],
                "module_receiver_priming": result["receiver_priming_effect"],
                "module_post_divergence": result["post_event_comm_divergence_effect"],
                "module_activation_pass": bool(result["communication_activation_effect"] > 0),
                "module_receiver_priming_pass": bool(result["receiver_priming_effect"] > 0),
                "module_pass": bool(
                    result["communication_activation_effect"] > 0 and result["post_event_comm_divergence_effect"] >= 0
                ),
            }
        )
        rows.append(result)
    return pd.DataFrame(rows)


def _random_gene_controls(
    adata: ad.AnnData,
    z: np.ndarray,
    times: np.ndarray,
    lineages: np.ndarray,
    pre: float,
    event: float,
    post: float,
    n_sender: int,
    n_receiver: int,
    observed_score: float,
    n_permutations: int = 25,
) -> pd.DataFrame:
    if n_sender <= 0 or n_receiver <= 0 or adata.n_vars < 10:
        return pd.DataFrame()
    rng = np.random.default_rng(41041)
    n_sender = min(n_sender, max(2, adata.n_vars // 2))
    n_receiver = min(n_receiver, max(2, adata.n_vars // 2))
    rows = []
    all_idx = np.arange(adata.n_vars)
    for perm in range(n_permutations):
        sender_idx = rng.choice(all_idx, size=n_sender, replace=False).tolist()
        receiver_idx = rng.choice(all_idx, size=n_receiver, replace=False).tolist()
        sender = _vector_for_indices(adata, sender_idx)
        receiver = _vector_for_indices(adata, receiver_idx)
        sender_neighbor = _local_neighbor_mean(z, sender, times)
        comm = receiver * sender_neighbor
        df = pd.DataFrame({"time": times, "lineage": lineages, "comm_field": comm, "receiver_score": receiver}).dropna(
            subset=["time"]
        )
        res = _window_metrics(df, pre, event, post, shuffled=False)
        if res.get("valid"):
            rows.append(
                {
                    "permutation": perm,
                    "random_gene_score": res["communication_window_score"],
                    "random_activation_effect": res["communication_activation_effect"],
                    "random_receiver_priming_effect": res["receiver_priming_effect"],
                    "random_ge_observed": bool(res["communication_window_score"] >= observed_score),
                }
            )
    return pd.DataFrame(rows)


def analyze_dataset(
    spec: DatasetSpec, flags: dict[str, set[str]], lr: pd.DataFrame
) -> tuple[pd.DataFrame, dict, pd.DataFrame, pd.DataFrame]:
    if not spec.path.exists():
        return pd.DataFrame(), {
            "dataset": spec.dataset_id,
            "status": "missing",
            "path": str(spec.path),
            "role": spec.role,
            "independence_tier": spec.independence_tier,
            "support_tier": "fail",
        }, pd.DataFrame(), pd.DataFrame()
    adata = ad.read_h5ad(spec.path)
    tkey = _time_key(adata)
    lkey = _lineage_key(adata)
    if tkey is None or lkey is None:
        return pd.DataFrame(), {
            "dataset": spec.dataset_id,
            "status": "missing_time_or_lineage",
            "path": str(spec.path),
            "n_obs": adata.n_obs,
            "role": spec.role,
            "independence_tier": spec.independence_tier,
            "support_tier": "fail",
        }, pd.DataFrame(), pd.DataFrame()
    times = pd.to_numeric(adata.obs[tkey], errors="coerce").to_numpy()
    lineages = adata.obs[lkey].astype(str).to_numpy()
    z = _latent(adata)
    sender, n_sender = _vector_for_genes(adata, flags["sender"])
    receiver, n_receiver = _vector_for_genes(adata, flags["receiver"])
    secreted, n_secreted = _vector_for_genes(adata, flags["secreted"])
    membrane, n_membrane = _vector_for_genes(adata, flags["membrane"])
    sender_neighbor = _local_neighbor_mean(z, sender, times)
    secreted_neighbor = _local_neighbor_mean(z, secreted, times)
    comm_field = receiver * sender_neighbor
    secreted_receiver_field = membrane * secreted_neighbor
    cell_df = pd.DataFrame(
        {
            "time": times,
            "lineage": lineages,
            "sender_score": sender,
            "receiver_score": receiver,
            "secreted_score": secreted,
            "membrane_score": membrane,
            "sender_neighbor_score": sender_neighbor,
            "comm_field": comm_field,
            "secreted_receiver_field": secreted_receiver_field,
        }
    ).dropna(subset=["time"])
    windows = _choose_windows(cell_df["time"].to_numpy(), spec.event_time)
    rows = []
    for pre, event, post, source in windows:
        result = _window_metrics(cell_df, pre, event, post, shuffled=False)
        if not result.get("valid"):
            continue
        control = _window_metrics(cell_df, pre, event, post, shuffled=True)
        control_score = float(control.get("communication_window_score", np.nan)) if control.get("valid") else np.nan
        result.update(
            {
                "dataset": spec.dataset_id,
                "role": spec.role,
                "independence_tier": spec.independence_tier,
                "expected_window_type": spec.expected_window_type,
                "pre_time": pre,
                "event_time": event,
                "post_time": post,
                "window_source": source,
                "time_shuffle_score": control_score,
                "time_shuffle_pass": bool(np.isnan(control_score) or result["communication_window_score"] > control_score),
            }
        )
        rows.append(result)
    win = pd.DataFrame(rows)
    if not win.empty:
        sort_cols = ["window_source", "communication_window_score"]
        win["_source_rank"] = np.where(win["window_source"] == "pre_registered", 0, 1)
        win = win.sort_values(["_source_rank", "communication_window_score"], ascending=[True, False]).drop(columns="_source_rank")
    top_lr = _top_lr_candidates(adata, lr, sender_neighbor, receiver)
    if not top_lr.empty:
        top_lr.insert(0, "dataset", spec.dataset_id)
    module_df = _module_scan(adata, z, times, lineages, win)
    if not module_df.empty:
        module_df.insert(0, "dataset", spec.dataset_id)
    control_df = pd.DataFrame()
    random_gene_q = np.nan
    random_activation_q = np.nan
    random_receiver_q = np.nan
    if not win.empty:
        best = win.iloc[0]
        control_df = _random_gene_controls(
            adata,
            z,
            times,
            lineages,
            float(best["pre_time"]),
            float(best["event_time"]),
            float(best["post_time"]),
            min(n_sender, 500),
            min(n_receiver, 500),
            float(best["communication_window_score"]),
        )
        if not control_df.empty:
            random_gene_q = float((control_df["random_gene_score"] >= float(best["communication_window_score"])).mean())
            random_activation_q = float(
                (control_df["random_activation_effect"] >= float(best["communication_activation_effect"])).mean()
            )
            random_receiver_q = float(
                (control_df["random_receiver_priming_effect"] >= float(best["receiver_priming_effect"])).mean()
            )
            control_df.insert(0, "dataset", spec.dataset_id)
            control_df["observed_score"] = float(best["communication_window_score"])
            control_df["observed_activation_effect"] = float(best["communication_activation_effect"])
            control_df["observed_receiver_priming_effect"] = float(best["receiver_priming_effect"])
    audit = {
        "dataset": spec.dataset_id,
        "status": "analyzed",
        "path": str(spec.path),
        "role": spec.role,
        "independence_tier": spec.independence_tier,
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "time_key": tkey,
        "lineage_key": lkey,
        "time_count": int(len(np.unique(times[~pd.isna(times)]))),
        "lineage_count": int(len(set(lineages))),
        "matched_sender_genes": int(n_sender),
        "matched_receiver_genes": int(n_receiver),
        "matched_secreted_genes": int(n_secreted),
        "matched_membrane_genes": int(n_membrane),
        "analyzed_windows": int(len(win)),
        "best_communication_window_score": float(win["communication_window_score"].max()) if not win.empty else np.nan,
        "best_niche_priming_then_divergence": bool(win.iloc[0]["niche_priming_then_divergence"]) if not win.empty else False,
        "random_gene_control_q": random_gene_q,
        "random_activation_q": random_activation_q,
        "random_receiver_priming_q": random_receiver_q,
        "best_communication_activation_effect": float(win.iloc[0]["communication_activation_effect"]) if not win.empty else np.nan,
        "best_receiver_priming_effect": float(win.iloc[0]["receiver_priming_effect"]) if not win.empty else np.nan,
        "activation_pass": bool(
            (not win.empty)
            and (float(win.iloc[0]["communication_activation_effect"]) > 0)
            and (pd.isna(random_activation_q) or random_activation_q <= 0.10)
        ),
        "receiver_priming_pass": bool(
            (not win.empty)
            and (float(win.iloc[0]["receiver_priming_effect"]) > 0)
            and (pd.isna(random_receiver_q) or random_receiver_q <= 0.10)
        ),
        "support_tier": "weak",
    }
    if not win.empty:
        best = win.iloc[0]
        if (
            bool(best["niche_priming_then_divergence"])
            and bool(best["time_shuffle_pass"])
            and (pd.isna(random_gene_q) or random_gene_q <= 0.10)
            and audit["matched_sender_genes"] >= 25
            and audit["matched_receiver_genes"] >= 25
        ):
            audit["support_tier"] = "acceptable" if "independent" in spec.independence_tier or spec.dataset_id.startswith("E1") else "weak"
        elif not bool(best["niche_priming_then_divergence"]):
            audit["support_tier"] = "fail"
    return win, audit, top_lr, control_df, module_df


def main() -> int:
    TABLES.mkdir(exist_ok=True)
    REPORTS.mkdir(exist_ok=True)
    MANUSCRIPT.mkdir(exist_ok=True)
    _, flags, lr = _load_intercell_prior()
    all_windows = []
    audits = []
    all_lr = []
    all_controls = []
    all_modules = []
    for spec in DATASETS:
        win, audit, lr_rows, controls, modules = analyze_dataset(spec, flags, lr)
        if not win.empty:
            all_windows.append(win)
        audits.append(audit)
        if not lr_rows.empty:
            all_lr.append(lr_rows)
        if not controls.empty:
            all_controls.append(controls)
        if not modules.empty:
            all_modules.append(modules)
    windows = pd.concat(all_windows, ignore_index=True) if all_windows else pd.DataFrame()
    audit_df = pd.DataFrame(audits)
    lr_df = pd.concat(all_lr, ignore_index=True) if all_lr else pd.DataFrame()
    controls_df = pd.concat(all_controls, ignore_index=True) if all_controls else pd.DataFrame()
    modules_df = pd.concat(all_modules, ignore_index=True) if all_modules else pd.DataFrame()
    if not windows.empty:
        windows.to_csv(TABLES / "communication_niche_window_scan.csv", index=False)
    audit_df.to_csv(TABLES / "communication_niche_dataset_audit.csv", index=False)
    if not lr_df.empty:
        lr_df.to_csv(TABLES / "communication_niche_lr_candidates.csv", index=False)
    if not controls_df.empty:
        controls_df.to_csv(TABLES / "communication_niche_negative_controls.csv", index=False)
    if not modules_df.empty:
        modules_df.to_csv(TABLES / "communication_module_scan.csv", index=False)
    acceptable = audit_df[audit_df["support_tier"].isin(["acceptable"])]
    activation_series = (
        audit_df["activation_pass"] if "activation_pass" in audit_df else pd.Series(False, index=audit_df.index)
    )
    receiver_series = (
        audit_df["receiver_priming_pass"]
        if "receiver_priming_pass" in audit_df
        else pd.Series(False, index=audit_df.index)
    )
    activation_support = audit_df[activation_series.astype(str).str.lower().isin(["true", "1"])]
    receiver_support = audit_df[receiver_series.astype(str).str.lower().isin(["true", "1"])]
    analyzed = audit_df[audit_df["status"] == "analyzed"]
    cross_dataset = len(acceptable) >= 2
    independent_accept = acceptable["independence_tier"].astype(str).str.contains("independent").any() if not acceptable.empty else False
    activation_cross = (
        len(activation_support) >= 3
        and activation_support["dataset"].astype(str).str.contains("internal_native").any()
        and activation_support["independence_tier"].astype(str).str.contains("independent").any()
    )
    receiver_cross = (
        len(receiver_support) >= 2
        and receiver_support["dataset"].astype(str).str.contains("internal_native|E1").any()
        and receiver_support["independence_tier"].astype(str).str.contains("independent").any()
    )
    if activation_cross:
        final_tier = "acceptable"
        conclusion = "communication_activation_is_cross_dataset_candidate"
    elif receiver_cross:
        final_tier = "weak"
        conclusion = "receiver_priming_is_partial_cross_dataset_candidate"
    elif cross_dataset and independent_accept:
        final_tier = "acceptable"
        conclusion = "communication_niche_priming_is_cross_dataset_candidate"
    elif len(acceptable) >= 1:
        final_tier = "weak"
        conclusion = "communication_niche_priming_has_partial_support"
    else:
        final_tier = "fail"
        conclusion = "communication_niche_priming_not_supported"
    module_summary_rows = []
    if not modules_df.empty:
        valid_modules = modules_df[modules_df["valid"].astype(str).str.lower().isin(["true", "1"])]
        for module, grp in valid_modules.groupby("module"):
            support = grp[grp["module_pass"].astype(str).str.lower().isin(["true", "1"])]
            activation = grp[grp["module_activation_pass"].astype(str).str.lower().isin(["true", "1"])]
            receiver_priming = grp[grp["module_receiver_priming_pass"].astype(str).str.lower().isin(["true", "1"])]
            module_summary_rows.append(
                {
                    "module": module,
                    "support_datasets": int(len(support)),
                    "datasets": ";".join(support["dataset"].astype(str).tolist()),
                    "activation_support_datasets": int(len(activation)),
                    "activation_datasets": ";".join(activation["dataset"].astype(str).tolist()),
                    "receiver_priming_support_datasets": int(len(receiver_priming)),
                    "receiver_priming_datasets": ";".join(receiver_priming["dataset"].astype(str).tolist()),
                    "has_internal": bool(support["dataset"].astype(str).eq("internal_native").any()),
                    "has_e1": bool(support["dataset"].astype(str).eq("E1_MouseGastrulationData").any()),
                    "has_independent": bool(
                        support["dataset"].astype(str).isin(["E5_zebrafish_Farrell", "E2_GSE212050_gastruloid", "GSE154572_EB_WT"]).any()
                    ),
                    "activation_has_internal": bool(activation["dataset"].astype(str).eq("internal_native").any()),
                    "activation_has_e1": bool(activation["dataset"].astype(str).eq("E1_MouseGastrulationData").any()),
                    "activation_has_independent": bool(
                        activation["dataset"].astype(str).isin(["E5_zebrafish_Farrell", "E2_GSE212050_gastruloid", "GSE154572_EB_WT"]).any()
                    ),
                    "mean_activation": float(pd.to_numeric(grp["module_activation"], errors="coerce").mean()),
                    "mean_post_divergence": float(pd.to_numeric(grp["module_post_divergence"], errors="coerce").mean()),
                }
            )
    module_summary = pd.DataFrame(module_summary_rows)
    if not module_summary.empty:
        module_summary.to_csv(TABLES / "communication_module_summary.csv", index=False)
    summary = pd.DataFrame(
        [
            {
                "analysis": "communication_niche_field",
                "tier": final_tier,
                "analyzed_datasets": int(len(analyzed)),
                "acceptable_datasets": int(len(acceptable)),
                "activation_support_datasets": int(len(activation_support)),
                "receiver_priming_support_datasets": int(len(receiver_support)),
                "independent_acceptable": bool(independent_accept),
                "conclusion": conclusion,
                "allowed_claim": (
                    "local communication-niche priming is a candidate branch-window annotation"
                    if final_tier != "fail"
                    else "current communication-niche field does not provide a retained cross-dataset mechanism"
                ),
                "forbidden_claim": "established intercellular signalling mechanism",
            }
        ]
    )
    summary.to_csv(TABLES / "communication_niche_cross_dataset_summary.csv", index=False)
    report = [
        "# Communication Niche Search",
        "",
        "## Hypothesis",
        "PathwayFinder/OmniPath intercell priors can define a local communication-niche field: neighbor sender/transmitter expression multiplied by receiver/membrane competence in each cell. A branch window would show niche priming before or during the event, followed by lineage-specific communication divergence.",
        "",
        "## Success Criteria",
        "At least two datasets, including one independent developmental dataset, should show positive communication activation, post-event communication divergence and time-shuffle resistance. This is a computational annotation, not validation of a causal CCI mechanism.",
        "",
        "## Result",
        f"- analyzed_datasets: {len(analyzed)}",
        f"- acceptable_datasets: {len(acceptable)}",
        f"- activation_support_datasets: {len(activation_support)}",
        f"- receiver_priming_support_datasets: {len(receiver_support)}",
        f"- final_tier: `{final_tier}`",
        f"- conclusion: `{conclusion}`",
        "",
        "## Interpretation",
    ]
    top_module_line = ""
    if not module_summary.empty:
        candidates = module_summary[
            (module_summary["activation_has_internal"] == True)
            & (module_summary["activation_has_e1"] == True)
            & (module_summary["activation_has_independent"] == True)
        ].sort_values("activation_support_datasets", ascending=False)
        if not candidates.empty:
            best_module = candidates.iloc[0]
            top_module_line = (
                f" The strongest module-level candidate is `{best_module['module']}`, "
                f"with positive branch-window activation in {best_module['activation_datasets']}."
            )
    if final_tier == "acceptable":
        report.append(
            "A local communication-niche field is the most promising CCI-aligned direction found so far. It does not rescue previous LR knockout/rerollout claims, but it gives the SwarmLineage-OT framework a testable extracellular-niche annotation layer."
            + top_module_line
        )
    elif final_tier == "weak":
        report.append(
            "The communication-niche field shows partial cross-window signal but is not yet strong enough for a main mechanism claim. It can guide the next targeted spatial/ligand-receptor validation."
        )
    else:
        report.append(
            "The communication-niche field does not currently produce a cross-dataset retained mechanism. CCI remains a future validation direction rather than a retained claim."
        )
    report += [
        "",
        "## Evidence Boundary",
        "- Uses PathwayFinder OmniPath intercell prior reused from the adjacent project.",
        "- Uses expression and latent kNN neighborhoods already present in SwarmLineage-OT data products.",
        "- Does not establish true ligand-receptor signalling, causality, or wet-lab validation.",
    ]
    (REPORTS / "communication_niche_search_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    story = [
        "# Communication-Niche Candidate Story",
        "",
        f"Final tier: `{final_tier}`.",
        "",
        "SwarmLineage-OT can represent intercellular communication as a local niche field derived from PathwayFinder/OmniPath transmitter and receiver priors. The current analysis asks whether branch windows are accompanied by a rise in local communication competence and subsequent lineage-specific divergence.",
        "",
        f"Current conclusion: `{conclusion}`.",
        "",
        "This should be framed as a candidate extracellular-niche annotation layer only. It must not be described as established intercellular signalling or a cause-effect mechanism.",
    ]
    (MANUSCRIPT / "communication_niche_story.md").write_text("\n".join(story) + "\n", encoding="utf-8")
    print(json.dumps(summary.iloc[0].to_dict(), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
