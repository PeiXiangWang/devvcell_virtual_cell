# SwarmLineage-OT

SwarmLineage-OT is a research prototype for converting single-cell optimal-transport pseudo-lineages into an executable finite-agent virtual-cell population model.

The current branch is intentionally conservative: if the full swarm model does not beat the strongest baseline, the reports say so. Fallback teachers are labelled `toy_sinkhorn_fallback` and must not be described as native moscot or WOT results.

## What Is Implemented

- strict held-out-time preprocessing with leakage audit;
- native-backend status tracking for moscot/WOT and explicit toy Sinkhorn fallback;
- separate ablation definitions for linear interpolation, OT interpolation, intrinsic neural dynamics, teacher-velocity dynamics and swarm modules;
- `SwarmLineageDynamics`, a trainable PyTorch module with intrinsic velocity, teacher velocity, swarm terms, stochastic birth/death hazards, adaptive diffusion, CCI modulation and memory-field support;
- multi-step finite-agent rollout with birth/death event logs;
- baseline execution matrix and scientific gap audit;
- synthetic quick fixture for CI and smoke testing.

## Quick Reproduction

Unix/WSL:

```bash
bash reproducibility/run_all.sh --quick-fixture
```

Windows PowerShell:

```powershell
powershell -ExecutionPolicy Bypass -File .\reproducibility\run_all.ps1 -QuickFixture
```

Main real-data path:

```bash
python -m src.data.preprocess --config configs/data.yaml
python -m src.ot_teacher.run_moscot --config configs/ot_teacher.yaml
python -m src.ot_teacher.run_wot --config configs/ot_teacher.yaml
python -m src.ot_teacher.build_teacher --config configs/ot_teacher.yaml
python -m src.train.train_model --config configs/model.yaml
python -m src.train.ablations --config configs/train.yaml
python -m src.train.evaluate --config configs/train.yaml
```

## Key Reports

- `reports/leakage_audit.md`
- `reports/scientific_gap_audit.md`
- `reports/baseline_execution_matrix.csv`
- `reports/module_contribution_audit.md`
- `manuscript/final_retained_results_and_methods.md`

## Archived Prior Work

Earlier DevGuard and DevSpectrum materials remain in `docs/`, `scripts/devguard`, `scripts/devspectrum`, and `src/devguard`/`src/devspectrum` for provenance. They are no longer the root README focus of this branch.
