#!/usr/bin/env bash
set -euo pipefail

python -m src.data.preprocess --config configs/data.yaml
python -m src.ot_teacher.run_moscot --config configs/ot_teacher.yaml
python -m src.ot_teacher.run_wot --config configs/ot_teacher.yaml
python -m src.ot_teacher.build_teacher --config configs/ot_teacher.yaml
python -m src.train.train_model --config configs/model.yaml
python -m src.train.ablations --config configs/train.yaml
python -m src.train.evaluate --config configs/train.yaml
python -m src.perturbation.simulate_lr_knockout --config configs/model.yaml
python -m src.perturbation.simulate_gene_perturb --config configs/model.yaml
python -m src.perturbation.evaluate_perturbation --config configs/model.yaml
python -m src.plotting.figures_main
python -m src.plotting.figures_extended
python -m src.utils.reproducibility

