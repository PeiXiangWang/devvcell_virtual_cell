#!/usr/bin/env bash
set -euo pipefail

QUICK_FLAG=""
if [[ "${1:-}" == "--quick-fixture" ]]; then
  QUICK_FLAG="--quick-fixture"
fi

python -m src.data.preprocess --config configs/data.yaml ${QUICK_FLAG}
python -m src.ot_teacher.run_moscot --config configs/ot_teacher.yaml ${QUICK_FLAG}
python -m src.ot_teacher.run_wot --config configs/ot_teacher.yaml ${QUICK_FLAG}
python -m src.ot_teacher.build_teacher --config configs/ot_teacher.yaml ${QUICK_FLAG}
python -m src.train.train_model --config configs/model.yaml ${QUICK_FLAG}
python -m src.train.ablations --config configs/train.yaml ${QUICK_FLAG}
python -m src.train.evaluate --config configs/train.yaml ${QUICK_FLAG}
python -m src.perturbation.simulate_lr_knockout --config configs/model.yaml ${QUICK_FLAG}
python -m src.perturbation.simulate_gene_perturb --config configs/model.yaml ${QUICK_FLAG}
python -m src.perturbation.evaluate_perturbation --config configs/model.yaml ${QUICK_FLAG}
python -m src.plotting.figures_main
python -m src.plotting.figures_extended
python -m src.utils.reproducibility
