# Leakage Audit

- split_mode: `strict_time_holdout`
- holdout_time: `15.0`
- preprocessing HVG/SVD fit cells: 7000 train cells only
- evaluation-only held-out cells present in output for scoring: 1000
- held-out cells are transformed by the train-fitted SVD but do not contribute to HVG selection or SVD fitting.
- teacher construction must additionally exclude `split_role == eval_holdout` cells.
