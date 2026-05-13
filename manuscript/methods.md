# Methods

See `manuscript/theory.md` and `manuscript/methods_theory.tex` for the mathematical formulation. The implemented pipeline uses AnnData preprocessing, PCA latent states, adjacent-stage entropic OT, PyTorch velocity fitting, stochastic finite-agent rollout, and agent-level ablations.

Evaluation is organized around teacher fidelity, emergent-law robustness and mechanistic usefulness. `M0b_ot_interpolation` is the OT teacher/reference interpolation. It is not treated as a baseline that SwarmLineage-OT must beat. Each gate is reported as fail, weak, acceptable or strong, with negative controls, seed stability and rollout support recorded for every discovery law.

Statistical summaries use paired seed-level comparisons, bootstrap confidence intervals, permutation controls and Benjamini-Hochberg correction where applicable.
