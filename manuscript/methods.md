# Methods

See `manuscript/theory.md` and `manuscript/methods_theory.tex` for the mathematical formulation. The implemented pipeline uses AnnData preprocessing, PCA latent states, adjacent-stage entropic OT, PyTorch velocity fitting, stochastic finite-agent rollout, and agent-level ablations.

Evaluation is organized around teacher fidelity, emergent-law robustness and mechanistic usefulness. `M0b_ot_interpolation` is the OT teacher/reference interpolation. It is not treated as a baseline that SwarmLineage-OT must beat.

Statistical summaries use paired seed-level comparisons, bootstrap confidence intervals and Benjamini-Hochberg correction where applicable.
