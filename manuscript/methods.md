# Methods

For each developmental atlas dataset, cells were standardized to `time_numeric`, `time_point`, `lineage` and `cell_type`, then stratified by time and lineage. PCA was fit within each external dataset only. Native moscot TemporalProblem extraction was attempted through the clean native environment; fallback status is recorded per dataset if native extraction fails.

The branch-window detector scans adjacent time triples and selects the window maximizing a pre-registered score combining transient lineage-separation contraction, post-event divergence, local velocity alignment and transition entropy. Negative controls shuffle time labels, velocity, lineage labels, teacher velocity or neighbor graph. Baselines include fate entropy alone, lineage separation alone, cell-type composition change and alignment alone.
