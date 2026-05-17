# Reviewer Attack Matrix

| attack | current answer | evidence | remaining gap | allowed claim |
|---|---|---|---|---|
| Does the signal generalize beyond E1? | Not yet at acceptable tier. Final EB/spatial sprint remains weak or blocked. | final external atlas table | independent annotated developmental dataset needed | internal/E1-supported hypothesis |
| Is GSE154572 a validation dataset? | No. It is an independent EB time-series feasibility row with cluster-proxy labels. | GSE154572 metadata audit | curated lineage/cell-type labels absent | weak stress test only |
| Is spatial condensation validated? | No. STDS0000074 was verified and one h5ad inspected, but cell-level multi-stage annotations were unavailable in the inspected object. | external data integrity audit | spatial cell-state matrix with stage and cell type needed | spatial validation remains future work |
| Is this clone fate prediction? | No. Clone-aware tests remain stress tests and do not establish fate-diversification prediction. | clone audit tables | richer clone/time data required | clone line is future work |
| Is the model causal or superior to OT? | No. The framework realizes an OT pseudo-lineage and audits order parameters; it is not a causal or superiority claim. | claim audit | experimental perturbation absent | computational order-parameter framework |
