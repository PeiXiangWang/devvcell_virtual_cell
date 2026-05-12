# Literature Positioning

This scan was last updated on 2026-05-11. It is used to prevent novelty overclaiming.

## Anchors

- **moscot**: multi-omics single-cell optimal transport for scalable temporal, spatial and multimodal mapping. Primary reference: [Mapping cells through time and space with moscot](https://www.nature.com/articles/s41586-024-08453-2).
- **Waddington-OT**: infers temporal couplings and ancestor-descendant/fate flows from destructive time-course scRNA-seq snapshots. Primary reference: [Schiebinger et al., Cell 2019](https://doi.org/10.1016/j.cell.2019.01.006).
- **TIGON**: dynamic unbalanced OT with growth and GRN reconstruction, so SwarmLineage-OT cannot claim to be the first growth-aware OT method. Primary reference: [Nature Machine Intelligence](https://www.nature.com/articles/s42256-023-00763-w).
- **COMMOT**: collective OT for spatial cell-cell communication, so this project cannot claim first OT-based CCI. Primary reference: [Nature Methods](https://www.nature.com/articles/s41592-022-01728-4).
- **CellRank 2**: unified multiview fate mapping and fate probabilities. Primary reference: [Nature Methods](https://www.nature.com/articles/s41592-024-02303-9).
- **scIMF and IADOT**: interaction-aware or mean-field dynamics are active prior art. scIMF models interacting multicellular dynamics with McKean-Vlasov SDEs ([PLOS Computational Biology, 2026](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1013916)); IADOT adds cell-cell interaction networks to dynamic OT and flow matching ([OpenReview, 2025/2026](https://openreview.net/forum?id=8H1L06TGNS)).
- **Virtual cell / perturbation prediction**: GEARS predicts transcriptional perturbation responses ([Nature Biotechnology](https://www.nature.com/articles/s41587-023-01905-6)); scGPT and scFoundation are single-cell foundation models ([scGPT](https://www.nature.com/articles/s41592-024-02201-0), [scFoundation](https://www.nature.com/articles/s41592-024-02305-7)).

## Novelty Matrix

| Capability | WOT | moscot | TIGON | COMMOT | CellRank2 | TrajectoryNet/MIOFlow | scIMF/IADOT | GEARS/scGPT/scFoundation | SwarmLineage-OT |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| OT lineage teacher | yes | yes | partial | no | partial | partial | partial | no | yes |
| executable agent-based virtual cell | no | no | no | no | no | no | partial mean-field | no | yes |
| local swarm rules | no | no | no | no | no | no | partial interactions | no | yes |
| birth/death hazard | growth prior | unbalanced OT | yes | no | no | no | partial | no | yes |
| adaptive diffusion | no | solver diffusion | dynamic model | no | no | flow/SDE variants | yes | no | yes |
| CCI-aware interactions | no | possible via costs | reported application | yes | no | no | yes | no | yes |
| finite-population stochastic simulation | no | no | no | no | Markov chain only | continuous flow | mostly mean-field | no | yes |
| perturbation counterfactuals | limited | possible mapping | limited | LR perturbation | no | limited | yes | yes expression | yes control-layer |
| interpretable lineage-constrained energy | no | no | no | no | no | limited | limited | no | yes |
| held-out temporal reconstruction | yes | yes | yes | no | yes | yes | yes | no | yes |
| external/lineage/perturbation validation | yes in paper | yes in paper | yes in paper | yes in paper | yes in paper | varies | emerging | yes perturbation | required, not yet complete |

## Non-Substitutable Claim

The defensible novelty target is not “OT with growth” or “OT with CCI.” The target is: converting OT-inferred developmental couplings into supervision for a finite-agent swarm virtual-cell simulator with explicit lineage-constrained birth, diffusion, local interaction rules, CCI modulation and perturbation control. This remains a hypothesis until native baselines, external data and validation pass the project gates.

