# Leakage Audit

- Internal strict holdout remains governed by configs/train.yaml and configs/data.yaml.
- E1 preprocessing was fit only on the external MouseGastrulationData component.
- E1 does not use internal teacher information.
- L1 and E2 analyses use separate external files and do not import internal labels or teacher couplings.
- Holdout bridge edges, where present in internal runs, must be interpreted as bridge edges rather than ordinary adjacent observed edges.
