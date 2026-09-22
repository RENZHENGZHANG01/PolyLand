# Applicability-domain outputs

This directory contains a structure-based applicability-domain analysis of
the complete post-novelty-filter screening pool.

The analysis uses 2048-bit Morgan fingerprints with radius 2 and Tanimoto
similarity. For each gas, two reference domains are evaluated:

- `ladder`: unique structures in the original ladder-polymer training split;
- `hybrid`: the ladder training structures plus the complete linear-polymer
  training data.

The similarity threshold for each domain is the fifth percentile of the
leave-one-out nearest-neighbor similarity distribution among its training
structures. A candidate is inside a domain when its nearest-training-neighbor
similarity is greater than or equal to this threshold.

## Files

- `candidate_applicability_domain.csv.gz`: candidate-level similarities,
  nearest training PIDs, polymer types, and domain flags for all five gases;
- `selected_candidate_applicability_domain.csv`: the same fields for the 22
  screening-selected candidates;
- `applicability_domain_thresholds.csv`: gas- and domain-specific cutoffs and
  training-set diagnostics;
- `applicability_domain_summary.csv`: domain coverage by gas, generation
  method, and selection status;
- `applicability_domain_summary.pdf`: vector summary figure;
- `applicability_domain_summary.png`: raster preview of the figure.

Run the workflow from the repository root after installing the `ml` optional
dependencies:

```bash
python -m pip install -e '.[ml]'
python scripts/applicability_domain.py
```

The analysis is a structural-domain diagnostic. It does not replace
model-specific uncertainty quantification or experimental validation.
