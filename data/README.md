# Data provenance and schema

All gas-permeability values are in **Barrer**. Processed modeling targets also
include `*_log10`, defined as `log10(permeability in Barrer)`.

## Raw inputs

| File | Rows | Description |
|---|---:|---|
| `raw/linear_permeability.csv` | 838 | Linear-polymer permeability measurements |
| `raw/ladder_permeability.csv` | 143 | Literature-curated ladder and semi-ladder polymer measurements |
| `raw/md_ffv_index.csv` | 563 | Polymer identifiers and molecular-dynamics metadata used by the analysis |
| `screening/polyland_screening_candidates.csv` | 77,246 | All post-novelty-filter candidates evaluated by the five-gas ensemble screen |
| `screening/selected_candidates.csv` | 22 | SLP identifiers and SMILES for candidates highlighted in the manuscript |
| `processed/ad_ladder_training_structures.csv` | 178 | Unique gas-specific ladder structures in the original model-training partitions, used to reproduce the applicability-domain analysis |

The `match_N2_Barrer` field in `md_ffv_index.csv` is a measurement-level key
used to distinguish records that share a polymer identifier.

## Processed tables

Running `python scripts/prepare_data.py` writes one linear and one ladder table
for each of O2, N2, H2, CH4, and CO2. Each table contains:

- `PID`: polymer identifier;
- `Type`: polymer class;
- `SMILES`: repeat-unit representation;
- `density_MD` and `density_MD_std`: molecular-dynamics density descriptors;
- `n_repeat_vdw`, `vdw`, and `FFV`: structural and free-volume descriptors;
- the gas permeability in Barrer;
- the corresponding `*_log10` modeling target.

Rows without a measured value for the selected gas are omitted from that
gas-specific table.

`processed/ad_ladder_training_structures.csv` records the unique ladder and
semi-ladder structures present in the original 80% training partition for
each gas. The applicability-domain workflow combines these entries with the
corresponding complete linear table to construct the hybrid domain.

## Candidate screening table

`screening/polyland_screening_candidates.csv` is a flat, machine-readable
table with one row per candidate entering the final screen. The primary fields
are:

- `candidate_id`: stable row identifier assigned by the table builder;
- `generation_method`: `DiT`, `LLM`, or `rule-based`;
- `SMILES`: polymer repeat-unit representation;
- `selected`: `True` for candidates highlighted in the manuscript;
- `selection_id`: SLP figure identifier for selected rows;
- `<gas>_log10_Barrer_ensemble_mean`: mean of the top-five model predictions;
- `<gas>_log10_Barrer_ensemble_sd`: standard deviation across those five point
  predictions;
- `<gas1>_<gas2>_selectivity`: predicted ideal selectivity computed as the
  permeability ratio.

Permeability in Barrer can be recovered as
`10 ** <gas>_log10_Barrer_ensemble_mean`.

The five separation tasks are O2/N2, H2/CH4, H2/N2, CO2/CH4, and CO2/N2.
The ensemble standard-deviation fields quantify between-model spread; they are
not quantile-based prediction intervals.

## Sources and reuse

The ladder-polymer literature references are listed in the PolyLand
manuscript. The linear data originate from the gas-permeability data used in
the cited POINT2/MSA workflow. The repository's MIT License covers the code;
third-party data and publications remain subject to their original terms.
