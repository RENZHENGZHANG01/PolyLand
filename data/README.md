# Data provenance and schema

All gas-permeability values are in **Barrer**. Processed modeling targets also
include `*_log10`, defined as `log10(permeability in Barrer)`.

## Raw inputs

| File | Rows | Description |
|---|---:|---|
| `raw/linear_permeability.csv` | 838 | Linear-polymer permeability measurements |
| `raw/ladder_permeability.csv` | 143 | Literature-curated ladder and semi-ladder polymer measurements |
| `raw/md_ffv_index.csv` | 563 | Polymer identifiers and molecular-dynamics metadata used by the analysis |

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

## Sources and reuse

The ladder-polymer literature references are listed in the PolyLand
manuscript. The linear data originate from the gas-permeability data used in
the cited POINT2/MSA workflow. The repository's MIT License covers the code;
third-party data and publications remain subject to their original terms.
