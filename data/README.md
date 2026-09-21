# Data provenance and schema

All gas permeability values are in **Barrer**. The processed modeling targets
also include `*_log10`, defined as `log10(permeability in Barrer)`.

## Raw inputs

| File | Rows | Purpose |
|---|---:|---|
| `raw/linear_permeability.csv` | 838 | Linear-polymer permeability source with stable `PID` values |
| `raw/ladder_permeability.csv` | 143 | Literature-curated ladder/semi-ladder source, including source links |
| `raw/md_ffv_index.csv` | 563 | Records with MD density/FFV metadata used to define the common analysis subset |

`md_ffv_index.csv` deliberately contains `match_N2_Barrer` but no CH4 column.
N2 is retained only as a measurement-level join key for duplicate ladder PIDs.
This prevents the historically contaminated CH4 values from entering the
repository while preserving a traceable match to the curated source.

## Processed tables

`scripts/prepare_data.py` writes `final_{linear|ladder}_data_{gas}.csv` for O2,
N2, H2, CH4, and CO2. Each table contains:

- `PID`, `Type`, and polymer repeat-unit `SMILES`;
- MD density, van der Waals, and FFV metadata;
- permeability in Barrer and its base-10 logarithm.

The processed rows are the intersection of the relevant permeability label
and the MD/FFV eligibility index. That selection is why Table 1 counts are
smaller than the number of non-null measurements in the full 143-row ladder
source.

The literature references for the ladder collection are listed in the
PolyLand manuscript. The linear data originate from the gas-permeability data
used in the cited POINT2/MSA workflow. The repository's MIT License covers the
code, not any additional rights in third-party publications.

