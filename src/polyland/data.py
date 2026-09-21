"""Build the gas-permeability tables used by the PolyLand analysis.

The MD/FFV index defines the subset used in Table 1 and in the predictive
experiments.  Ladder records are matched to the literature source with
``PID + SMILES + N2``.  The N2 value is part of the key because a small
number of polymer identifiers have multiple reported measurements.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

GASES = ("O2", "N2", "H2", "CH4", "CO2")
INDEX_COLUMNS = (
    "PID",
    "Type",
    "SMILES",
    "density_MD",
    "density_MD_std",
    "n_repeat_vdw",
    "vdw",
    "FFV",
)


def _read_sources(raw_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    linear = pd.read_csv(raw_dir / "linear_permeability.csv")
    ladder = pd.read_csv(raw_dir / "ladder_permeability.csv")
    index = pd.read_csv(raw_dir / "md_ffv_index.csv")

    required_linear = {"PID", "SMILES", *GASES}
    required_ladder = {"PID", "SMILES", "N2", *GASES}
    required_index = {*INDEX_COLUMNS, "match_N2_Barrer"}
    for name, frame, required in (
        ("linear_permeability.csv", linear, required_linear),
        ("ladder_permeability.csv", ladder, required_ladder),
        ("md_ffv_index.csv", index, required_index),
    ):
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(f"{name} is missing columns: {sorted(missing)}")
    return linear, ladder, index


def _write_gas_tables(frame: pd.DataFrame, kind: str, output_dir: Path) -> dict[str, Path]:
    written: dict[str, Path] = {}
    for gas in GASES:
        columns = [*INDEX_COLUMNS, gas]
        gas_frame = frame.loc[frame[gas].notna(), columns].copy()
        if (gas_frame[gas] <= 0).any():
            raise ValueError(f"{kind} {gas} contains non-positive permeability")
        gas_frame[f"{gas}_log10"] = np.log10(gas_frame[gas].astype(float))
        path = output_dir / f"final_{kind}_data_{gas}.csv"
        gas_frame.to_csv(path, index=False)
        written[f"{kind}_{gas}"] = path
    return written


def build_processed_tables(raw_dir: str | Path, output_dir: str | Path) -> dict[str, Path]:
    """Rebuild the ten per-gas linear/ladder tables.

    The linear records are joined by their unique ``PID``.  Ladder records
    are joined by ``PID``, ``SMILES`` and the reported N2 permeability to
    preserve distinct measurements for repeated polymer identifiers.
    """

    raw_dir = Path(raw_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    linear, ladder, index = _read_sources(raw_dir)

    index_linear = index.loc[~index["PID"].astype(str).str.startswith("LP")].copy()
    index_ladder = index.loc[index["PID"].astype(str).str.startswith("LP")].copy()

    linear_values = linear[["PID", *GASES]].copy()
    if linear_values["PID"].duplicated().any():
        raise ValueError("Linear source contains duplicate PID values")
    linear_merged = index_linear[list(INDEX_COLUMNS)].merge(
        linear_values,
        on="PID",
        how="left",
        validate="one_to_one",
    )

    ladder_values = ladder[["PID", "SMILES", "N2", "O2", "H2", "CH4", "CO2"]].copy()
    if ladder_values.duplicated(["PID", "SMILES", "N2"]).any():
        raise ValueError("Ladder source has an ambiguous PID/SMILES/N2 key")
    ladder_merged = index_ladder.rename(columns={"match_N2_Barrer": "N2"}).merge(
        ladder_values,
        on=["PID", "SMILES", "N2"],
        how="left",
        validate="many_to_one",
    )
    ladder_merged = ladder_merged[[*INDEX_COLUMNS, *GASES]]

    paired = ladder_merged.dropna(subset=["N2", "CH4"])
    if len(paired) and np.allclose(paired["N2"], paired["CH4"]):
        raise ValueError(
            "Ladder N2 and CH4 are unexpectedly identical; verify the input tables."
        )

    written = _write_gas_tables(linear_merged, "linear", output_dir)
    written.update(_write_gas_tables(ladder_merged, "ladder", output_dir))
    return written
