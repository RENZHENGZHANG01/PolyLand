"""Compute and format the descriptive statistics reported in Table 1."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .data import GASES

ROWS = ("Count", "Mean", "Median", "Min", "Max")


def table1_statistics(processed_dir: str | Path) -> pd.DataFrame:
    """Return Table 1 as a DataFrame with a two-level gas/type column index."""

    processed_dir = Path(processed_dir)
    columns: dict[tuple[str, str], list[float]] = {}
    for gas in GASES:
        for kind in ("Linear", "Ladder"):
            path = processed_dir / f"final_{kind.lower()}_data_{gas}.csv"
            frame = pd.read_csv(path)
            values = pd.to_numeric(frame[gas], errors="raise")
            columns[(gas, kind)] = [
                int(values.count()),
                float(values.mean()),
                float(values.median()),
                float(values.min()),
                float(values.max()),
            ]
    result = pd.DataFrame(columns, index=ROWS)
    result.columns = pd.MultiIndex.from_tuples(result.columns, names=["Gas", "Dataset"])
    return result


def manuscript_rounding(stats: pd.DataFrame) -> pd.DataFrame:
    """Apply the display precision used in the manuscript's Table 1."""

    rounded = stats.copy().astype(object)
    for column in rounded.columns:
        rounded.loc["Count", column] = int(stats.loc["Count", column])
        rounded.loc["Mean", column] = round(float(stats.loc["Mean", column]))
        rounded.loc["Median", column] = round(float(stats.loc["Median", column]), 2)
        minimum = float(stats.loc["Min", column])
        rounded.loc["Min", column] = round(minimum, 5 if minimum < 0.01 else 3)
        rounded.loc["Max", column] = round(float(stats.loc["Max", column]), 3)
    return rounded

