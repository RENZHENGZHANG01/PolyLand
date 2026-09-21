#!/usr/bin/env python3
"""Analyze empirical FFV/permeability relationships for candidate polymers."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from polyland.md_ffv import consistency_table


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--gas", choices=["O2", "N2", "H2", "CH4", "CO2"], required=True)
    parser.add_argument("--prediction-column", required=True)
    parser.add_argument("--degree", type=int, choices=[1, 2, 3], default=2)
    parser.add_argument("--processed-dir", type=Path, default=Path("data/processed"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    candidates = pd.read_csv(args.candidates)
    result = consistency_table(
        args.processed_dir,
        candidates,
        args.gas,
        args.degree,
        args.prediction_column,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    print(f"Wrote {len(result)} FFV-analysis rows to {args.output}")


if __name__ == "__main__":
    main()
