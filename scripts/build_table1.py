#!/usr/bin/env python3
"""Recompute Table 1 from the processed per-gas tables."""

from __future__ import annotations

import argparse
from pathlib import Path

from polyland.table1 import manuscript_rounding, table1_statistics


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", type=Path, default=Path("data/processed"))
    parser.add_argument("--output", type=Path, default=Path("results/table1_stats.csv"))
    args = parser.parse_args()
    stats = table1_statistics(args.processed_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    stats.to_csv(args.output)
    print(manuscript_rounding(stats).to_string())
    print(f"\nFull-precision statistics written to {args.output}")


if __name__ == "__main__":
    main()

