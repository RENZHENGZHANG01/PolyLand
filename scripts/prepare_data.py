#!/usr/bin/env python3
"""Rebuild processed PolyLand permeability tables from curated sources."""

from __future__ import annotations

import argparse
from pathlib import Path

from polyland.data import build_processed_tables


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, default=Path("data/raw"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed"))
    args = parser.parse_args()
    for name, path in sorted(build_processed_tables(args.raw_dir, args.output_dir).items()):
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()

