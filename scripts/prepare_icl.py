#!/usr/bin/env python3
"""Build optimization-aware ICL pairs and structured prompt files."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from polyland.icl import SEPARATION_TASKS, build_pairs, prompt_record, write_prompts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", type=Path, default=Path("data/processed"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/icl"))
    parser.add_argument("--k", type=int, choices=[1, 5, 10], default=5)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for target, competitor in SEPARATION_TASKS:
        pairs = build_pairs(
            args.processed_dir,
            target,
            competitor,
            k=args.k,
            optimized_only=True,
            rank_by_similarity=True,
            random_state=args.seed,
        )
        stem = f"{target}_{competitor}_filtered_sim_k{args.k}"
        pairs.to_csv(args.output_dir / f"{stem}.csv", index=False)
        records: dict[str, dict[str, object]] = {}
        for pid, group in pairs.groupby("PID", sort=True):
            first = group.iloc[0]
            records[str(pid)] = prompt_record(
                str(first["Original_SMILES"]), target, competitor, group.head(args.k)
            )
        write_prompts(records, args.output_dir / f"{stem}.json")
        print(f"{stem}: {len(records)} prompts")


if __name__ == "__main__":
    main()

