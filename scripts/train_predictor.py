#!/usr/bin/env python3
"""Train one PolyLand permeability model on a shared ladder holdout."""

from __future__ import annotations

import argparse
from pathlib import Path

from polyland.predictive import (
    load_split,
    save_run,
    train_fingerprint_model,
    train_graph_model,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gas", choices=["O2", "N2", "H2", "CH4", "CO2"], required=True)
    parser.add_argument("--training", choices=["linear", "ladder", "hybrid"], required=True)
    parser.add_argument("--model", choices=["qrf", "mlp", "gin", "gcn", "grea"], required=True)
    parser.add_argument(
        "--fingerprint",
        choices=["Morgan", "RDKit", "MACCS", "TopologicalTorsion", "AtomPair"],
        default="Morgan",
    )
    parser.add_argument("--radius", type=int, default=2)
    parser.add_argument("--n-bits", type=int, default=2048)
    parser.add_argument("--n-trials", type=int, default=200)
    parser.add_argument("--epochs", type=int, default=500)
    parser.add_argument("--processed-dir", type=Path, default=Path("data/processed"))
    parser.add_argument("--output-root", type=Path, default=Path("results/models"))
    args = parser.parse_args()

    train, test = load_split(args.processed_dir, args.gas, args.training)
    if args.model in {"qrf", "mlp"}:
        model, predictions, metrics = train_fingerprint_model(
            train,
            test,
            args.gas,
            args.model,
            args.fingerprint,
            args.radius,
            args.n_bits,
        )
    else:
        model, predictions, metrics = train_graph_model(
            train, test, args.gas, args.model, args.n_trials, args.epochs
        )
    output = args.output_root / args.gas / args.training / f"{args.model}_{args.fingerprint}"
    save_run(model, predictions, metrics, args.model, output)
    print(metrics)
    print(f"Run written to {output}")


if __name__ == "__main__":
    main()

