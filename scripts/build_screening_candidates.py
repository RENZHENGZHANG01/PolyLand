#!/usr/bin/env python3
"""Build the complete PolyLand post-novelty-filter screening table."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


GASES = ("O2", "N2", "H2", "CH4", "CO2")
SEPARATION_TASKS = (
    ("O2", "N2"),
    ("H2", "CH4"),
    ("H2", "N2"),
    ("CO2", "CH4"),
    ("CO2", "N2"),
)
SOURCE_LABELS = {
    "DiT_Feb162025": "DiT",
    "rule_based_allcombined_unique": "rule-based",
    "LLM_v2_combined_novel_smiles_drop_duplicate": "LLM",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing final_<gas>_results.csv for all five gases.",
    )
    parser.add_argument(
        "--selected",
        type=Path,
        required=True,
        help="CSV with selection_id and SMILES for the highlighted candidates.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/screening/polyland_screening_candidates.csv"),
        help="Destination CSV.",
    )
    return parser.parse_args()


def load_gas_table(path: Path, gas: str) -> pd.DataFrame:
    table = pd.read_csv(path)
    required = {"synthetic_data", "SMILES", "best_models_mean", "best_models_std"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    if table[["synthetic_data", "SMILES"]].duplicated().any():
        raise ValueError(f"{path} contains duplicate source/SMILES keys")
    if table[list(required)].isna().any().any():
        raise ValueError(f"{path} contains missing required values")

    return table[["synthetic_data", "SMILES", "best_models_mean", "best_models_std"]].rename(
        columns={
            "best_models_mean": f"{gas}_log10_Barrer_ensemble_mean",
            "best_models_std": f"{gas}_log10_Barrer_ensemble_sd",
        }
    )


def build_table(input_dir: Path, selected_path: Path) -> pd.DataFrame:
    gas_tables = {
        gas: load_gas_table(input_dir / f"final_{gas}_results.csv", gas)
        for gas in GASES
    }

    key_columns = ["synthetic_data", "SMILES"]
    reference_keys = gas_tables[GASES[0]][key_columns]
    for gas in GASES[1:]:
        if not reference_keys.equals(gas_tables[gas][key_columns]):
            raise ValueError(f"{gas} candidate keys or row order do not match {GASES[0]}")

    result = reference_keys.copy()
    result.insert(0, "candidate_id", [f"POLYLAND-{i:06d}" for i in range(1, len(result) + 1)])
    result.insert(1, "generation_method", result["synthetic_data"].map(SOURCE_LABELS))
    if result["generation_method"].isna().any():
        unknown = sorted(result.loc[result["generation_method"].isna(), "synthetic_data"].unique())
        raise ValueError(f"Unrecognized candidate sources: {unknown}")
    result = result.rename(columns={"synthetic_data": "source_table"})

    for gas in GASES:
        mean_col = f"{gas}_log10_Barrer_ensemble_mean"
        sd_col = f"{gas}_log10_Barrer_ensemble_sd"
        result[mean_col] = gas_tables[gas][mean_col].to_numpy()
        result[sd_col] = gas_tables[gas][sd_col].to_numpy()

    for numerator, denominator in SEPARATION_TASKS:
        prefix = f"{numerator}_{denominator}"
        log_selectivity = (
            result[f"{numerator}_log10_Barrer_ensemble_mean"]
            - result[f"{denominator}_log10_Barrer_ensemble_mean"]
        )
        result[f"{prefix}_selectivity"] = np.power(10.0, log_selectivity)

    selected = pd.read_csv(selected_path, dtype={"selection_id": "string", "SMILES": "string"})
    if set(selected.columns) != {"selection_id", "SMILES"}:
        raise ValueError("Selected-candidate CSV must contain exactly selection_id and SMILES")
    if selected["selection_id"].duplicated().any() or selected["SMILES"].duplicated().any():
        raise ValueError("Selected-candidate IDs and SMILES must be unique")

    selected_map = selected.set_index("SMILES")["selection_id"]
    result["selection_id"] = result["SMILES"].map(selected_map)
    result["selected"] = result["selection_id"].notna()
    matched = int(result["selected"].sum())
    if matched != len(selected):
        missing = selected.loc[~selected["SMILES"].isin(result["SMILES"]), "selection_id"].tolist()
        raise ValueError(
            f"Matched {matched} of {len(selected)} selected candidates; missing IDs: {missing}"
        )

    numeric_columns = result.select_dtypes(include=["number"]).columns
    if not np.isfinite(result[numeric_columns].to_numpy()).all():
        raise ValueError("Output contains non-finite numeric values")

    identity = [
        "candidate_id",
        "generation_method",
        "SMILES",
        "selected",
        "selection_id",
    ]
    result = result.drop(columns=["source_table"])
    remaining = [column for column in result.columns if column not in identity]
    return result[identity + remaining]


def main() -> None:
    args = parse_args()
    table = build_table(args.input_dir, args.selected)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, index=False, float_format="%.10g")

    source_counts = table["generation_method"].value_counts().to_dict()
    print(f"Wrote {len(table):,} candidates to {args.output}")
    print(f"Selected candidates: {int(table['selected'].sum())}")
    print(f"Generation methods: {source_counts}")


if __name__ == "__main__":
    main()
