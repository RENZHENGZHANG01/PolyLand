"""Optimization-aware in-context-learning data preparation."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from .data import GASES

SEPARATION_TASKS = (("O2", "N2"), ("H2", "N2"), ("H2", "CH4"), ("CO2", "N2"), ("CO2", "CH4"))


def _fingerprint(smiles: str):
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("Install PolyLand with the 'ml' extra to prepare ICL pairs") from exc
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)


def build_pairs(
    processed_dir: str | Path,
    target_gas: str,
    competing_gas: str,
    k: int = 5,
    optimized_only: bool = True,
    rank_by_similarity: bool = True,
    random_state: int = 42,
) -> pd.DataFrame:
    """Pair each source polymer with up to *k* candidate demonstrations."""

    from rdkit import DataStructs

    processed_dir = Path(processed_dir)
    target = pd.read_csv(processed_dir / f"final_ladder_data_{target_gas}.csv")
    competitor = pd.read_csv(processed_dir / f"final_ladder_data_{competing_gas}.csv")
    paired = target[["PID", "SMILES", f"{target_gas}_log10"]].merge(
        competitor[["PID", f"{competing_gas}_log10"]], on="PID", how="inner"
    )
    fps = {s: _fingerprint(s) for s in paired["SMILES"].unique()}
    rows: list[dict[str, object]] = []
    for _, source in paired.iterrows():
        candidates = paired.loc[paired["SMILES"] != source["SMILES"]].copy()
        if optimized_only:
            candidates = candidates.loc[
                (candidates[f"{target_gas}_log10"] > source[f"{target_gas}_log10"])
                & (candidates[f"{competing_gas}_log10"] < source[f"{competing_gas}_log10"])
            ]
        if candidates.empty or fps.get(source["SMILES"]) is None:
            continue
        candidates["Similarity"] = candidates["SMILES"].map(
            lambda value: DataStructs.TanimotoSimilarity(
                fps[source["SMILES"]], fps[value]
            )
            if fps.get(value) is not None
            else 0.0
        )
        if rank_by_similarity:
            candidates = candidates.nlargest(k, "Similarity")
        else:
            candidates = candidates.sample(
                n=min(k, len(candidates)), random_state=random_state
            )
        for _, candidate in candidates.iterrows():
            rows.append(
                {
                    "PID": source["PID"],
                    "Original_SMILES": source["SMILES"],
                    f"Original_{target_gas}_log10": source[f"{target_gas}_log10"],
                    f"Original_{competing_gas}_log10": source[f"{competing_gas}_log10"],
                    "Optimized_SMILES": candidate["SMILES"],
                    f"Optimized_{target_gas}_log10": candidate[f"{target_gas}_log10"],
                    f"Optimized_{competing_gas}_log10": candidate[f"{competing_gas}_log10"],
                    "Similarity": candidate["Similarity"],
                }
            )
    return pd.DataFrame(rows)


def prompt_record(
    smiles: str,
    target_gas: str,
    competing_gas: str,
    examples: pd.DataFrame,
) -> dict[str, object]:
    """Create one structured ICL prompt record."""

    return {
        "context": (
            "Ladder polymers are rigid, double-stranded macromolecules. "
            f"Optimize a repeat unit for {target_gas}/{competing_gas} separation. "
            "An asterisk marks a polymerization point."
        ),
        "examples": [
            {
                "SMILES_original": row["Original_SMILES"],
                "SMILES_optimized": row["Optimized_SMILES"],
            }
            for _, row in examples.iterrows()
        ],
        "question": {
            "instruction": (
                f"Increase {target_gas} permeability while reducing {competing_gas} "
                "permeability. Retain a ladder or semi-ladder repeat unit and exactly "
                "2 or 4 polymerization-point asterisks. Return a new, valid SMILES."
            ),
            "SMILES_test": smiles,
        },
        "expected_output_format": "Optimized SMILES: {SMILES}",
    }


def write_prompts(records: dict[str, dict[str, object]], path: str | Path) -> None:
    Path(path).write_text(json.dumps(records, indent=2) + "\n", encoding="utf-8")

