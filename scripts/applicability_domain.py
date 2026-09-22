#!/usr/bin/env python3
"""Structure-based applicability-domain analysis for screened candidates.

The analysis uses Morgan fingerprints (radius 2, 2048 bits) and Tanimoto
similarity.  For each gas, it treats the original ladder-only and hybrid
(linear + ladder) training pools as separate reference domains.  A domain
threshold is derived from the fifth percentile of leave-one-out nearest-
neighbor similarities within the corresponding training pool.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit import RDLogger


GASES = ("O2", "N2", "H2", "CH4", "CO2")
SOURCE_LABELS = {
    "DiT": "DiT",
    "LLM": "LLM",
    "rule-based": "Rule-based",
}


@dataclass(frozen=True)
class ReferenceDomain:
    gas: str
    name: str
    fingerprints: list
    pids: list[str]
    polymer_types: list[str]
    threshold: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--candidates",
        type=Path,
        default=Path("data/screening/polyland_screening_candidates.csv"),
    )
    parser.add_argument(
        "--linear-dir",
        type=Path,
        default=Path("data/processed"),
        help="Directory containing final_linear_data_<gas>.csv files.",
    )
    parser.add_argument(
        "--ladder-split-dir",
        type=Path,
        default=None,
        help="Directory containing final_ladder_data_<gas>_train.csv files.",
    )
    parser.add_argument(
        "--ladder-training-structures",
        type=Path,
        default=Path("data/processed/ad_ladder_training_structures.csv"),
        help=(
            "Gas-specific ladder structures used for training. Used when "
            "--ladder-split-dir is omitted."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/applicability_domain"),
    )
    parser.add_argument("--radius", type=int, default=2)
    parser.add_argument("--n-bits", type=int, default=2048)
    parser.add_argument("--threshold-quantile", type=float, default=0.05)
    return parser.parse_args()


def canonicalize(smiles: str) -> tuple[str | None, object | None]:
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None, None
    return Chem.MolToSmiles(mol, canonical=True), mol


def prepare_unique_reference(
    frames: Iterable[pd.DataFrame], generator
) -> tuple[list, list[str], list[str], list[str]]:
    seen: set[str] = set()
    fingerprints: list = []
    pids: list[str] = []
    polymer_types: list[str] = []
    canonical_smiles: list[str] = []

    for frame in frames:
        for row in frame.itertuples(index=False):
            canonical, mol = canonicalize(getattr(row, "SMILES"))
            if canonical is None or canonical in seen:
                continue
            seen.add(canonical)
            fingerprints.append(generator.GetFingerprint(mol))
            pids.append(str(getattr(row, "PID")))
            polymer_types.append(str(getattr(row, "Type")))
            canonical_smiles.append(canonical)
    return fingerprints, pids, polymer_types, canonical_smiles


def leave_one_out_max_similarity(fingerprints: list) -> np.ndarray:
    if len(fingerprints) < 2:
        raise ValueError("At least two unique training structures are required")
    values = np.empty(len(fingerprints), dtype=float)
    for i, fp in enumerate(fingerprints):
        similarities = DataStructs.BulkTanimotoSimilarity(fp, fingerprints)
        similarities[i] = -1.0
        values[i] = max(similarities)
    return values


def build_domains(
    gas: str,
    linear_dir: Path,
    ladder_split_dir: Path | None,
    ladder_training_structures: Path,
    generator,
    threshold_quantile: float,
) -> tuple[ReferenceDomain, ReferenceDomain, list[dict]]:
    if ladder_split_dir is not None:
        ladder = pd.read_csv(ladder_split_dir / f"final_ladder_data_{gas}_train.csv")
    else:
        ladder_all = pd.read_csv(ladder_training_structures)
        ladder = ladder_all[ladder_all["gas"].eq(gas)].copy()
    linear = pd.read_csv(linear_dir / f"final_linear_data_{gas}.csv")

    ladder_values = prepare_unique_reference([ladder], generator)
    hybrid_values = prepare_unique_reference([ladder, linear], generator)
    records: list[dict] = []
    domains: list[ReferenceDomain] = []

    for name, values in (("ladder", ladder_values), ("hybrid", hybrid_values)):
        fingerprints, pids, polymer_types, _ = values
        loo = leave_one_out_max_similarity(fingerprints)
        threshold = float(np.quantile(loo, threshold_quantile, method="linear"))
        domains.append(
            ReferenceDomain(
                gas=gas,
                name=name,
                fingerprints=fingerprints,
                pids=pids,
                polymer_types=polymer_types,
                threshold=threshold,
            )
        )
        records.append(
            {
                "gas": gas,
                "domain": name,
                "n_unique_training_structures": len(fingerprints),
                "threshold_quantile": threshold_quantile,
                "similarity_threshold": threshold,
                "training_loo_similarity_min": float(np.min(loo)),
                "training_loo_similarity_median": float(np.median(loo)),
                "training_loo_similarity_max": float(np.max(loo)),
            }
        )
    return domains[0], domains[1], records


def fingerprint_candidates(
    smiles_values: pd.Series, generator
) -> tuple[list, np.ndarray, list[str | None], list[str]]:
    canonical_to_index: dict[str, int] = {}
    unique_fingerprints: list = []
    row_to_unique = np.full(len(smiles_values), -1, dtype=int)
    canonical_rows: list[str | None] = []
    invalid: list[str] = []

    for row_index, smiles in enumerate(smiles_values.astype(str)):
        canonical, mol = canonicalize(smiles)
        canonical_rows.append(canonical)
        if canonical is None:
            invalid.append(smiles)
            continue
        unique_index = canonical_to_index.get(canonical)
        if unique_index is None:
            unique_index = len(unique_fingerprints)
            canonical_to_index[canonical] = unique_index
            unique_fingerprints.append(generator.GetFingerprint(mol))
        row_to_unique[row_index] = unique_index
    return unique_fingerprints, row_to_unique, canonical_rows, invalid


def nearest_neighbor_scores(
    candidate_fingerprints: list, domain: ReferenceDomain
) -> tuple[np.ndarray, np.ndarray]:
    similarities = np.empty(len(candidate_fingerprints), dtype=float)
    nearest_indices = np.empty(len(candidate_fingerprints), dtype=int)
    for i, fp in enumerate(candidate_fingerprints):
        values = DataStructs.BulkTanimotoSimilarity(fp, domain.fingerprints)
        nearest_index = int(np.argmax(values))
        similarities[i] = values[nearest_index]
        nearest_indices[i] = nearest_index
        if (i + 1) % 10000 == 0:
            print(
                f"{domain.gas} {domain.name}: scored {i + 1:,}/"
                f"{len(candidate_fingerprints):,} unique candidates",
                flush=True,
            )
    return similarities, nearest_indices


def expand_scores_to_rows(
    output: pd.DataFrame,
    row_to_unique: np.ndarray,
    domain: ReferenceDomain,
    unique_similarities: np.ndarray,
    unique_neighbors: np.ndarray,
) -> None:
    prefix = f"{domain.gas}_{domain.name}"
    valid = row_to_unique >= 0
    row_similarity = np.full(len(output), np.nan, dtype=float)
    row_neighbor = np.full(len(output), -1, dtype=int)
    row_similarity[valid] = unique_similarities[row_to_unique[valid]]
    row_neighbor[valid] = unique_neighbors[row_to_unique[valid]]

    output[f"{prefix}_max_tanimoto"] = row_similarity
    output[f"{prefix}_nearest_pid"] = [
        domain.pids[i] if i >= 0 else "" for i in row_neighbor
    ]
    output[f"{prefix}_nearest_polymer_type"] = [
        domain.polymer_types[i] if i >= 0 else "" for i in row_neighbor
    ]
    output[f"{prefix}_in_domain"] = valid & (row_similarity >= domain.threshold)


def summarize(output: pd.DataFrame, thresholds: pd.DataFrame) -> pd.DataFrame:
    group_masks: list[tuple[str, pd.Series]] = [
        ("All candidates", pd.Series(True, index=output.index)),
        ("Selected candidates", output["selected"].astype(bool)),
    ]
    for source, label in SOURCE_LABELS.items():
        group_masks.append((label, output["generation_method"].eq(source)))

    records: list[dict] = []
    for threshold in thresholds.itertuples(index=False):
        prefix = f"{threshold.gas}_{threshold.domain}"
        score_column = f"{prefix}_max_tanimoto"
        flag_column = f"{prefix}_in_domain"
        for label, mask in group_masks:
            scores = output.loc[mask, score_column].dropna()
            flags = output.loc[mask & output[score_column].notna(), flag_column]
            records.append(
                {
                    "gas": threshold.gas,
                    "domain": threshold.domain,
                    "candidate_group": label,
                    "n_candidates": int(mask.sum()),
                    "n_valid_smiles": int(len(scores)),
                    "n_in_domain": int(flags.sum()),
                    "fraction_in_domain": float(flags.mean()) if len(flags) else np.nan,
                    "max_tanimoto_q05": float(scores.quantile(0.05)) if len(scores) else np.nan,
                    "max_tanimoto_median": float(scores.median()) if len(scores) else np.nan,
                    "max_tanimoto_q95": float(scores.quantile(0.95)) if len(scores) else np.nan,
                    "similarity_threshold": float(threshold.similarity_threshold),
                }
            )
    return pd.DataFrame.from_records(records)


def plot_summary(summary: pd.DataFrame, output_dir: Path) -> None:
    row_order = ["All candidates", "DiT", "LLM", "Rule-based", "Selected candidates"]
    gas_order = list(GASES)
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.6), constrained_layout=True)
    image = None

    for ax, domain, title in zip(
        axes,
        ("ladder", "hybrid"),
        ("Ladder-only training domain", "Hybrid training domain"),
    ):
        subset = summary[summary["domain"].eq(domain)]
        matrix = (
            subset.pivot(
                index="candidate_group", columns="gas", values="fraction_in_domain"
            )
            .reindex(index=row_order, columns=gas_order)
            .to_numpy()
        )
        image = ax.imshow(matrix * 100.0, vmin=0, vmax=100, cmap="YlGnBu", aspect="auto")
        ax.set_xticks(range(len(gas_order)), gas_order)
        ax.set_yticks(range(len(row_order)), row_order)
        ax.set_title(title, fontsize=11, pad=10)
        ax.set_xlabel("Gas permeability model")
        for row in range(matrix.shape[0]):
            for col in range(matrix.shape[1]):
                value = matrix[row, col] * 100.0
                color = "white" if value >= 62 else "black"
                ax.text(col, row, f"{value:.1f}%", ha="center", va="center", color=color, fontsize=8.5)
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)

    assert image is not None
    colorbar = fig.colorbar(image, ax=axes, fraction=0.035, pad=0.03)
    colorbar.set_label("Candidates inside applicability domain (%)")
    fig.savefig(output_dir / "applicability_domain_summary.pdf", bbox_inches="tight")
    fig.savefig(output_dir / "applicability_domain_summary.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    if not 0.0 < args.threshold_quantile < 1.0:
        raise ValueError("--threshold-quantile must be between 0 and 1")

    RDLogger.DisableLog("rdApp.warning")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=args.radius, fpSize=args.n_bits
    )

    candidates = pd.read_csv(args.candidates)
    required = {"candidate_id", "generation_method", "SMILES", "selected", "selection_id"}
    missing = required.difference(candidates.columns)
    if missing:
        raise ValueError(f"Candidate table is missing required columns: {sorted(missing)}")

    output = candidates[
        ["candidate_id", "generation_method", "SMILES", "selected", "selection_id"]
    ].copy()
    if output["selected"].dtype != bool:
        output["selected"] = (
            output["selected"].astype(str).str.strip().str.lower().eq("true")
        )

    candidate_fps, row_to_unique, canonical_rows, invalid = fingerprint_candidates(
        output["SMILES"], generator
    )
    output.insert(3, "canonical_SMILES", canonical_rows)
    output.insert(4, "valid_smiles", row_to_unique >= 0)
    print(
        f"Candidates: {len(output):,} rows; {len(candidate_fps):,} unique canonical "
        f"structures; {len(invalid):,} invalid SMILES",
        flush=True,
    )

    threshold_records: list[dict] = []
    for gas in GASES:
        ladder_domain, hybrid_domain, records = build_domains(
            gas,
            args.linear_dir,
            args.ladder_split_dir,
            args.ladder_training_structures,
            generator,
            args.threshold_quantile,
        )
        threshold_records.extend(records)
        for domain in (ladder_domain, hybrid_domain):
            unique_scores, unique_neighbors = nearest_neighbor_scores(candidate_fps, domain)
            expand_scores_to_rows(
                output,
                row_to_unique,
                domain,
                unique_scores,
                unique_neighbors,
            )

    thresholds = pd.DataFrame.from_records(threshold_records)
    thresholds.insert(5, "fingerprint", "Morgan")
    thresholds.insert(6, "radius", args.radius)
    thresholds.insert(7, "n_bits", args.n_bits)
    summary = summarize(output, thresholds)
    selected = output[output["selected"]].copy()

    output.to_csv(
        args.output_dir / "candidate_applicability_domain.csv.gz",
        index=False,
        compression="gzip",
    )
    selected.to_csv(args.output_dir / "selected_candidate_applicability_domain.csv", index=False)
    thresholds.to_csv(args.output_dir / "applicability_domain_thresholds.csv", index=False)
    summary.to_csv(args.output_dir / "applicability_domain_summary.csv", index=False)
    plot_summary(summary, args.output_dir)

    print("\nThresholds")
    print(thresholds[["gas", "domain", "n_unique_training_structures", "similarity_threshold"]].to_string(index=False))
    print("\nSelected candidate coverage")
    print(
        summary[summary["candidate_group"].eq("Selected candidates")][
            ["gas", "domain", "n_in_domain", "n_valid_smiles", "fraction_in_domain"]
        ].to_string(index=False)
    )
    print(f"\nResults written to {args.output_dir}")


if __name__ == "__main__":
    main()
