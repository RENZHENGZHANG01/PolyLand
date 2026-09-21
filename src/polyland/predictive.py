"""Predictive-model utilities for the three PolyLand training settings."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .fingerprints import smiles_to_fingerprints


def load_split(
    processed_dir: str | Path,
    gas: str,
    training: str,
    test_size: float = 0.2,
    random_state: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Create a shared ladder holdout and one of the Li/La/Li+La train sets."""

    from sklearn.model_selection import train_test_split

    processed_dir = Path(processed_dir)
    linear = pd.read_csv(processed_dir / f"final_linear_data_{gas}.csv")
    ladder = pd.read_csv(processed_dir / f"final_ladder_data_{gas}.csv")
    ladder_train, ladder_test = train_test_split(
        ladder, test_size=test_size, random_state=random_state
    )
    if training == "linear":
        train = linear
    elif training == "ladder":
        train = ladder_train
    elif training == "hybrid":
        train = pd.concat([linear, ladder_train], ignore_index=True)
    else:
        raise ValueError("training must be one of: linear, ladder, hybrid")
    return train.reset_index(drop=True), ladder_test.reset_index(drop=True)


def _metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    from sklearn.metrics import mean_squared_error, r2_score

    return {
        "rmse": float(np.sqrt(mean_squared_error(y_true, y_pred))),
        "r2": float(r2_score(y_true, y_pred)),
    }


def train_fingerprint_model(
    train: pd.DataFrame,
    test: pd.DataFrame,
    gas: str,
    model_name: str,
    fingerprint: str = "Morgan",
    radius: int = 2,
    n_bits: int = 2048,
    random_state: int = 42,
) -> tuple[Any, pd.DataFrame, dict[str, float]]:
    """Train QRF or MLP-D and return model, holdout predictions, and metrics."""

    target = f"{gas}_log10"
    x_train = smiles_to_fingerprints(train["SMILES"].tolist(), fingerprint, radius, n_bits)
    x_test = smiles_to_fingerprints(test["SMILES"].tolist(), fingerprint, radius, n_bits)
    y_train = train[target].to_numpy(dtype=float)
    y_test = test[target].to_numpy(dtype=float)

    if model_name == "qrf":
        try:
            from quantile_forest import RandomForestQuantileRegressor
        except ImportError as exc:  # pragma: no cover
            raise RuntimeError("Install PolyLand with the 'ml' extra to train QRF") from exc
        model = RandomForestQuantileRegressor(random_state=random_state, n_jobs=-1)
        model.fit(x_train, y_train)
        mean = np.asarray(model.predict(x_test)).reshape(-1)
        low = np.asarray(model.predict(x_test, quantiles=0.05)).reshape(-1)
        high = np.asarray(model.predict(x_test, quantiles=0.95)).reshape(-1)
    elif model_name == "mlp":
        try:
            import tensorflow as tf
        except ImportError as exc:  # pragma: no cover
            raise RuntimeError("Install PolyLand with the 'mlp' extra to train MLP-D") from exc
        tf.keras.utils.set_random_seed(random_state)

        class MCDropout(tf.keras.layers.Dropout):
            def call(self, inputs: Any, training: bool | None = None) -> Any:
                return super().call(inputs, training=True)

        model = tf.keras.Sequential(
            [
                tf.keras.layers.Input(shape=(x_train.shape[1],)),
                tf.keras.layers.Dense(512, activation="relu"),
                MCDropout(0.2),
                tf.keras.layers.Dense(128, activation="relu"),
                MCDropout(0.2),
                tf.keras.layers.Dense(1),
            ]
        )
        model.compile(optimizer="adam", loss="mean_squared_error")
        model.fit(x_train, y_train, epochs=100, batch_size=32, validation_split=0.1)
        samples = np.stack(
            [model(x_test, training=True).numpy().reshape(-1) for _ in range(100)]
        )
        mean = samples.mean(axis=0)
        low, high = np.percentile(samples, [5, 95], axis=0)
    else:
        raise ValueError("model_name must be 'qrf' or 'mlp'")

    predictions = test[["PID", "SMILES"]].copy()
    predictions["observed_log10_Barrer"] = y_test
    predictions["predicted_log10_Barrer"] = mean
    predictions["lower_95"] = low
    predictions["upper_95"] = high
    return model, predictions, _metrics(y_test, mean)


def train_graph_model(
    train: pd.DataFrame,
    test: pd.DataFrame,
    gas: str,
    model_name: str,
    n_trials: int = 200,
    epochs: int = 500,
) -> tuple[Any, pd.DataFrame, dict[str, float]]:
    """Train the torch-molecule GIN/GCN/GREA models used in the paper."""

    try:
        from torch_molecule import GREAMolecularPredictor, GNNMolecularPredictor
        from torch_molecule.utils.search import ParameterSpec, ParameterType
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("Install PolyLand with the 'graph' extra for graph models") from exc

    target = f"{gas}_log10"
    y_train = train[target].to_numpy(dtype=float)
    y_test = test[target].to_numpy(dtype=float)
    smiles_train = train["SMILES"].tolist()
    smiles_test = test["SMILES"].tolist()

    common = dict(
        num_tasks=1,
        task_type="regression",
        model_name=f"{model_name}_{gas}",
        batch_size=512,
        epochs=epochs,
        evaluate_criterion="r2",
        evaluate_higher_better=True,
        verbose=True,
    )
    if model_name == "grea":
        model = GREAMolecularPredictor(**common)
        search = {
            "gnn_type": ParameterSpec(
                ParameterType.CATEGORICAL, ["gin-virtual", "gcn-virtual", "gin", "gcn"]
            ),
            "norm_layer": ParameterSpec(
                ParameterType.CATEGORICAL, ["batch_norm", "layer_norm", "size_norm"]
            ),
            "num_layer": ParameterSpec(ParameterType.INTEGER, value_range=(2, 6)),
            "emb_dim": ParameterSpec(ParameterType.INTEGER, value_range=(256, 512)),
            "learning_rate": ParameterSpec(ParameterType.FLOAT, value_range=(1e-4, 1e-2)),
            "drop_ratio": ParameterSpec(ParameterType.FLOAT, value_range=(0.05, 0.5)),
            "gamma": ParameterSpec(ParameterType.FLOAT, value_range=(0.25, 0.75)),
            "augmented_feature": ParameterSpec(
                ParameterType.CATEGORICAL, ["maccs,morgan", "maccs", "morgan", None]
            ),
        }
    elif model_name in {"gin", "gcn"}:
        model = GNNMolecularPredictor(**common)
        search = {
            "gnn_type": ParameterSpec(
                ParameterType.CATEGORICAL,
                [model_name, f"{model_name}-virtual"],
            ),
            "norm_layer": ParameterSpec(
                ParameterType.CATEGORICAL, ["batch_norm", "layer_norm", "size_norm"]
            ),
            "num_layer": ParameterSpec(ParameterType.INTEGER, value_range=(2, 5)),
            "emb_dim": ParameterSpec(ParameterType.INTEGER, value_range=(256, 512)),
            "learning_rate": ParameterSpec(ParameterType.FLOAT, value_range=(1e-4, 1e-2)),
            "drop_ratio": ParameterSpec(ParameterType.FLOAT, value_range=(0.05, 0.5)),
            "augmented_feature": ParameterSpec(
                ParameterType.CATEGORICAL, ["maccs,morgan", "maccs", "morgan", None]
            ),
        }
    else:
        raise ValueError("model_name must be one of: gin, gcn, grea")

    model.autofit(
        X_train=smiles_train,
        y_train=y_train,
        X_val=smiles_train,
        y_val=y_train,
        search_parameters=search,
        n_trials=n_trials,
    )
    output = model.predict(smiles_test)
    mean = np.asarray(output["prediction"]).reshape(-1)
    variance = np.asarray(output.get("variance", np.zeros_like(mean))).reshape(-1)
    interval = 1.96 * np.sqrt(np.maximum(variance, 0))
    predictions = test[["PID", "SMILES"]].copy()
    predictions["observed_log10_Barrer"] = y_test
    predictions["predicted_log10_Barrer"] = mean
    predictions["lower_95"] = mean - interval
    predictions["upper_95"] = mean + interval
    return model, predictions, _metrics(y_test, mean)


def save_run(
    model: Any,
    predictions: pd.DataFrame,
    metrics: dict[str, float],
    model_name: str,
    output_dir: str | Path,
) -> None:
    """Persist one trained model and its holdout results."""

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(output_dir / "holdout_predictions.csv", index=False)
    (output_dir / "metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if model_name in {"gin", "gcn", "grea"}:
        model.save_model(path=str(output_dir / "model.pt"))
    elif model_name == "mlp":
        model.save(output_dir / "model.keras")
    else:
        import joblib

        joblib.dump(model, output_dir / "model.pkl")

