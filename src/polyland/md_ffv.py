"""Empirical FFV/permeability models used as a consistency check."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

GAS_THRESHOLDS = {
    "CH4": (3.0, 5.4),
    "CO2": (4.0, 6.0),
    "H2": (4.0, 6.0),
    "N2": (3.0, 6.0),
    "O2": (4.0, 6.0),
}


def load_empirical_data(processed_dir: str | Path, gas: str) -> pd.DataFrame:
    """Load linear and ladder records and apply the paper's FFV outlier limits."""

    processed_dir = Path(processed_dir)
    frames = []
    for kind in ("linear", "ladder"):
        frame = pd.read_csv(processed_dir / f"final_{kind}_data_{gas}.csv")
        frame["Dataset"] = kind
        frames.append(frame)
    data = pd.concat(frames, ignore_index=True)
    data = data.loc[data["FFV"].notna() & (data["FFV"] > 0)].copy()
    data["inverse_FFV"] = 1.0 / data["FFV"]
    max_log, max_inverse = GAS_THRESHOLDS[gas]
    return data.loc[
        (data[f"{gas}_log10"] <= max_log) & (data["inverse_FFV"] <= max_inverse)
    ].reset_index(drop=True)


def fit_polynomial_bayesian(data: pd.DataFrame, gas: str, degree: int):
    """Fit Bayesian ridge regression of log10(permeability) on 1/FFV."""

    if degree not in {1, 2, 3}:
        raise ValueError("degree must be 1, 2, or 3")
    from sklearn.linear_model import BayesianRidge
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import PolynomialFeatures, StandardScaler

    model = make_pipeline(
        PolynomialFeatures(degree=degree, include_bias=False),
        StandardScaler(),
        BayesianRidge(),
    )
    model.fit(data[["inverse_FFV"]], data[f"{gas}_log10"])
    return model


def consistency_table(
    processed_dir: str | Path,
    candidates: pd.DataFrame,
    gas: str,
    degree: int,
    prediction_column: str,
) -> pd.DataFrame:
    """Compare ML predictions with an empirical FFV trend.

    The returned residual is a diagnostic agreement measure.  It is not an
    independent validation because both quantities are model-derived.
    """

    required = {"SMILES", "FFV", prediction_column}
    missing = required.difference(candidates.columns)
    if missing:
        raise ValueError(f"Candidate table is missing columns: {sorted(missing)}")
    training = load_empirical_data(processed_dir, gas)
    model = fit_polynomial_bayesian(training, gas, degree)
    result = candidates.copy()
    result["inverse_FFV"] = 1.0 / result["FFV"].astype(float)
    result[f"{gas}_empirical_log10"] = model.predict(result[["inverse_FFV"]])
    result[f"{gas}_consistency_residual"] = (
        result[prediction_column] - result[f"{gas}_empirical_log10"]
    )
    return result

