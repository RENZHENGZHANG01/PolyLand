# PolyLand

Research code and datasets for **PolyLand: ML-Powered Design of Ladder
Polymers for Gas Separation**.

PolyLand provides a reproducible workflow for preparing gas-permeability
datasets, training polymer property predictors, constructing in-context
learning examples, and generating candidate ladder polymers. The repository
supports five gases: O2, N2, H2, CH4, and CO2.

## Features

- curated linear- and ladder-polymer permeability datasets;
- deterministic data preparation and summary-table generation;
- QRF, MLP-D, GIN, GCN, and GREA prediction workflows;
- Morgan, MACCS, and polyBERT molecular representations;
- optimization-aware in-context-learning example selection;
- OpenAI and Hugging Face generation entry points;
- automated data and regression tests.

## Repository layout

```text
data/
  raw/          Curated source tables and polymer metadata
  processed/    Per-gas modeling tables produced by the data pipeline
docs/           Workflow and reproducibility notes
results/        Small deterministic outputs
scripts/        Command-line entry points
src/polyland/   Reusable data, modeling, ICL, and analysis modules
tests/          Automated tests
```

## Installation

Python 3.10 or later is recommended.

```bash
git clone https://github.com/RENZHENGZHANG01/PolyLand.git
cd PolyLand
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

On Windows PowerShell, activate the environment with:

```powershell
.venv\Scripts\Activate.ps1
```

Install the optional dependencies required for a particular workflow:

```bash
# Molecular fingerprints and QRF
python -m pip install -e '.[ml]'

# TensorFlow MLP-D
python -m pip install -e '.[ml,mlp]'

# GIN, GCN, and GREA
python -m pip install -e '.[ml,graph]'

# LLM generation clients
python -m pip install -e '.[ml,llm]'

# Development and testing tools
python -m pip install -e '.[dev]'
```

## Quick start

Run the complete data preparation and summary-table workflow from the
repository root:

```bash
python scripts/prepare_data.py
python scripts/build_table1.py
python -m unittest discover -s tests -v
```

The prepared datasets are written to `data/processed/`, and the summary table
is written to `results/table1_stats.csv`.

## Data

Permeability values are reported in Barrer. Each processed gas table contains
the polymer identifier, repeat-unit SMILES, polymer class, molecular-dynamics
metadata, the measured permeability, and its base-10 logarithm.

The ten modeling tables follow this naming convention:

```text
data/processed/final_linear_data_<gas>.csv
data/processed/final_ladder_data_<gas>.csv
```

See [`data/README.md`](data/README.md) for column definitions and provenance.

## Train a permeability predictor

`train_predictor.py` supports three training-set configurations:

- `linear`: linear polymers only;
- `ladder`: ladder polymers only;
- `hybrid`: linear and ladder polymers.

Example QRF run for CO2 permeability:

```bash
python scripts/train_predictor.py \
  --gas CO2 \
  --training hybrid \
  --model qrf \
  --fingerprint Morgan
```

Example GREA run:

```bash
python scripts/train_predictor.py \
  --gas CO2 \
  --training hybrid \
  --model grea \
  --n-trials 200
```

Supported gases are `O2`, `N2`, `H2`, `CH4`, and `CO2`. Permeability targets
are modeled as `log10(Barrer)`. Model artifacts and holdout predictions are
written below `results/models/`.

## Prepare ICL examples and generate candidates

Create similarity-ranked in-context-learning examples:

```bash
python scripts/prepare_icl.py --k 5 --seed 42
```

Generate candidates with an OpenAI model:

```bash
export OPENAI_API_KEY=your_key_here
python scripts/generate_llm.py \
  --prompts results/icl/CO2_CH4_filtered_sim_k5.json \
  --provider openai \
  --model gpt-4o-mini \
  --output results/llm/co2_ch4.csv
```

To use a Hugging Face model, set `--provider huggingface` and provide the model
ID with `--model`. Some gated models require `HUGGINGFACE_HUB_TOKEN`.
Credentials are read from the environment and should never be committed.

## Testing

Run the test suite with:

```bash
python -m unittest discover -s tests -v
```

GitHub Actions runs the data pipeline and tests for every push and pull
request to `main`.

## Citation

If you use PolyLand in your work, cite the accompanying manuscript. Citation
metadata are available in [`CITATION.cff`](CITATION.cff).

## License

The code is released under the MIT License. See [`NOTICE`](NOTICE) for POINT2
attribution and information about literature-derived data.
