# Workflow reference

The following entry points cover the main PolyLand workflows.

| Workflow | Entry point | Primary output |
|---|---|---|
| Prepare modeling tables | `scripts/prepare_data.py` | `data/processed/*.csv` |
| Build dataset summary | `scripts/build_table1.py` | `results/table1_stats.csv` |
| Train permeability models | `scripts/train_predictor.py` | `results/models/` |
| Prepare ICL examples | `scripts/prepare_icl.py` | `results/icl/` |
| Generate candidates | `scripts/generate_llm.py` | User-selected CSV output |
| Build the full screening table | `scripts/build_screening_candidates.py` | `data/screening/polyland_screening_candidates.csv` |
| Analyze FFV relationships | `scripts/md_ffv_consistency.py` | User-selected CSV output |

Run all commands from the repository root. Random seeds and data-selection
rules are exposed through the command-line interfaces so experiments can be
repeated consistently.

Exact numerical agreement for neural and graph models can depend on hardware,
dependency versions, and hyperparameter-search execution. The deterministic
data pipeline and regression tests are also run by GitHub Actions.

The screening-table builder expects `final_<gas>_results.csv` for O2, N2, H2,
CH4, and CO2. It verifies that all five files contain the same ordered
source/SMILES keys, that the keys are unique, and that every selected SMILES
matches exactly one screening row before writing the combined CSV.
