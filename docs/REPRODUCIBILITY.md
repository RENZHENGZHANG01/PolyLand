# Workflow reference

The following entry points cover the main PolyLand workflows.

| Workflow | Entry point | Primary output |
|---|---|---|
| Prepare modeling tables | `scripts/prepare_data.py` | `data/processed/*.csv` |
| Build dataset summary | `scripts/build_table1.py` | `results/table1_stats.csv` |
| Train permeability models | `scripts/train_predictor.py` | `results/models/` |
| Prepare ICL examples | `scripts/prepare_icl.py` | `results/icl/` |
| Generate candidates | `scripts/generate_llm.py` | User-selected CSV output |
| Analyze FFV relationships | `scripts/md_ffv_consistency.py` | User-selected CSV output |

Run all commands from the repository root. Random seeds and data-selection
rules are exposed through the command-line interfaces so experiments can be
repeated consistently.

Exact numerical agreement for neural and graph models can depend on hardware,
dependency versions, and hyperparameter-search execution. The deterministic
data pipeline and regression tests are also run by GitHub Actions.
