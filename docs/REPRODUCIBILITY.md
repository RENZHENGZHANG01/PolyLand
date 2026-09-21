# Reproducibility map

| Manuscript component | Repository entry point | Status |
|---|---|---|
| Dataset curation / Table 1 | `scripts/prepare_data.py`, `scripts/build_table1.py` | Reproducible from included tables |
| QRF and MLP-D | `scripts/train_predictor.py` | Code and inputs included |
| GIN, GCN, GREA | `scripts/train_predictor.py` | Code and inputs included; installs `torch-molecule` |
| LLM ICL pairing/prompts | `scripts/prepare_icl.py` | Code and inputs included |
| LLM generation | `scripts/generate_llm.py` | Requires model access/API credentials |
| MD/FFV consistency check | `scripts/md_ffv_consistency.py` | Code and training inputs included; candidate FFV/predictions supplied by caller |
| Reaction-template generation | — | Archived project contained outputs, not source |
| Graph DiT training/generation | — | Archived project contained outputs, not source |

The two unavailable generators are stated explicitly so that this repository
does not imply end-to-end reproducibility where the source was not recoverable.

## Suggested Code Availability statement

> Code and curated input tables supporting the data preparation, Table 1,
> permeability models, LLM in-context-learning workflow, and MD/FFV
> consistency analysis are available at
> https://github.com/RENZHENGZHANG01/PolyLand. The repository documents the
> scope of the archived code, dependencies, data provenance, and the known
> reproducibility boundary for generator components whose source was not
> present in the project archive.

