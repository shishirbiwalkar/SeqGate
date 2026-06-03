# Contributing

## Local setup

```bash
# 1. Install tools (conda recommended)
conda create -n seqgate python=3.10
conda activate seqgate
pip install -r requirements.txt

# 2. Install bioinformatics tools (system-level)
# Required: FastQC, MultiQC, BWA, STAR, GATK, featureCounts (subread)
# On macOS via Homebrew/bioconda; on Linux via bioconda/apt

# 3. Configure paths
cp config/pipeline_config.yaml config/pipeline_config.local.yaml
# Edit pipeline_config.local.yaml — set reference genome paths, tool paths

# 4. Download test data
bash scripts/download_data.sh
```

## Running the pipelines

```bash
# DNA pipeline (QC → BWA alignment → GATK variant calling)
python scripts/dna_pipeline.py

# RNA pipeline (QC → STAR alignment → featureCounts)
python scripts/rna_pipeline.py

# Visualise results
python scripts/plot_results.py
```

## Adding a pipeline step

1. Add the function to the relevant pipeline script (`dna_pipeline.py` or `rna_pipeline.py`).
2. Any configurable parameter (paths, thresholds, tool flags) goes in `config/pipeline_config.yaml`, not hardcoded.
3. Add an entry to the `docs/` dataset table if you add new test data.
4. Update `README.md` if the step changes the expected workflow.

## Config

`config/pipeline_config.yaml` is the single source of truth for all paths and parameters. Never hardcode absolute paths in scripts — read them from config. The `.local.yaml` variant is git-ignored so personal machine paths don't leak.

## Code style

- Python scripts: PEP 8, descriptive variable names.
- Shell scripts: `set -euo pipefail` at the top.
- Print progress to stdout with clear section headers (see existing scripts for pattern).

## Test data

Minimal test datasets live in `data/` after running `scripts/download_data.sh`. These are subsets designed to complete quickly (seconds, not hours). Do not commit raw FASTQ files.
