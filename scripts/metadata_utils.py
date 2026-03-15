"""
Metadata utilities shared by DNA and RNA workflows.

Keeps manifest schema checks in one place so multi-omics projects stay consistent.
"""

from pathlib import Path
from typing import Iterable

import pandas as pd


def _require_columns(df: pd.DataFrame, required: Iterable[str], context: str) -> None:
    missing = set(required).difference(df.columns)
    if missing:
        raise ValueError(f"{context} manifest missing required columns: {sorted(missing)}")


def load_dna_manifest(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    _require_columns(df, {"sample_id", "fastq_r1", "fastq_r2"}, "DNA")
    return df


def load_rna_manifest(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    _require_columns(df, {"sample_id", "fastq_r1", "fastq_r2"}, "RNA")
    return df


def join_dna_rna_manifests(dna_df: pd.DataFrame, rna_df: pd.DataFrame) -> pd.DataFrame:
    #align DNA and RNA samples on sample_id to support multi-omics analyses
    merged = dna_df.merge(rna_df, on="sample_id", suffixes=("_dna", "_rna"))
    return merged

