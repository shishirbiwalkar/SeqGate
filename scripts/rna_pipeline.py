#!/usr/bin/env python
"""
RNA workflow orchestrator for Genomics_Workflows.

Transforms raw FASTQ into:
  - QC reports (FastQC/MultiQC)
  - STAR-aligned BAM files
  - Gene-level count matrix via featureCounts
"""

import argparse
import subprocess
import textwrap
from pathlib import Path

import yaml
import pandas as pd

from metadata_utils import load_rna_manifest


def run(cmd: str) -> None:
    #thin wrapper to print and run shell commands
    print(f"\n[cmd] {cmd}")
    subprocess.run(cmd, shell=True, check=True)


def load_config(config_path: Path) -> dict:
    with config_path.open() as fh:
        return yaml.safe_load(fh)


def prepare_output_dirs(output_root: Path) -> None:
    for sub in ["fastqc", "bam", "counts", "logs", "tmp"]:
        (output_root / sub).mkdir(parents=True, exist_ok=True)


def build_fastqc_cmd(fastq_r1: Path, fastq_r2: Path, out_dir: Path, threads: int) -> str:
    return (
        f"fastqc --threads {threads} "
        f"--outdir {out_dir} "
        f"{fastq_r1} {fastq_r2}"
    )


def build_star_cmd(sample: str, fastq_r1: Path, fastq_r2: Path, genome_dir: Path, out_prefix: Path, threads: int) -> str:
    cmd = textwrap.dedent(
        f"""
        STAR \\
          --runThreadN {threads} \\
          --genomeDir {genome_dir} \\
          --readFilesIn {fastq_r1} {fastq_r2} \\
          --readFilesCommand zcat \\
          --outSAMtype BAM SortedByCoordinate \\
          --outFileNamePrefix {out_prefix}
        """
    ).strip()
    return cmd


def build_featurecounts_cmd(bam_paths, gtf: Path, out_file: Path, threads: int, stranded: int = 0) -> str:
    bam_list = " ".join(str(p) for p in bam_paths)
    cmd = textwrap.dedent(
        f"""
        featureCounts \\
          -T {threads} \\
          -p \\
          -s {stranded} \\
          -a {gtf} \\
          -o {out_file} \\
          {bam_list}
        """
    ).strip()
    return cmd


def rna_workflow(config_path: Path) -> None:
    cfg = load_config(config_path)
    input_cfg = cfg["input"]
    ref_cfg = cfg["reference"]
    qc_cfg = cfg["qc"]
    align_cfg = cfg["alignment"]

    manifest_path = Path(input_cfg["rna_manifest"])
    output_root = Path(input_cfg["output_dir"]).joinpath("rna")

    prepare_output_dirs(output_root)
    manifest = load_rna_manifest(manifest_path)

    bam_paths = []

    for _, row in manifest.iterrows():
        sample = str(row["sample_id"])
        fastq_r1 = Path(row["fastq_r1"])
        fastq_r2 = Path(row["fastq_r2"])

        fastqc_dir = output_root.joinpath("fastqc")
        fastqc_cmd = build_fastqc_cmd(
            fastq_r1=fastq_r1,
            fastq_r2=fastq_r2,
            out_dir=fastqc_dir,
            threads=qc_cfg.get("fastqc_threads", 4),
        )
        run(fastqc_cmd)

        star_prefix = output_root.joinpath("bam", f"{sample}.")
        star_cmd = build_star_cmd(
            sample=sample,
            fastq_r1=fastq_r1,
            fastq_r2=fastq_r2,
            genome_dir=Path(ref_cfg["star_index_dir"]),
            out_prefix=star_prefix,
            threads=align_cfg.get("star_threads", 4),
        )
        run(star_cmd)

        bam_path = output_root.joinpath("bam", f"{sample}.Aligned.sortedByCoord.out.bam")
        bam_paths.append(bam_path)

    counts_dir = output_root.joinpath("counts")
    counts_dir.mkdir(parents=True, exist_ok=True)
    counts_path = counts_dir.joinpath("gene_counts.txt")

    stranded_map = {"unstranded": 0, "stranded": 1, "reverse": 2}
    first_row = manifest.iloc[0]
    strandedness = str(first_row.get("strandedness", "unstranded"))
    stranded_code = stranded_map.get(strandedness, 0)

    fc_cmd = build_featurecounts_cmd(
        bam_paths=bam_paths,
        gtf=Path(ref_cfg["gtf_annotation"]),
        out_file=counts_path,
        threads=align_cfg.get("star_threads", 4),
        stranded=stranded_code,
    )
    run(fc_cmd)

    multiqc_cmd = f"multiqc {output_root} -o {output_root}"
    run(multiqc_cmd)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run RNA QC→STAR→featureCounts workflow from config and manifest")
    parser.add_argument(
        "--config",
        required=True,
        help="Path to pipeline_config.yaml",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    rna_workflow(Path(args.config))

