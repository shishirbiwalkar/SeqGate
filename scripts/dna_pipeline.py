#!/usr/bin/env python
"""
DNA workflow orchestrator for Genomics_Workflows.

Transforms raw FASTQ into analysis-ready BAM and GVCF using:
  - FastQC/MultiQC for QC
  - BWA-MEM for alignment
  - samtools + GATK for post-processing and variant calling

This script does not execute heavy tools inside Python; instead it assembles
explicit shell commands so runs are transparent and reproducible.
"""

import argparse
import subprocess
import textwrap
from pathlib import Path

import yaml
import pandas as pd


def run(cmd: str) -> None:
    #thin wrapper to print and run shell commands
    print(f"\n[cmd] {cmd}")
    subprocess.run(cmd, shell=True, check=True)


def load_config(config_path: Path) -> dict:
    with config_path.open() as fh:
        return yaml.safe_load(fh)


def load_manifest(manifest_path: Path) -> pd.DataFrame:
    df = pd.read_csv(manifest_path, sep="\t")
    required = {"sample_id", "fastq_r1", "fastq_r2"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Manifest missing required columns: {sorted(missing)}")
    return df


def prepare_output_dirs(output_root: Path) -> None:
    for sub in ["fastqc", "bam", "gvcf", "logs", "tmp"]:
        (output_root / sub).mkdir(parents=True, exist_ok=True)


def build_fastqc_cmd(fastq_r1: Path, fastq_r2: Path, out_dir: Path, threads: int) -> str:
    return (
        f"fastqc --threads {threads} "
        f"--outdir {out_dir} "
        f"{fastq_r1} {fastq_r2}"
    )


def build_bwa_cmd(sample: str, fastq_r1: Path, fastq_r2: Path, ref_prefix: Path, out_bam: Path, threads: int, rg_template: str) -> str:
    rg = rg_template.format(sample=sample)
    bam_unsorted = out_bam.with_suffix(".unsorted.bam")
    cmd = textwrap.dedent(
        f"""
        bwa mem -t {threads} -R '{rg}' {ref_prefix} {fastq_r1} {fastq_r2} \
        | samtools view -bS - \
        | samtools sort -@ {threads} -o {bam_unsorted} -
        && gatk MarkDuplicates \
           -I {bam_unsorted} \
           -O {out_bam} \
           -M {out_bam.with_suffix('.metrics.txt')} \
           --CREATE_INDEX true \
        && rm -f {bam_unsorted}
        """
    ).strip()
    return cmd


def build_haplotypecaller_cmd(sample: str, ref_fa: Path, bam: Path, out_gvcf: Path, cfg: dict) -> str:
    gatk_cfg = cfg["gatk"]
    emit_mode = gatk_cfg.get("emit_mode", "GVCF")
    mbq = gatk_cfg.get("min_base_quality", 20)
    mmq = gatk_cfg.get("min_mapping_quality", 20)
    threads = gatk_cfg.get("haplotypecaller_threads", 4)
    intervals = gatk_cfg.get("intervals")

    interval_part = f"-L {intervals}" if intervals else ""

    cmd = textwrap.dedent(
        f"""
        gatk HaplotypeCaller \
          -R {ref_fa} \
          -I {bam} \
          -O {out_gvcf} \
          -ERC {emit_mode} \
          --native-pair-hmm-threads {threads} \
          --min-base-quality-score {mbq} \
          --minimum-mapping-quality {mmq} \
          {interval_part}
        """
    ).strip()
    return cmd


def build_multiqc_cmd(output_root: Path) -> str:
    return f"multiqc {output_root} -o {output_root}"


def dna_workflow(config_path: Path) -> None:
    cfg = load_config(config_path)
    input_cfg = cfg["input"]
    ref_cfg = cfg["reference"]
    qc_cfg = cfg["qc"]
    align_cfg = cfg["alignment"]

    manifest_path = Path(input_cfg["dna_manifest"])
    output_root = Path(input_cfg["output_dir"]).joinpath("dna")

    prepare_output_dirs(output_root)
    manifest = load_manifest(manifest_path)

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

        bam_dir = output_root.joinpath("bam")
        bam_path = bam_dir.joinpath(f"{sample}.bam")
        bam_dir.mkdir(parents=True, exist_ok=True)
        bwa_cmd = build_bwa_cmd(
            sample=sample,
            fastq_r1=fastq_r1,
            fastq_r2=fastq_r2,
            ref_prefix=Path(ref_cfg["bwa_index_prefix"]),
            out_bam=bam_path,
            threads=align_cfg.get("bwa_threads", 4),
            rg_template=align_cfg["read_group_template"],
        )
        run(bwa_cmd)

        gvcf_dir = output_root.joinpath("gvcf")
        gvcf_dir.mkdir(parents=True, exist_ok=True)
        gvcf_path = gvcf_dir.joinpath(f"{sample}.g.vcf.gz")
        hc_cmd = build_haplotypecaller_cmd(
            sample=sample,
            ref_fa=Path(ref_cfg["genome_fa"]),
            bam=bam_path,
            out_gvcf=gvcf_path,
            cfg=cfg,
        )
        run(hc_cmd)

    multiqc_cmd = build_multiqc_cmd(output_root)
    run(multiqc_cmd)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run DNA QC→BWA→GATK workflow from config and manifest")
    parser.add_argument(
        "--config",
        required=True,
        help="Path to pipeline_config.yaml",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    dna_workflow(Path(args.config))

