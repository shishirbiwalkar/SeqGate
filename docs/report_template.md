#Genomics_Workflows Analytical Report

##Overview

This report summarizes DNA and RNA workflows that transform raw sequencing data into interpretable outputs:
- DNA: QC→BWA-MEM alignment→GATK HaplotypeCaller GVCFs
- RNA: QC→STAR alignment→featureCounts gene-level counts

##Inputs

- DNA samples manifest: `data/metadata/sample_manifest_dna.tsv`
- RNA samples manifest: `data/metadata/sample_manifest_rna.tsv`
- Reference genome and annotations under `data/reference/`

##DNA Workflow Highlights

- Mapping rates and duplicate rates per sample
- Coverage depth distribution and percent bases above target depth
- Count of called variants per sample (from GVCFs) after basic filtering (to be plugged in downstream)

##RNA Workflow Highlights

- Alignment rates per sample
- Gene-level count matrix (`results/rna/counts/gene_counts.txt`)
- Library complexity and library size summaries

##Figures

Expected figures under `results/figures/`:
- `dna_coverage_summary.png`: mean/median coverage and percentage bases above thresholds
- `dna_mapping_rates.png`: mapping and duplicate rates per sample
- `rna_library_sizes.png`: distribution of library sizes per sample
- `rna_pca_qc.png`: PCA of normalized counts (optional downstream R/Bioconductor step)

##Notes

- Workflows are config-driven via `config/pipeline_config.yaml`.
- All heavy tools (BWA, GATK, STAR, featureCounts, FastQC, MultiQC) are invoked via explicit shell commands for transparency.

