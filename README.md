#SeqGate

**DNA & RNA sequencing: assess quality, then proceed.**

Author: **Shishir Biwalkar**

##What this project does

This pipeline **analyzes DNA and RNA sequencing data** to determine whether they are **good enough to proceed with downstream analysis**. It runs quality control (FastQC, MultiQC), alignment (BWA for DNA, STAR for RNA), and produces **analysis-ready outputs** (GVCFs for variant calling, gene count matrix for expression). You get a clear go/no-go view from QC figures and metrics, then the same workflow delivers the inputs needed for variant discovery and differential expression.

##What this project demonstrates

Production-style genomics workflows (aligned to typical bioinformatics roles):

- **DNA**: QC → BWA-MEM alignment → GATK variant calling → GVCFs (and optional annotation).
- **RNA**: QC → STAR alignment → featureCounts → gene-level count matrix for DESeq2/edgeR.
- **Metadata**: Sample manifests and utilities keep DNA/RNA samples aligned for multi-omics.

Configuration-driven and reproducible; all steps use explicit shell commands.

##What type of data you can test

| Assay | Data type | Format | Notes |
|-------|-----------|--------|--------|
| **DNA** | Whole-genome or targeted DNA-seq | Paired-end FASTQ (`.fastq.gz`, `.fq.gz`) | One pair (R1, R2) per sample. Reference: same organism as data (e.g. *Arabidopsis*, human). |
| **RNA** | Bulk RNA-seq (mRNA) | Paired-end FASTQ (`.fastq.gz`, `.fq.gz`) | One pair per sample. Reference: genome FASTA + GTF for the organism. Strandedness in manifest (`unstranded` / `stranded` / `reverse`) for featureCounts. |

**Supported inputs:** Illumina-style paired-end FASTQs; paths and sample IDs in TSV manifests. Single-end is not wired in. Compressed (gzip) or uncompressed FASTQ are fine.

##What data we used for testing

| Purpose | Data | Source | What we used |
|---------|------|--------|--------------|
| **DNA pipeline / QC** | ERR1760144 | ENA, *Arabidopsis halleri* (PRJEB18647) | First **1,000 read pairs** (subset) from the real run so the project runs without the full ~20GB download. Full run available via `scripts/download_data.sh`. |
| **RNA pipeline / QC** | SRR1039508, SRR1039509 | SRA/ENA, human airway smooth muscle (PRJNA229998, Himes et al.) | **Minimal paired-end FASTQs** (500 pairs per run) from `scripts/create_minimal_rna_fastq.py` for quick testing; or full runs (~1.6GB each) via `scripts/download_data.sh`. |

QC figures and pipeline runs in this repo were generated using the **subset (DNA) and minimal (RNA)** test data above. For production-style runs (variant calling, differential expression), use the full downloads and the correct reference genome and GTF.

##Original datasets

Full study and run details, citations, and reference genome notes are in [docs/DATASETS.md](docs/DATASETS.md).

**Get the data** (run from project root):

```bash
bash scripts/download_data.sh
```

This downloads the FASTQs from ENA FTP into `data/raw/fastq/`. The manifests point to these files. For full DNA data (ERR1760144, ~20GB total), the download can take a long time; a **small subset** of the DNA run is also present in `data/raw/fastq/` (ERR1760144_1/2, first 1k read pairs) so the project runs without waiting for the full download. For RNA, `scripts/create_minimal_rna_fastq.py` can create minimal paired-end FASTQs (500 pairs per run) if you have not run the full download. For the DNA workflow you need an *Arabidopsis* reference; for the RNA workflow you need a human reference (e.g. GRCh38 or hg19) and GTF—see `docs/DATASETS.md`.

**Generate figures** (after data is in `data/raw/fastq/`):

```bash
pip install -r requirements.txt   # matplotlib (and optional pandas, pyyaml)
python scripts/plot_results.py
```

Outputs appear in `results/figures/`.

##Folder layout

```text
Genomics_Workflows/
  config/
    pipeline_config.yaml      #central config for DNA+RNA
  data/
    raw/
      fastq/                  #FASTQs after running scripts/download_data.sh
    metadata/
      sample_manifest_dna.tsv
      sample_manifest_rna.tsv
    reference/
      genome.fa, genome.fa.fai, genome.dict, BWA index, STAR index, genes.gtf
  results/
    dna/                      #created by dna_pipeline.py
    rna/                      #created by rna_pipeline.py
    figures/                  #summary plots
    qc/                       #additional QC exports if needed
  scripts/
    download_data.sh           #fetch original ENA/SRA FASTQs
    dna_pipeline.py            #QC→BWA→GATK workflow
    rna_pipeline.py           #QC→STAR→featureCounts workflow
    metadata_utils.py         #manifest validation and DNA/RNA join
  docs/
    DATASETS.md               #original studies, run IDs, references
    report_template.md        #outline for analytical report
  notebooks/
    (optional analysis notebooks)
```

##Configuration

The central config file is `config/pipeline_config.yaml`:

- `input`: locations of FASTQ directories, output directory, and DNA/RNA manifests.
- `reference`: reference genome FASTA, dictionary/index files, STAR genome index, GTF annotation.
- `qc`: FastQC thread count and coverage targets.
- `alignment`: thread counts and read group template for BWA; thread count for STAR.
- `gatk`: key knobs for HaplotypeCaller (threads, emit mode, base and mapping quality filters, optional intervals).
- `annotation`: placeholders for snpEff/Ensembl VEP configuration.
- `report`: where the human-readable analytical report and figures will be written.

Edit this file once to point at your environment; the same configuration drives both DNA and RNA workflows.

##DNA workflow (QC→BWA→GATK)

The DNA orchestrator lives in `scripts/dna_pipeline.py`. It:

1. Reads `config/pipeline_config.yaml` and `data/metadata/sample_manifest_dna.tsv`.
2. Runs **FastQC** for each paired-end library and writes reports to `results/dna/fastqc/`.
3. Runs **BWA-MEM** with a proper read group, then pipes into **samtools sort** and **GATK MarkDuplicates** to produce analysis-ready, indexed BAM files under `results/dna/bam/`.
4. Runs **GATK HaplotypeCaller** in GVCF mode per sample, writing to `results/dna/gvcf/`.
5. Runs **MultiQC** over the DNA results directory to aggregate QC into a single HTML report.

Run it as:

```bash
cd "/Users/shishirbiwalkar/BioInformatics Project/Genomics_Workflows"
python scripts/dna_pipeline.py --config config/pipeline_config.yaml
```

This mirrors GATK best practices at a high level: **alignment→duplicate marking→per-sample GVCF calling**, leaving joint genotyping and downstream refinement for a separate modeling layer.

##RNA workflow (QC→STAR→counts)

The RNA orchestrator lives in `scripts/rna_pipeline.py`. It:

1. Reads `config/pipeline_config.yaml` and `data/metadata/sample_manifest_rna.tsv`.
2. Runs **FastQC** on each RNA library, writing reports to `results/rna/fastqc/`.
3. Aligns reads with **STAR** against the configured genome index, producing sorted BAMs under `results/rna/bam/`.
4. Uses **featureCounts** with the GTF annotation to build a **gene-level count matrix** at `results/rna/counts/gene_counts.txt`.
5. Runs **MultiQC** over the RNA results directory for alignment/QC aggregation.

Run it as:

```bash
cd "/Users/shishirbiwalkar/BioInformatics Project/Genomics_Workflows"
python scripts/rna_pipeline.py --config config/pipeline_config.yaml
```

The generated counts matrix is ready for downstream **R/Bioconductor** differential expression and QC workflows (DESeq2, edgeR, limma-voom, etc.).

##Metadata and multi-omics alignment

Sample manifests (`sample_manifest_dna.tsv` and `sample_manifest_rna.tsv`) encode:

- `sample_id`: logical sample identifier shared across DNA and RNA.
- `fastq_r1`, `fastq_r2`: paths to raw reads under `data/raw/fastq/`.
- Library, platform, run identifiers, and RNA `strandedness` (used to set the `-s` flag in featureCounts).

The helper module `scripts/metadata_utils.py` provides:

- Schema validation for DNA and RNA manifests.
- A join function to align DNA and RNA data on `sample_id` for downstream multi-omics modeling (e.g., variant burden vs expression changes).

##Analytical reports and visualizations

**Meaningful figures** are generated from the actual data in `data/raw/fastq/` (and, when available, from MultiQC JSON after running the pipelines). From the project root, run:

```bash
python scripts/plot_results.py
```

This writes into `results/figures/`:

- **total_reads_per_sample.png** - Sequencing depth per sample (DNA vs RNA).
- **mean_quality_per_sample.png** -Per-sample mean base quality (Phred).
- **gc_content_per_sample.png** — GC % per sample.
- **assay_summary.png** - Total reads and sample count by assay (DNA vs RNA).
- **depth_and_library_size.png** — Depth / library size (reads in millions).

The script uses real read counts from the FASTQ files listed in the manifests; if you have run FastQC/MultiQC, it will use those metrics where available. The `docs/report_template.md` file outlines the full report structure (mapping rates, coverage, variant counts, etc.) for when the full pipelines have been run.

##Why this matches the job description

- **End-to-end workflows**: QC → alignment → variant calling (DNA) and QC → alignment → counts (RNA); you assess whether data are fit to proceed, then produce analysis-ready outputs.
- **Established tools**: BWA, GATK, STAR, featureCounts, FastQC, and MultiQC.
- **Reproducibility**: Single config, explicit commands, structured manifests.
- **Data and metadata management**: Manifests and utilities keep DNA/RNA aligned for multi-omics.
- **Ready for downstream analysis**: GVCFs and gene counts feed variant discovery and differential expression.

-
