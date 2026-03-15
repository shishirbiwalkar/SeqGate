#Original datasets

Author: **Shishir Biwalkar**

This project uses **public sequencing data** from ENA/SRA. All paths and run IDs below are used by the manifests and the download script.

##DNA workflow

| Run ID      | Study / project        | Organism           | Layout   | Use in pipeline |
|-------------|------------------------|--------------------|----------|------------------|
| **ERR1760144** | PRJEB18647 (Arabidopsis halleri) | *Arabidopsis halleri* | Paired-end | Single sample in `sample_manifest_dna.tsv` |

- **Source**: [ENA ERR1760144](https://www.ebi.ac.uk/ena/browser/view/ERR1760144), [Bioinformatics Workbook GATK DNA-Seq tutorial](https://bioinformaticsworkbook.org/dataAnalysis/VariantCalling/gatk-dnaseq-best-practices-workflow.html).
- **Reference**: Use an *Arabidopsis thaliana* (or *A. halleri*) reference genome and build BWA index and GATK dictionary/fai. The pipeline expects `config/reference` paths to point to that reference.

##RNA workflow

| Run ID        | Study / project        | Organism    | Layout   | Use in pipeline |
|---------------|------------------------|-------------|----------|------------------|
| **SRR1039508** | PRJNA229998 (Himes et al. airway smooth muscle) | *Homo sapiens* | Paired-end | Sample 1 in `sample_manifest_rna.tsv` |
| **SRR1039509** | PRJNA229998 (same study) | *Homo sapiens* | Paired-end | Sample 2 in `sample_manifest_rna.tsv` |

- **Source**: [SRA SRR1039508](https://www.ncbi.nlm.nih.gov/sra/SRR1039508), [SRR1039509](https://www.ncbi.nlm.nih.gov/sra/SRR1039509); study [PRJNA229998](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA229998), GEO GSE52778.
- **Citation**: Himes BE, Jiang X, Wagner P et al. RNA-Seq transcriptome profiling identifies CRISPLD2 as a glucocorticoid responsive gene that modulates cytokine function in airway smooth muscle cells. *PLoS One* 2014;9(6):e99625.
- **Reference**: Human genome (e.g. GRCh38 or hg19). The Himes study used hg19; use the same reference and matching GTF for STAR and featureCounts.

##Downloading the data

From the project root, run:

```bash
bash scripts/download_data.sh
```

This fetches the FASTQ files from ENA FTP into `data/raw/fastq/` and names them by run ID (e.g. `ERR1760144_1.fastq.gz`, `SRR1039508_1.fastq.gz`). The manifests point to these paths.

##File paths after download

- DNA: `data/raw/fastq/ERR1760144_1.fastq.gz`, `data/raw/fastq/ERR1760144_2.fastq.gz`
- RNA: `data/raw/fastq/SRR1039508_1.fastq.gz`, `data/raw/fastq/SRR1039508_2.fastq.gz`, `data/raw/fastq/SRR1039509_1.fastq.gz`, `data/raw/fastq/SRR1039509_2.fastq.gz`
