#Reference genomes and indexes

This directory must contain the reference files used by the pipelines. The **original datasets** use different organisms for DNA and RNA (see [docs/DATASETS.md](../docs/DATASETS.md)):

- **DNA workflow** (manifest: ERR1760144): Use an *Arabidopsis thaliana* (or *A. halleri*) reference: FASTA, `samtools faidx`, `gatk CreateSequenceDictionary`, and BWA index. Point `config/pipeline_config.yaml` `reference.genome_fa` and `reference.bwa_index_prefix` here.
- **RNA workflow** (manifests: SRR1039508, SRR1039509): Use a human reference (e.g. GRCh38 or hg19) and matching GTF. Build a STAR index and point `reference.star_index_dir` and `reference.gtf_annotation` here.

You can use one reference layout for DNA and another for RNA; just ensure the config paths match the organism for each workflow.
