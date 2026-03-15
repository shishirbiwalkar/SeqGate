#!/usr/bin/env bash
#Download original ENA/SRA FASTQ files into data/raw/fastq.
#Run from project root: bash scripts/download_data.sh
#Author: Shishir Biwalkar (project)

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
FASTQ_DIR="$PROJECT_ROOT/data/raw/fastq"
ENA_BASE="ftp://ftp.sra.ebi.ac.uk/vol1/fastq"

mkdir -p "$FASTQ_DIR"
cd "$FASTQ_DIR"

#DNA: ERR1760144 (Arabidopsis halleri, paired-end; GATK tutorial dataset)
echo "Downloading DNA run ERR1760144..."
curl -# -L -o ERR1760144_1.fastq.gz "$ENA_BASE/ERR176/004/ERR1760144/ERR1760144_1.fastq.gz"
curl -# -L -o ERR1760144_2.fastq.gz "$ENA_BASE/ERR176/004/ERR1760144/ERR1760144_2.fastq.gz"

#RNA: SRR1039508, SRR1039509 (Himes et al. human airway smooth muscle, paired-end)
echo "Downloading RNA run SRR1039508..."
curl -# -L -o SRR1039508_1.fastq.gz "$ENA_BASE/SRR103/008/SRR1039508/SRR1039508_1.fastq.gz"
curl -# -L -o SRR1039508_2.fastq.gz "$ENA_BASE/SRR103/008/SRR1039508/SRR1039508_2.fastq.gz"
echo "Downloading RNA run SRR1039509..."
curl -# -L -o SRR1039509_1.fastq.gz "$ENA_BASE/SRR103/009/SRR1039509/SRR1039509_1.fastq.gz"
curl -# -L -o SRR1039509_2.fastq.gz "$ENA_BASE/SRR103/009/SRR1039509/SRR1039509_2.fastq.gz"

echo "Done. FASTQs are in $FASTQ_DIR"
ls -la "$FASTQ_DIR"
