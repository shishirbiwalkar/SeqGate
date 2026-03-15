#!/usr/bin/env python
#Create minimal valid paired-end FASTQ for RNA demo when full download is not run.
#Writes data/raw/fastq/SRR1039508_1.fastq.gz and _2.fastq.gz (and SRR1039509) with 500 read pairs.
#Author: Shishir Biwalkar (project)

import gzip
from pathlib import Path

BASE = "GGATTTGAAGACGGCATACGAGAT"
QUAL = "I" * 75  #Phred 40

def write_paired_fastq(run: str, fastq_dir: Path, n_reads: int = 500) -> None:
    fastq_dir.mkdir(parents=True, exist_ok=True)
    r1_path = fastq_dir / f"{run}_1.fastq.gz"
    r2_path = fastq_dir / f"{run}_2.fastq.gz"
    seq1 = (BASE + "A" * 75)[:75]
    seq2 = ("C" * 20 + "GTAC" * 14)[:75]
    with gzip.open(r1_path, "wt") as f1, gzip.open(r2_path, "wt") as f2:
        for i in range(n_reads):
            f1.write(f"@{run}.{i}/1\n{seq1}\n+\n{QUAL}\n")
            f2.write(f"@{run}.{i}/2\n{seq2}\n+\n{QUAL}\n")
    print(f"Wrote {r1_path} and {r2_path} ({n_reads} pairs)")

def main() -> None:
    root = Path(__file__).resolve().parent.parent
    fastq_dir = root / "data" / "raw" / "fastq"
    for run in ["SRR1039508", "SRR1039509"]:
        write_paired_fastq(run, fastq_dir)

if __name__ == "__main__":
    main()
