#!/usr/bin/env bash

# Verzeichnisse
FASTQ_DIR=/data/users/awidmer/RNA_Seq/data/raw_data/reads_Lung
OUT_DIR=/data/users/awidmer/RNA_Seq/results/sam2
INDEX_BASENAME=/data/users/awidmer/RNA_Seq/results/index/index
SPLICE_SITE_FILE=/data/users/awidmer/RNA_Seq/results/exons_and_splice_sites/splice_sites.txt

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SLURM_SCRIPT="${SCRIPT_DIR}/mapping2.slurm"

mkdir -p "$OUT_DIR"

for file_1 in "$FASTQ_DIR"/*_1.fastq.gz; do
    file_2=${file_1%_1.fastq.gz}_2.fastq.gz
    sample_name=$(basename "${file_1%_1.fastq.gz}")

    echo "Submitting mapping job for: $sample_name"
    sbatch "$SLURM_SCRIPT" \
        "$file_1" \
        "$file_2" \
        "$sample_name" \
        "$INDEX_BASENAME" \
        "$OUT_DIR" \
        "$SPLICE_SITE_FILE"
done
