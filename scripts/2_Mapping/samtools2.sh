#!/usr/bin/env bash

SAM_DIR=/data/users/awidmer/RNA_Seq/results/sam2
OUT_DIR=/data/users/awidmer/RNA_Seq/results/bam2

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SLURM_SCRIPT="${SCRIPT_DIR}/samtools2.slurm"

mkdir -p "$OUT_DIR"

for sam in "$SAM_DIR"/*.sam; do
    sample=$(basename "$sam" .sam)
    out="${OUT_DIR}/${sample}.bam"

    echo "Submitting $sample"
    sbatch "$SLURM_SCRIPT" "$sam" "$out"
done
