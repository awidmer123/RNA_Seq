#!/usr/bin/env bash

# PATH
SAM_DIR=/data/users/awidmer/RNA_Seq/results/sam
OUT_DIR=/data/users/awidmer/RNA_Seq/results/bam
REF_FILE=/data/users/awidmer/RNA_Seq/data/raw_data/mouse_genome/Mus_musculus.GRCm39.dna.primary_assembly.fa

# MAP EACH FASTQ TO REF
for file in `ls -1 $SAM_DIR/*.sam`; do
    sbatch ${BASH_SOURCE[0]%/samtools.sh}/samtools.slurm $file $REF_FILE ${OUT_DIR}/$(basename ${file%.sam}).bam
done