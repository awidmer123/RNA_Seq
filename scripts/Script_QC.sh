#!/usr/bin/env bash
#SBATCH --job-name=sample_1_lung
#SBATCH --output=/data/users/awidmer/RNA_Seq/results/QC/fastqc_%j.out
#SBATCH --error=/data/users/awidmer/RNA_Seq/logs/QC/fastqc_%j.err
#SBATCH --time=02:10:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --partition=pibu_el8

#Set paths
CONTAINER="/containers/apptainer/fastqc-0.12.1.sif"
INPUT_DIR="/data/users/awidmer/RNA_Seq/data/raw_data/reads_Lung"
OUTPUT_DIR="/data/users/awidmer/RNA_Seq/results/QC"

#Run FastQC using Apptainer Container
apptainer exec ${CONTAINER} fastqc \
    ${INPUT_DIR}/*.fastq.gz \
    -o ${OUTPUT_DIR} \
    -t $SLURM_CPUS_PER_TASK

echo "Task complete!"