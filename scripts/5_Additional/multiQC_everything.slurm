#!/usr/bin/env bash

#SBATCH --job-name=sample_1_lung
#SBATCH --mail-user=andri.widmer@students.unibe.ch
#SBATCH --output=/data/users/awidmer/RNA_Seq/logs/3_FeatureCounts/fastqc_featurecounts_%j.out
#SBATCH --error=/data/users/awidmer/RNA_Seq/logs/3_FeatureCounts/fastqc_featurecounts_%j.err
#SBATCH --time=02:10:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --partition=pibu_el8

#Set paths
CONTAINER="/containers/apptainer/multiqc-1.19.sif"
INPUT_DIR="/data/users/awidmer/RNA_Seq/results"
OUTPUT_DIR="/data/users/awidmer/RNA_Seq/results"

#Run FastQC using Apptainer Container
apptainer exec ${CONTAINER} multiqc \
    ${INPUT_DIR}\
    -o ${OUTPUT_DIR} \

echo "Task complete!"