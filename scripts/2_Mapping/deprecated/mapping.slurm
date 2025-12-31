#!/usr/bin/env bash

#SBATCH --mail-type=end,error
#SBATCH --mail-user=andri.widmer@students.unibe.ch
#SBATCH --output=/data/users/awidmer/RNA_Seq/logs/2_Mapping/mapping/slurm-alignment-mapping-%j.out
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=rna_seq_alignement_mapping

FILE_1=$1
FILE_2=$2
SAMPLE_NAME=$3
INDEX_BASENAME=$4
OUT_DIR=$5
SPLICE_SITE_FILE=$6

CONTAINER=/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif

apptainer exec $CONTAINER hisat2 -p ${SLURM_CPUS_PER_TASK} \
-x ${INDEX_BASENAME} \
-1 $FILE_1 -2 $FILE_2 \
-S ${OUT_DIR}/${SAMPLE_NAME}.sam \
--summary-file ${OUT_DIR}/sam/${SAMPLE_NAME}_summary.txt \
--rna-strandness RF \
--known-splicesite-infile $SPLICE_SITE_FILE