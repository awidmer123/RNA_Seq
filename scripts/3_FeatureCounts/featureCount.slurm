#!/usr/bin/env bash

#SBATCH --mail-type=end,error
#SBATCH --mail-user=andri.widmer@students.unibe.ch
#SBATCH --output=/data/users/awidmer/RNA_Seq/logs/3_FeatureCounts/slurm-featureCount-%j.out
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=featureCounts

CONTAINER=/containers/apptainer/subread_2.0.6.sif
ANNOTATION=data/raw_data/mouse_genome/Mus_musculus.GRCm39.115.gtf
OUT_FILE=results/counts/counts.txt

apptainer exec $CONTAINER featureCounts \
  -p --countReadPairs \
  -B -C \
  -s 2 \
  -t exon \
  -g gene_id \
  -a ${ANNOTATION} \
  -o ${OUT_FILE} \
  results/bam/*.bam


