#!/usr/bin/env bash

#SBATCH --mail-type=end,error
#SBATCH --mail-user=andri.widmer@students.unibe.ch
#SBATCH --output=/data/users/awidmer/RNA_Seq/logs/2_Mapping/hisat2/slurm-alignement-gtf-process-%j.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=rna_seq_hisat2_process

CONTAINER=/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif
# MAKE SURE TO HAVE THE UNZIPPED VERSION OF ANNOTATION (.gtf and not .gtf.gz)

GTF_FILE=/data/users/awidmer/RNA_Seq/data/raw_data/mouse_genome/Mus_musculus.GRCm39.115.gtf
OUT_DIR=/data/users/awidmer/RNA_Seq/results/exons_and_splice_sites

apptainer exec $CONTAINER hisat2_extract_splice_sites.py $GTF_FILE > $OUT_DIR/splice_sites.txt
apptainer exec $CONTAINER hisat2_extract_exons.py  $GTF_FILE > $OUT_DIR/exons.txt