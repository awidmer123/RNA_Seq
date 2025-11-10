#!/usr/bin/env bash

#SBATCH --mail-type=end,error
#SBATCH --mail-user=andri.widmer@students.unibe.ch
#SBATCH --output=/data/users/awidmer/RNA_Seq/logs/2_Mapping/index/slurm-alignment-index-%j.out
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=rna_seq_mapping_index

# PATH
# MAKE SURE TO HAVE THE UNZIPPED VERSION OF REFERENCE (.fa and not .fa.gz)
REF_FILE=/data/users/awidmer/RNA_Seq/data/raw_data/mouse_genome/Mus_musculus.GRCm39.dna.primary_assembly.fa
INDEX_BASENAME=/data/users/awidmer/RNA_Seq/results/index/index

# CONTAINERS PATH
CONTAINER=/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif

apptainer exec $CONTAINER hisat2-build $REF_FILE $INDEX_BASENAME -p ${SLURM_CPUS_PER_TASK}