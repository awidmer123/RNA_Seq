#!/usr/bin/env bash

#SBATCH --mail-type=end,error
#SBATCH --mail-user=andri.widmer@students.unibe.ch
#SBATCH --output=/data/users/awidmer/RNA_Seq/logs/2_Mapping/sam_out/slurm-samtools-%j.out
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=rna_seq_alignment_samtools

SAM_FILE=$1
REF_FILE=$2
OUT_FILE=$3
UNSORTED_FILE=${OUT_FILE%.bam}.unsorted.bam
CONTAINER=/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif

#index reference for samtools (make the fai)
apptainer exec $CONTAINER samtools faidx $REF_FILE

#convert sam to bam
apptainer exec $CONTAINER samtools view -hb -t ${REF_FILE}.fai ${SAM_FILE} > ${UNSORTED_FILE}

#sort
apptainer exec $CONTAINER samtools sort -@ ${SLURM_CPUS_PER_TASK} -o $OUT_FILE -m 8G ${UNSORTED_FILE}
rm ${UNSORTED_FILE}

#index bam (make the bai)
apptainer exec $CONTAINER samtools index -@ ${SLURM_CPUS_PER_TASK} $OUT_FILE