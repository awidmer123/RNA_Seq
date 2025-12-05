#!/usr/bin/env bash

#SBATCH --mail-type=end,error
#SBATCH --mail-user=andri.widmer@students.unibe.ch
#SBATCH --output=/data/users/awidmer/RNA_Seq/logs/2_Mapping/sam_out/slurm-samtools2-%j.out
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=rna_seq_samtools_all

SAM_FILE="$1"
OUT_FILE="$2"
UNSORTED_FILE="${OUT_FILE%.bam}.unsorted.bam"

CONTAINER=/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif

echo "Processing $SAM_FILE → $OUT_FILE"

# 1) SAM -> BAM conversion
apptainer exec --bind /data/ "$CONTAINER" \
    samtools view -hbS "$SAM_FILE" > "$UNSORTED_FILE"

# 2) Sorting
apptainer exec --bind /data/ "$CONTAINER" \
    samtools sort \
        -@ "$SLURM_CPUS_PER_TASK" \
        -m 3G \
        -o "$OUT_FILE" \
        "$UNSORTED_FILE"

rm -f "$UNSORTED_FILE"

# 3) BAM index
apptainer exec --bind /data/ "$CONTAINER" \
    samtools index -@ "$SLURM_CPUS_PER_TASK" "$OUT_FILE"

echo "Finished: $OUT_FILE"
