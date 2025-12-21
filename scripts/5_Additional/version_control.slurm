#!/usr/bin/env bash
#SBATCH --job-name=tool_versions
#SBATCH --output=/data/users/awidmer/RNA_Seq/logs/5_Versions/tool_versions_%j.out
#SBATCH --error=/data/users/awidmer/RNA_Seq/logs/5_Versions/tool_versions_%j.err
#SBATCH --time=00:05:00
#SBATCH --mem=256M
#SBATCH --cpus-per-task=1
#SBATCH --partition=pibu_el8

# Output file
OUTFILE="/data/users/awidmer/RNA_Seq/results/versions/tool_versions.txt"

# Truncate output file
echo "Tool versions used in RNA-Seq analysis" > "$OUTFILE"
echo "Generated on: $(date)" >> "$OUTFILE"
echo "------------------------------------" >> "$OUTFILE"
echo >> "$OUTFILE"

# Helper function
run_version () {
    local label="$1"
    local container="$2"
    local cmd="$3"

    echo "### ${label}" >> "$OUTFILE"
    echo "Container: ${container}" >> "$OUTFILE"

    if apptainer exec "$container" bash -c "$cmd" >> "$OUTFILE" 2>&1; then
        echo >> "$OUTFILE"
    else
        echo "ERROR: Could not retrieve version for ${label}" >> "$OUTFILE"
        echo >> "$OUTFILE"
    fi
}

# Apptainer version
echo "### Apptainer" >> "$OUTFILE"
apptainer --version >> "$OUTFILE" 2>&1 || echo "Apptainer version not available" >> "$OUTFILE"
echo >> "$OUTFILE"

# FastQC
run_version \
  "FastQC" \
  "/containers/apptainer/fastqc-0.12.1.sif" \
  "fastqc --version"

# HISAT2
run_version \
  "HISAT2" \
  "/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif" \
  "hisat2 --version"

# Samtools
run_version \
  "Samtools" \
  "/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif" \
  "samtools --version | head -n 1"

# featureCounts (Subread)
run_version \
  "featureCounts (Subread)" \
  "/containers/apptainer/subread_2.0.6.sif" \
  "featureCounts -v"

# MultiQC
run_version \
  "MultiQC" \
  "/containers/apptainer/multiqc-1.19.sif" \
  "multiqc --version"

echo "Done. Versions written to ${OUTFILE}"
