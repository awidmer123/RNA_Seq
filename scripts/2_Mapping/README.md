# 2_Mapping

This directory contains scripts for read alignment and post-processing in the
RNA-seq workflow. The final pipeline performs splice site extraction, genome
indexing, paired-end read alignment using HISAT2, and conversion of SAM files
into sorted and indexed BAM files using samtools within an Apptainer container.

Earlier development iterations have been moved to a dedicated `deprecated/`
subdirectory and are not part of the final analysis.

## Contents and running order

### 1. `hisat2_process.slurm`
SLURM batch script extracting splice sites and exon information from the
uncompressed GTF annotation file using HISAT2 helper scripts.

**Output:**  
- Splice site and exon files for downstream alignment

---

### 2. `index.slurm`
SLURM script building the HISAT2 genome index from the reference genome FASTA
file.

**Output:**  
- HISAT2 genome index files

---

### 3. `mapping2_all`
Wrapper script iterating over paired-end FASTQ files and submitting one mapping
job per sample.

---

### 4. `mapping2.slurm`
SLURM batch script performing paired-end read alignment using HISAT2 with the
pre-built genome index and splice site information.

**Input:**  
- Paired-end FASTQ files  
- HISAT2 genome index  
- Splice site/exon files

**Output:**  
- SAM alignment files  
- Alignment summary statistics

---

### 5. `samtools2.sh`
Wrapper script submitting samtools jobs for all generated SAM files.

---

### 6. `samtools2.slurm`
SLURM batch script converting SAM files to BAM format, sorting them, and
generating BAM index files.

**Output:**  
- Sorted BAM files  
- BAM index (`.bai`) files

---

## Notes
- Scripts should be executed in the order listed above.
- This directory represents the final, reproducible mapping pipeline used in
  the analysis.
- Output BAM files are used as input for read quantification in the
  `3_FeatureCounts` directory.
