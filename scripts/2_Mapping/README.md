## Mapping — HISAT2 & Samtools Pipeline

This directory contains the scripts used to run the RNA-seq mapping workflow. The pipeline extracts splice/exon information, builds the HISAT2 genome index, aligns paired-end FASTQ reads, and converts the resulting SAM files into sorted/indexed BAM files using samtools inside an Apptainer container.

## Scripts overview

gtf_process.slurm
Extract splice sites and exons from the unzipped GTF file (HISAT2 helpers).

index.slurm
Build the HISAT2 genome index from the reference .fa file.

mapping.sh
Loop over paired FASTQ files and submit one mapping job per sample.

mapping.slurm
Run HISAT2 alignment using the genome index and splice site file.
Outputs <sample>.sam and a summary file.

samtools.sh
Loop over all .sam files and submit one samtools conversion job per file.

samtools.slurm
Convert .sam → .bam, sort, and create .bai index using samtools.


## Running order of involved scripts

1. sbatch index.slurm
2. sbatch hisat2_process.slurm
3. ./mapping.sh
4. ./samtools.sh

