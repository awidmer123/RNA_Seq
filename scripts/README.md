
# Scripts

This directory contains all scripts used in the RNA-Seq workflow of this project.
To keep the repository reproducible and easy to navigate, scripts are organised into
five subdirectories reflecting the chronological order of the pipeline.

Each folder corresponds to one analysis stage, from raw read QC to downstream
differential expression and enrichment analyses, plus a small set of utility scripts
for version tracking and reporting.

---

## 1_QC — FastQC (raw read quality control)

Scripts in this folder are used to assess the quality of the **raw FASTQ files** prior
to alignment.

Typical tasks include:
- Running **FastQC** on all paired-end FASTQ files
- Producing per-sample HTML reports (per-base quality, GC content, adapter signals, etc.)

Output (typical):
- `*_fastqc.html` and `*_fastqc.zip` reports per sample

Notes:
- This folder contains *only* FastQC-related scripts (no MultiQC aggregation here).

---

## 2_Mapping — HISAT2 alignment + SAMtools processing

This folder contains scripts for **mapping reads to the reference genome** and for
standard post-alignment BAM processing.

Typical tasks include:
- Aligning paired-end reads using **HISAT2**
- Converting SAM → BAM (`samtools view`)
- Sorting BAM files by genomic coordinates (`samtools sort`)
- Indexing sorted BAM files (`samtools index`)
- Inspecting alignment summary metrics (e.g., overall alignment rate, proper pairs)

Output (typical):
- Sorted and indexed BAM files (`*.sorted.bam`, `*.sorted.bam.bai`)
- Alignment logs and summary statistics

Notes:
- These scripts are designed for HPC execution (Slurm) and are typically run per sample
  (often via job arrays).

---

## 3_FeatureCounts — gene-level quantification

Scripts in this folder perform **gene-level read counting** based on the aligned BAM files.

Typical tasks include:
- Running **featureCounts** with an Ensembl GTF annotation
- Applying appropriate settings for paired-end and strand-specific libraries
- Producing a count matrix with one column per sample

Output (typical):
- `counts.txt` (gene × sample count table)
- featureCounts summary output (assigned/unassigned read statistics)

Notes:
- The resulting count table is the main input for downstream R analysis (DESeq2).

---

## 4_Analysis — DESeq2, enrichGO, and “important genes”

This folder contains R scripts for **downstream statistical analysis and biological interpretation**.

Typical tasks include:
- Importing the featureCounts count table and sample metadata
- Differential expression analysis using **DESeq2**
  - contrasts such as WT case vs WT control, DKO case vs DKO control, and WT case vs DKO case
- Exploratory analysis and visualisation
  - PCA, heatmaps, volcano plots, summary tables
- Functional enrichment analysis using **clusterProfiler / enrichGO**
  - GO Biological Process enrichment on significant DEGs
  - optional redundancy reduction (e.g., `simplify()`)
- Identification and tracking of **important immune / interferon-related genes**
  - extraction of selected candidates and reporting in tables (e.g., “important genes” table)

Output (typical):
- DESeq2 result tables (raw + filtered)
- Figures (PCA/heatmap/volcano/GO barplots)
- Candidate gene tables and summaries used in the report

Notes:
- This folder contains the scripts that directly generate the key results reported in the final PDF.

---

## 5_Additional — utilities, MultiQC, and version tracking

This folder contains **supporting scripts** that improve reproducibility and reporting,
but do not belong to a single core workflow stage.

Typical contents include:
- **R package / session version tracking**
  - scripts that record package versions used in R analyses (for reproducibility)
- **MultiQC execution (Slurm)**
  - aggregation of QC metrics across samples into a single report (typically run on HPC)
- **General version/environment control scripts (Slurm)**
  - scripts that document software versions (tools, containers, environment info)

Output (typical):
- MultiQC report files
- Version logs (software + R packages)

Notes:
- This folder is intentionally separated to avoid mixing “core pipeline” scripts with utilities.

---

## General notes

- Computationally intensive steps were executed on the IBU cluster (Slurm).
- Where applicable, tools were run in containerised environments (Apptainer) to ensure stable versions.
- The directory structure is designed to mirror the workflow and make the analysis easy to reproduce.
