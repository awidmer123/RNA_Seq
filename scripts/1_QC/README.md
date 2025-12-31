# 1_QC

This directory contains scripts used for initial quality control of raw RNA-seq
data prior to read alignment. The quality control step is intended to identify
potential issues such as low base quality, adapter contamination, or abnormal
sequence composition before downstream processing.

## Contents and running order

### 1. `script_1_QC.sh`
Shell script performing quality control on raw FASTQ files using **FastQC** and
aggregating the results with **MultiQC**.  
The script is designed to be executed once all raw sequencing files are
available.

**Input:**  
- Raw FASTQ files

**Output:**  
- Individual FastQC reports  
- Aggregated MultiQC quality control report

---

## Notes
- This step should be completed before read alignment and quantification.
- The MultiQC report provides a global overview of sequencing quality across all
  samples.
- No files in this directory are modified after QC; all downstream steps rely
  on the original raw data.
