# 5_Additional

This directory contains auxiliary scripts used for quality control aggregation, version tracking, and reproducibility documentation. The files in this folder do not represent core analysis steps but support validation and transparency of the overall RNA-seq workflow.

## Contents

### `multiqc_everything.slurm`
SLURM batch script to run **MultiQC** across the entire project directory structure.  
It aggregates quality control reports from multiple pipeline stages (e.g. FastQC, alignment, featureCounts) into a single summary report, providing a global overview of data quality and processing performance.

**Output:**  
- Consolidated MultiQC report (HTML and data files)

---

### `Version_control.slurm`
SLURM script used to record software versions of key command-line tools employed in the RNA-seq workflow (e.g. aligners, QC tools).  
This ensures full reproducibility by documenting the exact tool versions used during analysis.

**Output:**  
- Versions list

---

### `R_pack_Versions.R`
R script that captures the versions of all relevant R packages used in the downstream analysis (e.g. DESeq2, ggplot2, clusterProfiler).  
This script complements the system-level version tracking and is intended to be run once after the R environment is finalized.

**Output:**  
- R package version list

---

## Notes
- Scripts in this directory are independent of the main analysis flow.
- They can be executed at any time once the corresponding tools and R environment are available.
- The directory primarily serves reproducibility, documentation, and quality control purposes.
