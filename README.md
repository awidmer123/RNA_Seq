# RNA_Seq

This repository contains the complete RNA-Seq analysis workflow for the project  
**“The Role of Interferon Signaling in the Pulmonary Immune Response”**, conducted as part of the RNA sequencing practical course (2025).

The analysis focuses on lung tissue from wild-type (WT) and interferon receptor double knockout (DKO) mice under infected (“case”) and uninfected (“control”) conditions, and follows a reproducible, stepwise RNA-Seq workflow from raw data quality control to downstream statistical and functional analyses.

---
## Repository structure
---

## scripts/

Contains all scripts used throughout the RNA-Seq workflow.  
Scripts are organised into five subdirectories reflecting the chronological order of the pipeline:

- **1_QC** – FastQC (raw read quality control)  
- **2_Mapping** – Read alignment (HISAT2) and BAM processing (SAMtools)  
- **3_FeatureCounts** – Gene-level read quantification  
- **4_Analysis** – DESeq2, enrichGO, and analysis of important immune-related genes  
- **5_Additional** – Utility scripts (MultiQC, R package versions, version control)

A detailed description of the contents and purpose of each subdirectory is provided in  
`scripts/README.md`.

---

## report/

Contains the final analysis report:
- The **Quarto source file (`.qmd`)**, including code, text, and figure generation
- The **rendered PDF report**, produced from the Quarto document

This folder represents the final, human-readable output of the project.

---

## resources/

Contains supporting documentation for the analysis, including:
- A **workflow overview** describing the RNA-Seq analysis pipeline and processing steps

This folder serves as a reference for understanding the overall structure and logic of the workflow.

---

## results/

Contains output files and metadata generated during the analysis, including:
- **MultiQC reports** summarising quality control metrics across samples
- **Version information**, documenting software and package versions used in the workflow

This folder supports transparency and reproducibility but does not contain raw data.

---

## Reproducibility and execution

- Computationally intensive steps were executed on the IBU high-performance computing cluster using Slurm.
- Where applicable, software tools were run in containerised environments (Apptainer) to ensure stable and reproducible versions.
- Downstream statistical analyses were performed locally in RStudio.
- Software and package versions are recorded in the `results/versions/` directory.

---

## Data origin

The RNA-Seq data analysed in this project originate from the study by Singhania et al. (2019).  
Subsamples were provided by the university for re-analysis within the practical course framework.

---

## Author

Andri L. Widmer  
University of Fribourg
