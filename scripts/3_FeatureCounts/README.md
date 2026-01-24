# 3_FeatureCount

This directory contains scripts for quantifying gene-level read counts from aligned RNA-seq data and for integrating the resulting summaries into the overall quality control reporting.

## Contents and running order

### 1. `featureCount.slurm`
SLURM batch script running **featureCounts** to assign aligned reads to genomic features (genes) based on a provided annotation file (GTF).  
The script generates raw count matrices that serve as the primary input for all downstream differential expression analyses.

**Input:**  
- Aligned BAM files  
- Gene annotation file (GTF)

**Output:**  
- Gene-level count tables  
- featureCounts summary files (assignment statistics)


## Notes
- Scripts in this directory should be executed after read alignment has been completed.
- The generated count tables represent the main quantitative input for DESeq2 analysis.
