# multiQC

This directory contains aggregated quality control results generated with
MultiQC. The reports summarize quality metrics collected across multiple
pipeline stages, providing a global overview of data quality and processing
performance.

## Contents

### `multiqc_report_1.html`
Interactive HTML report summarizing quality control metrics from tools such as
FastQC, alignment statistics, and featureCounts.  
This report is intended for visual inspection and documentation of overall data
quality.

### `multiqc_data_1/`
Directory containing the underlying data files and logs used by MultiQC to
generate the HTML report.

---

## Notes
- The MultiQC report consolidates QC results from several pipeline steps.
- No downstream analyses directly depend on these files; they serve quality
  assessment and documentation purposes only.
