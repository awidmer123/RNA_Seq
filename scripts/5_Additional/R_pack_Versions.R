# collect_R_versions.R
out_file <- "results/R_package_versions.txt"
dir.create("results", showWarnings = FALSE)

pkgs <- c(
  "BiocManager",
  "DESeq2",
  "EnhancedVolcano",
  "pheatmap",
  "ggplot2",
  "clusterProfiler",
  "org.Mm.eg.db"
)

lines <- c(
  "R and Bioconductor versions used in RNA-seq analysis",
  paste("Generated on:", format(Sys.time())),
  "--------------------------------------------------",
  "",
  "R version:",
  R.version.string,
  "",
  "Bioconductor version:",
  if (requireNamespace("BiocManager", quietly = TRUE)) as.character(BiocManager::version()) else "BiocManager not installed",
  "",
  "Package versions:"
)

for (p in pkgs) {
  if (requireNamespace(p, quietly = TRUE)) {
    lines <- c(lines, paste0(p, ": ", as.character(packageVersion(p))))
  } else {
    lines <- c(lines, paste0(p, ": NOT INSTALLED"))
  }
}

writeLines(lines, con = out_file)
cat("Wrote:", out_file, "\n")
