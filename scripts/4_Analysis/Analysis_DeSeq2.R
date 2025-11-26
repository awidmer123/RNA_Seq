
#MAIN GOAL: finding the genes that lead to differences between WT and our sampels etc.........


#install necessary package
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages("pheatmap")
library(ggplot2)
library(DESeq2)
library(pheatmap)

#reading in the necessary files
counts <- read.table("Desktop/UNIFR/2_Bioinf/RNA_Sequencing/counts.txt", header = TRUE, row.names = 1)
groups_info <- read.table("Desktop/UNIFR/2_Bioinf/RNA_Sequencing/groups_info.txt")


#Adjusting the table and formats of the input files
counts_cut <- counts[ , 6:ncol(counts) ]
groups_info <- groups_info[-7,] #just because in the first run i didn't get the Srr file from group ***24.
names(groups_info) <- c("label", "condition")

#creating the DESeq2 object, run the command DESeq and runn rlog()
counts_dds <- DESeqDataSetFromMatrix(countData = counts_cut, colData = groups_info, design = ~condition)
counts_dds_DESeq <- DESeq(counts_dds)
counts_dds_rlog <- rlog(counts_dds_DESeq)


#Finding suitable ways to show data
# 1. pca
# 2. heatmap

# 1. performing the pca
#log transformation of the data for the pca
plotPCA(object = counts_dds_rlog, intgroup = "condition")

# 2. heat map
mat <- assay(counts_dds_rlog)
head(rownames(mat))
# calculate variance per gene
gene_vars <- apply(mat, 1, var)

# keep the top 500 most variable genes
top_genes <- order(gene_vars, decreasing = TRUE)[1:20]
mat_top <- mat[top_genes, ]

#building a proper annotation data frame
annotation_df <- as.data.frame(colData(counts_dds_DESeq)[, "condition", drop = FALSE])

# Very important: rownames must match the samples in mat_top
rownames(annotation_df) <- colnames(mat_top)



pheatmap(
  mat_top,
  annotation_col = annotation_df,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = TRUE,
  fontsize_row = 6
)







