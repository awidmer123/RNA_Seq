
### 5. Exploratory data analysis

#MAIN GOAL: finding the genes that lead to differences between WT and our sampels etc.........


#install necessary package
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages("pheatmap")
library(ggplot2)
library(DESeq2)
library(pheatmap)
library(ggplot2)

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




### 6. Differential expression analysis

WT_con_VS_WT_case <- results(counts_dds_DESeq, contrast = c("condition", "Lung_WT_Control", "Lung_WT_Case"))

DKO_con_VS_DKO_case <- results(counts_dds_DESeq, contrast = c("condition", "Lung_DKO_Control", "Lung_DKO_Case"))


#creating volcano plot for the differential gene expression 

# Take one results object (switch the input to make it for the other groups as well.)
res <- WT_con_VS_WT_case

# Convert to data.frame and keep gene names
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Remove rows with NA padj
res_df <- res_df[!is.na(res_df$padj), ]

# Define significance + direction
res_df$regulation <- "NS"
res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange > 0] <- "Up"
res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange < 0] <- "Down"

# Volcano plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c("Down" = "blue", "NS" = "grey70", "Up" = "red")) +
  labs(
    title = "Volcano plot: Lung_WT_Control vs Lung_WT_Case",
    x = "log2 fold change (Case vs Control)",
    y = "-log10(adjusted p-value)",
    color = "Regulation"
  ) +
  theme_minimal()



# DE genes (padj < 0.05)
de_genes <- res_df[res_df$padj < 0.05, ]

nrow(de_genes)                      # number DE-Gene
sum(de_genes$log2FoldChange > 0)    # Up-regulated
sum(de_genes$log2FoldChange < 0)    # Down-regulated

