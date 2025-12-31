
### Exploratory data analysis

#MAIN GOAL: finding the genes that lead to differences between WT and our sampels etc.........


#loading necessary packages for general data anlysis

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
counts <- read.table("counts2.txt", header = TRUE, row.names = 1) 
groups_info <- read.table("groups_info.txt")


#Adjusting the table and formats of the input files
counts_cut <- counts[ , 6:ncol(counts) ]
names(groups_info) <- c("label", "condition")

#creating the DESeq2 object, run the command DESeq and runn rlog()
counts_dds <- DESeqDataSetFromMatrix(countData = counts_cut, colData = groups_info, design = ~condition)
counts_dds_DESeq <- DESeq(counts_dds)
counts_dds_rlog <- rlog(counts_dds_DESeq)

#---------------------------------------------------------------

#Finding suitable ways to show data

# 1. pca
# 2. heatmap
# 3. volcano plots

#---------------------------------------------------------------

# 1. performing the pca
#log transformation of the data for the pca

pca_plot <- plotPCA(counts_dds_rlog, intgroup = "condition") +
  theme_bw() +
  labs(
    title = "PCA of rlog-transformed counts",
    x = "PC1",
    y = "PC2"
  )
pca_plot
#we have to logtransform the numbers, to make the variance independent of the sample mean. ;)
#so since we got genes that are really highly expressed and some low, that the "absolute change is reduced to relative change"


#sometimes the adjusted p-value will be NA even though theres a number for the normal p value. thats because the "sample size" sometimes are too low.......

#---------------------------------------------------------------

# 2. heat map
mat <- assay(counts_dds_rlog)

# calculate variance per gene
gene_vars <- apply(mat, 1, var)

# keep the top 50 most variable genes
top_genes <- order(gene_vars, decreasing = TRUE)[1:200]
mat_top <- mat[top_genes, ]

#building a proper annotation data frame
annotation_df <- as.data.frame(colData(counts_dds_DESeq)[, "condition", drop = FALSE])

# Very important: rownames must match the samples in mat_top
rownames(annotation_df) <- colnames(mat_top)



heatmap <- pheatmap(
  mat_top,
  annotation_col = annotation_df,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = FALSE,
  fontsize_row = 3
)
#
heatmap

#---------------------------------------------------------------
### 6. Differential expression analysis


# WT control vs WT case
WT_con_VS_WT_case <- results(counts_dds_DESeq, contrast = c("condition", "Lung_WT_Case", "Lung_WT_Control"))

# DKO control vs DKO case
DKO_con_VS_DKO_case <- results(counts_dds_DESeq, contrast = c("condition", "Lung_DKO_Case", "Lung_DKO_Control"))

# WT case vs DKO case
WT_case_VS_DKO_case <- results(counts_dds_DESeq, contrast = c("condition", "Lung_DKO_Case", "Lung_WT_Case"))

# WT control vs DKO control
WT_con_VS_DKO_con <- results(counts_dds_DESeq, contrast = c("condition", "Lung_DKO_Control", "Lung_WT_Control"))

# --- Map ENSMUSG -> gene symbol (English comments, <-) ---
lung_map <- read.csv("Lung modules.csv", header = FALSE, stringsAsFactors = FALSE)
colnames(lung_map) <- c("ensembl_id", "gene_name", "module")

ens2sym <- setNames(lung_map$gene_name, lung_map$ensembl_id)

# Helper: return a label vector matching rownames(res)
get_labels <- function(res, map_vec) {
  rn <- rownames(res)
  if (all(grepl("^ENSMUSG", rn))) {
    lab <- unname(map_vec[rn])
    lab[is.na(lab)] <- rn[is.na(lab)]  # fallback if unmapped
    return(lab)
  } else {
    return(rn) # already gene symbols (or something else readable)
  }
}

# --- Define the EXACT genes you want to label ---
# the gene labels are exctracted from the output of the script "Important_genes.R".

lab_WT <- c("Tap1","Stat2","Irf7","Gbp5","Gbp2","Ifit2","Psmb8","Zbp1","Gbp10","Gbp3")

lab_DKO <- c("Oasl1","Oas2","Oas3","Cxcl10","Mx1","Gbp10","Cxcl9","Rsad2","Acta1","Krt13")

lab_WTca_vs_DKOca <- c("Tap1","Stat2","Irf7","Ifit1","Gbp5","Gbp2","Gbp3","Oas2","Ifit2","Oasl1")

# --- Build volcanoes with mapping + selectLab ---

# setting up a plot theme
panel_theme_volcano <- theme_minimal(base_size = 20) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"   # ok für volcano (wenn’s passt)
  )

# WT
Volcano_WT <- EnhancedVolcano(
  WT_con_VS_WT_case,
  lab = get_labels(WT_con_VS_WT_case, ens2sym),
  x = "log2FoldChange",
  y = "padj",
  title = NULL,
  selectLab = lab_WT,
  drawConnectors = TRUE,
  max.overlaps = Inf
) + panel_theme_volcano

# DKO 
Volcano_DKO <- EnhancedVolcano(
  DKO_con_VS_DKO_case,
  lab = get_labels(DKO_con_VS_DKO_case, ens2sym),
  x = "log2FoldChange",
  y = "padj",
  title = NULL,
  selectLab = lab_DKO,
  drawConnectors = TRUE,
  max.overlaps = Inf
) + panel_theme_volcano

# WT case vs DKO case
Volcano_WTca_vs_DKOca <- EnhancedVolcano(
  WT_case_VS_DKO_case,
  lab = get_labels(WT_case_VS_DKO_case, ens2sym),
  x = "log2FoldChange",
  y = "padj",
  title = NULL,
  selectLab = lab_WTca_vs_DKOca,
  drawConnectors = TRUE,
  max.overlaps = Inf
) + panel_theme_volcano

Volcano_WT
Volcano_DKO
Volcano_WTca_vs_DKOca
#---------------------------------------------------------------

# DE genes (padj < 0.05)
#extracting numbers for DKO
res_DKO <- as.data.frame(DKO_con_VS_DKO_case)
de_DKO <- res_DKO[!is.na(res_DKO$padj) & res_DKO$padj < 0.05, ]

#extracting numbers for case vs case
res_WT_DKO <- as.data.frame(WT_case_VS_DKO_case)
de_WT_DKO <- res_WT_DKO[!is.na(res_WT_DKO$padj) & res_WT_DKO$padj < 0.05, ]

#extracting numbers for control vs control
res_WT_DKO_con <- as.data.frame(WT_con_VS_DKO_con)
de_WT_DKO_con <- res_WT_DKO_con[!is.na(res_WT_DKO_con$padj) & res_WT_DKO_con$padj < 0.05, ]

#extracting numbers for control vs control
res_WT_DKO_con <- as.data.frame(WT_con_VS_DKO_con)
de_WT_DKO_con <- res_WT_DKO_con[!is.na(res_WT_DKO_con$padj) & res_WT_DKO_con$padj < 0.05, ]

#creating summary table
summary_DE <- data.frame(
  Comparison = c("WT Control vs Case", "DKO Control vs Case", "WT Case vs DKO Case", "WT Control vs DKO Control"),
  DE_genes = c(nrow(de_WT), nrow(de_DKO), nrow(de_WT_DKO), nrow(de_WT_DKO_con)),
  Upregulated = c(
    sum(de_WT$log2FoldChange > 0),
    sum(de_DKO$log2FoldChange > 0),
    sum(de_WT_DKO$log2FoldChange > 0),
    sum(de_WT_DKO_con$log2FoldChange > 0)
  ),
  Downregulated = c(
    sum(de_WT$log2FoldChange < 0),
    sum(de_DKO$log2FoldChange < 0),
    sum(de_WT_DKO$log2FoldChange < 0),
    sum(de_WT_DKO_con$log2FoldChange < 0)
  )
)
summary_DE