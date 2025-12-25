
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

#data for differential gene analysis
WT_con_VS_WT_case <- results(counts_dds_DESeq, contrast = c("condition", "Lung_WT_Control", "Lung_WT_Case"))

DKO_con_VS_DKO_case <- results(counts_dds_DESeq, contrast = c("condition", "Lung_DKO_Control", "Lung_DKO_Case"))

WT_con_VS_DKO_case <- results(counts_dds_DESeq, contrast = c("condition", "Lung_WT_Control", "Lung_DKO_Case"))


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

# --- Define the EXACT genes you want to label (from your Table 2 screenshot) ---
lab_WT <- c("Tap1","Stat2","Irf7","Gbp5","Gbp2","Ifit2","Psmb8","Zbp1","Gbp10","Gbp3")

lab_DKO <- c("Oasl1","Oas2","Oas3","Cxcl10","Mx1","Gbp10","Cxcl9","Rsad2","Acta1","Krt13")

lab_WTca_vs_DKOca <- c("Tap1","Stat2","Irf7","Ifit1","Gbp5","Gbp2","Gbp3","Oas2","Ifit2","Oasl1")

# --- Build volcanoes with mapping + selectLab ---
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
#extracting numbers
res_WT <- as.data.frame(WT_con_VS_WT_case)
de_WT <- res_WT[!is.na(res_WT$padj) & res_WT$padj < 0.05, ]

res_DKO <- as.data.frame(DKO_con_VS_DKO_case)
de_DKO <- res_DKO[!is.na(res_DKO$padj) & res_DKO$padj < 0.05, ]

#creating summary table
summary_DE <- data.frame(
  Comparison = c("WT Control vs Case", "DKO Control vs Case"),
  DE_genes = c(nrow(de_WT), nrow(de_DKO)),
  Upregulated = c(
    sum(de_WT$log2FoldChange > 0),
    sum(de_DKO$log2FoldChange > 0)
  ),
  Downregulated = c(
    sum(de_WT$log2FoldChange < 0),
    sum(de_DKO$log2FoldChange < 0)
  )
)

summary_DE
