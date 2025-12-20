
#At the moment this file operates with data created by the other file ("Analysis_DeSeq2"), that could be found in the same repo.
#This is bad practice and will be adjusted asap
#For the moment getting results and progressing in the project has higher priority.


#-------------------------------------------------------
#loading BiocManager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(BiocManager)

#loading packages
BiocManager::install("clusterProfiler", update = FALSE, ask = FALSE)
BiocManager::install("org.Mm.eg.db")
install.packages("enrichplot")

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

#-------------------------------------------------------


### 7. Overrepresentation analysis (GO)


# All genes, that were tested in DESeq2 analysisi
universe_genes <- rownames(counts_dds_DESeq)

# WT: Lung_WT_Control vs Lung_WT_Case
res_WT <- WT_con_VS_WT_case
res_DKO <- DKO_con_VS_DKO_case

# Significant genes (padj < 0.05, ni NAs)
res_WT_sig <- res_WT[!is.na(res_WT$padj) & res_WT$padj < 0.05, ]
res_DKO_sig <- res_DKO[!is.na(res_DKO$padj) & res_DKO$padj < 0.05, ]

de_genes_WT <- rownames(res_WT_sig)  # Ensembl-IDs der DE-Gene

de_genes_DKO <- rownames(res_DKO_sig)

ego_WT_BP <- enrichGO(
  gene          = de_genes_WT,      # DE genes
  universe      = universe_genes,   # all tested genes
  OrgDb         = org.Mm.eg.db,     # Mouse Annotation
  keyType       = "ENSEMBL",        # format of IDs
  ont           = "BP",             # "BP", "MF", "CC" or "ALL"
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE              # mapping Ensembl -> Gene Symbols
)
ego_DKO_BP <- enrichGO(
  gene          = de_genes_DKO,      # DE genes
  universe      = universe_genes,   # all tested genes
  OrgDb         = org.Mm.eg.db,     # Mouse Annotation
  keyType       = "ENSEMBL",        # format of IDs
  ont           = "BP",             # "BP", "MF", "CC" or "ALL"
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE              # mapping Ensembl -> Gene Symbols
)

# Sort by adjusted p-value (most significant first)
ego_res <- ego_WT_BP@result[order(ego_WT_BP@result$p.adjust), ]
ego_res_DKO <- ego_DKO_BP@result[order(ego_DKO_BP@result$p.adjust),]

# Keep only top 50 (adjust if you want 100)
ego_WT_BP_top <- ego_WT_BP
ego_WT_BP_top@result <- ego_res[1:50, ]

ego_DKO_BP_top <- ego_DKO_BP
ego_DKO_BP_top@result <- ego_res_DKO[1:50, ]

#simplify GO term matrix to filter out redundancy
ego_WT_BP_s <- simplify(
  ego_WT_BP_top,
  measure = "Wang",
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

ego_DKO_BP_s <- simplify(
  ego_DKO_BP_top,
  measure = "Wang",
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

barplot_ego_WT_BP_s <- barplot(ego_WT_BP_s, showCategory = 10) +
  panel_theme_go +
  theme(axis.text.y = element_text(size = 20))

barplot_ego_DKO_BP_s <- barplot(ego_DKO_BP_s, showCategory = 10) +
  panel_theme_go +
  theme(axis.text.y = element_text(size = 20))
