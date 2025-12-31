
#-------------------------------------------------------------------------
# performing an enrichGO analysis
#-------------------------------------------------------------------------

#At the moment this file operates with data created by the other file ("Analysis_DeSeq2"), that could be found in the same repo.

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


# All genes, that were tested in DESeq2 analysis
universe_genes <- rownames(counts_dds_DESeq)

# Getting the groups for comparison
res_WT <- WT_con_VS_WT_case
res_DKO <- DKO_con_VS_DKO_case
res_WT_DKO <- WT_case_VS_DKO_case
res_WT_DKO_con <- WT_con_VS_DKO_con

# Significant genes (padj < 0.05, ni NAs)
res_WT_sig <- res_WT[!is.na(res_WT$padj) & res_WT$padj < 0.05, ]
res_DKO_sig <- res_DKO[!is.na(res_DKO$padj) & res_DKO$padj < 0.05, ]
res_WT_DKO_sig <- res_WT_DKO[!is.na(res_WT_DKO$padj) & res_WT_DKO$padj < 0.05, ]
res_WT_DKO_con_sig <- res_WT_DKO_con[!is.na(res_WT_DKO_con$padj) & res_WT_DKO_con$padj < 0.05, ]

# row names/gene names of differentially expressed genes
de_genes_WT <- rownames(res_WT_sig)  # Ensembl-IDs der DE-Gene
de_genes_DKO <- rownames(res_DKO_sig)
de_genes_WT_DKO <- rownames(res_WT_DKO_sig)
de_genes_WT_DKO_con <- rownames(res_WT_DKO_con_sig)


#performing analysis for each group with BP parameter (Biological Process)

# WT case vs control
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

# DKO case vs control
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

#WT case vs DKO case
ego_WT_DKO_BP <- enrichGO(
  gene          = de_genes_WT_DKO,      # DE genes
  universe      = universe_genes,   # all tested genes
  OrgDb         = org.Mm.eg.db,     # Mouse Annotation
  keyType       = "ENSEMBL",        # format of IDs
  ont           = "BP",             # "BP", "MF", "CC" or "ALL"
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE              # mapping Ensembl -> Gene Symbols
)

# WT control vs DKO control
#This specific plot does not show up in the final report, since its interpretation does not support understanding of the data in the view of the author.
ego_WT_DKO_con_BP <- enrichGO(
  gene          = de_genes_WT_DKO_con,      # DE genes
  universe      = universe_genes,   # all tested genes
  OrgDb         = org.Mm.eg.db,     # Mouse Annotation
  keyType       = "ENSEMBL",        # format of IDs
  ont           = "BP",             # "BP"
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE              # mapping Ensembl -> Gene Symbols
)

# Sort by adjusted p-value (most significant first)
ego_res <- ego_WT_BP@result[order(ego_WT_BP@result$p.adjust), ]
ego_res_DKO <- ego_DKO_BP@result[order(ego_DKO_BP@result$p.adjust),]
ego_res_WT_DKO <- ego_WT_DKO_BP@result[order(ego_WT_DKO_BP@result$p.adjust),]
ego_res_WT_DKO_con <-ego_WT_DKO_con_BP@result[order(ego_WT_DKO_con_BP@result$p.adjust),]

# Keep only top 50 (adjust if you want 100)
#WT
ego_WT_BP_top <- ego_WT_BP
ego_WT_BP_top@result <- ego_res[1:50, ]
#DKO
ego_DKO_BP_top <- ego_DKO_BP
ego_DKO_BP_top@result <- ego_res_DKO[1:50, ]
#WT case vs DKO case
ego_DKO_WT_BP_top <- ego_WT_DKO_BP
ego_DKO_WT_BP_top@result <- ego_res_WT_DKO[1:50,]
#WT control vs DKO control
ego_DKO_WT_con_BP_top <- ego_WT_DKO_con_BP
ego_DKO_WT_con_BP_top@result <- ego_res_WT_DKO_con[1:50,]

#simplify GO term matrix to filter out redundancy
# WT
ego_WT_BP_s <- simplify(
  ego_WT_BP_top,
  measure = "Wang",
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)
# DKO
ego_DKO_BP_s <- simplify(
  ego_DKO_BP_top,
  measure = "Wang",
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)
# WT case vs DKO case
ego_WT_DKO_BP_s <- simplify(
  ego_DKO_WT_BP_top,
  measure = "Wang",
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)
# WT control vs DKO control
ego_WT_DKO_con_BP_s <- simplify(
  ego_DKO_WT_con_BP_top,
  measure = "Wang",
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

# barplots of simplified GO terms matrix

#panel theme for go terms
panel_theme_go <- theme_minimal(base_size = 20) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    panel.grid.minor = element_blank(),
    legend.position = "right",   # <-- wichtig
    legend.direction = "vertical"
  )

# WT
barplot_ego_WT_BP_s <- barplot(ego_WT_BP_s, showCategory = 10) +
  panel_theme_go +
  theme(axis.text.y = element_text(size = 20))
# DKO
barplot_ego_DKO_BP_s <- barplot(ego_DKO_BP_s, showCategory = 10) +
  panel_theme_go +
  theme(axis.text.y = element_text(size = 20))
# WT case vs DKO case
barplot_ego_WT_DKO_BP_s <- barplot(ego_WT_DKO_BP_s, showCategory = 10) +
  panel_theme_go +
  theme(axis.text.y = element_text(size = 20))
# WT control vs DKO control
barplot_ego_WT_DKO_con_BP_s <- barplot(ego_WT_DKO_con_BP_s, showCategory = 10) +
  panel_theme_go +
  theme(axis.text.y = element_text(size = 20))

