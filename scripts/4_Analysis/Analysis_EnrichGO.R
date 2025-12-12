
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

# Alle Gene, die in DESeq2 getestet wurden (Universe)
universe_genes <- rownames(counts_dds_DESeq)
length(universe_genes)


# WT: Lung_WT_Control vs Lung_WT_Case
res_WT <- WT_con_VS_WT_case

# Signifikante Gene (padj < 0.05, keine NAs)
res_WT_sig <- res_WT[!is.na(res_WT$padj) & res_WT$padj < 0.05, ]

de_genes_WT <- rownames(res_WT_sig)  # Ensembl-IDs der DE-Gene
length(de_genes_WT)

head(de_genes_WT)

ego_WT_BP <- enrichGO(
  gene          = de_genes_WT,      # DE-Gene
  universe      = universe_genes,   # alle getesteten Gene
  OrgDb         = org.Mm.eg.db,     # Maus-Annotation
  keyType       = "ENSEMBL",        # Format der IDs
  ont           = "BP",             # "BP", "MF", "CC" oder "ALL"
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE              # mappt Ensembl -> Gene Symbols
)

ego_WT_BP

#dotplot of the GO analysis
dotplot(ego_WT_BP, showCategory = 10) +
  ggtitle("GO BP enrichment: WT Control vs WT Case (padj < 0.05)")

#barplot of the GO analysis
barplot(ego_WT_BP, showCategory = 10) +
  ggtitle("GO BP enrichment: WT Control vs WT Case")

# Sort by adjusted p-value (most significant first)
ego_res <- ego_WT_BP@result[order(ego_WT_BP@result$p.adjust), ]

# Keep only top 50 (adjust if you want 100)
ego_WT_BP_top <- ego_WT_BP
ego_WT_BP_top@result <- ego_res[1:50, ]

#simplify GO term matrix to filter out redundancy
ego_WT_BP_s <- simplify(
  ego_WT_BP_top,
  measure = "Wang",
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

barplot(ego_WT_BP_s, showCategory = 10) +
  ggtitle("GO BP enrichment (simplified): WT Control vs WT Case")

