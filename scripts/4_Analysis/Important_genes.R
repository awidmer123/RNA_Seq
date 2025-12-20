
#setting up working directory
setwd("/home/andri/Desktop/UNIFR/2_Bioinf/RNA_Sequencing")

#loading necessary packages
library(dplyr)
library(knitr)
library(kableExtra)

#loading and modifying data
module_df <- read.table("Lung modules.csv",sep = ",")
colnames(module_df) <- c("ensembl_id", "Gene.name", "Module")

#------------------FUNCTIONS-----------------------------------------------

# 1. defining function that annotates modules
res_to_annot <- function(res_obj, module_df) {
  as.data.frame(res_obj) %>%
    tibble::rownames_to_column("ensembl_id") %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%     
    left_join(module_df, by = "ensembl_id")
}

# 2. defining function that picks interesting candidates
pick_candidates <- function(df, label, padj_cutoff = 0.05, lfc_cutoff = 1, top_n = 10) {
  df %>%
    mutate(is_ifn = gene_name %in% ifn_immune_genes) %>%
    filter(padj < padj_cutoff, abs(log2FoldChange) >= lfc_cutoff) %>%
    arrange(desc(is_ifn), padj, desc(abs(log2FoldChange))) %>%
    mutate(contrast = label) %>%
    select(contrast, ensembl_id, gene_name, module, baseMean, log2FoldChange, padj, is_ifn) %>%
    head(top_n)
}
#--------------------------------------------------------------------------

#transmute module dataframe
module_df <- module_df %>%
  transmute(
    ensembl_id = ensembl_id,
    gene_name = Gene.name,
    module = Module
  )

#using "res_to_annot" function
wt_annot  <- res_to_annot(WT_con_VS_WT_case, module_df)
dko_annot <- res_to_annot(DKO_con_VS_DKO_case, module_df)


#classical from literature derived genes that are involved in immune response and ifn related
ifn_immune_genes <- c(
  "Ifit1","Ifit2","Ifit3","Isg15","Rsad2","Mx1","Mx2",
  "Oas1a","Oas1g","Oas2","Oas3","Oasl1","Oasl2",
  "Irf7","Irf9","Stat1","Stat2","Usp18","Ddx58","Ifih1","Zbp1",
  "Cxcl9","Cxcl10","Tap1","Tap2","B2m","Psmb8","Psmb9","Nlrc5",
  "Gbp2","Gbp3","Gbp5","Gbp6","Gbp7","Gbp10",
  "Il1b"
)


#using the pick_candidates function to select for intersting hits
wt_hits  <- pick_candidates(wt_annot,  "WT: infected vs control")
dko_hits <- pick_candidates(dko_annot, "DKO: infected vs control")

wt_hits
dko_hits

dko_hits |>
  kable(
    format = "latex",
    booktabs = TRUE,
    caption = "Selected significantly differentially expressed genes in DKO case vs control."
  ) |>
  kable_styling(
    latex_options = c("hold_position", "striped", "scale_down"),
    font_size = 9
  )

wt_hits |>
  kable(
    format = "latex",
    booktabs = TRUE,
    caption = "Selected significantly differentially expressed genes in WT case vs control."
  ) |>
  kable_styling(
    latex_options = c("hold_position", "striped", "scale_down"),
    font_size = 9
  )




