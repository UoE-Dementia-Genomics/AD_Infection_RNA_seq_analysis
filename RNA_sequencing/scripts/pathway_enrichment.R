########################
# @author: Giulia Pegoraro
#
# Over-representation analysis (ORA) on WGCNA modules
# and limma / ANOVA-derived gene sets
#
# This script:
#  - Provides helper functions to strip Ensembl version suffixes
#    and convert Ensembl IDs via Ensembl (biomaRt)
#  - Performs GO BP enrichment (clusterProfiler::enrichGO) with
#    background universe and simplification by semantic similarity
#  - Performs KEGG pathway enrichment (clusterProfiler::enrichKEGG)
#  - Runs ORA for:
#      * limma contrasts (AD vs controls, with/without infection)
#      * WGCNA modules (all colours, plus specific tan/red/green modules)
#      * WGCNA hub genes
#      * Additional t-test / p-value based gene sets
#  - Visualises results with barplots and combined GO+KEGG panels
#  - Explores MSigDB-based enrichment (C2, C7, etc.) using enricher
########################

####################################################################
########################## LIBRARIES ###############################
####################################################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(purrr)
library(ggplot2)
library(cowplot)
library(msigdbr)
library(biomaRt)

####################################################################
########################## FUNCTIONS ###############################
####################################################################

############################################
# Function: removeVersion
#
# Description:
#   Remove the version number from Ensembl gene IDs, e.g.
#   "ENSG00000123456.7" -> "ENSG00000123456".
#
# Arguments:
#   genes : character vector of Ensembl IDs with version.
#
# Returns:
#   list of Ensembl IDs without version (one element per input).
############################################
removeVersion <- function(genes){
  gene_list <- strsplit((genes),"\\.")
  gene_list <- lapply(seq_along(gene_list), function(x) gene_list[[x]][[1]])
}

############################################
# Function: convertFromEnsembl
#
# Description:
#   Convert Ensembl gene IDs with version to a specified attribute
#   (e.g. HGNC symbol, Entrez ID) using Ensembl via biomaRt.
#
# Arguments:
#   targetList          : character vector of Ensembl IDs with version.
#   conversionAttribute : character, Ensembl attribute to retrieve
#                         (default: "hgnc_symbol").
#
# Returns:
#   data frame with columns:
#     ensembl_gene_id_version
#     <conversionAttribute>
############################################
convertFromEnsembl <- function(targetList, conversionAttribute = "hgnc_symbol"){
  ensembl <- useEnsembl(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl",
    version = 104
  )
  
  targetList <- getBM(
    attributes = c("ensembl_gene_id_version", conversionAttribute),
    filters    = "ensembl_gene_id_version",
    values     = targetList,
    mart       = ensembl
  )
  
  return(targetList)
}

############################################
# Function: convertToSymbol
#
# Description:
#   Map Ensembl gene IDs (with version) to HGNC gene symbols using biomaRt.
#
# Arguments:
#   targetList : character vector of Ensembl gene IDs with version.
#
# Returns:
#   data frame with columns:
#     hgnc_symbol
#     ensembl_gene_id_version
############################################
convertToSymbol <- function(targetList){
  ensembl = useEnsembl(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl",
    version = 104
  )
  
  targetList <- getBM(
    attributes = c("hgnc_symbol", "ensembl_gene_id_version"),
    filters    = "ensembl_gene_id_version",
    values     = targetList,
    mart       = ensembl
  )
  
  return(targetList)
}
############################################
# Function: GO_enrichment
#
# Description:
#   Perform GO Biological Process (BP) enrichment using enrichGO
#   on a list of Ensembl IDs with version, then simplify results
#   based on semantic similarity and generate a barplot.
#
# Arguments:
#   genes       : character vector of Ensembl IDs with version (gene set).
#   title_plot  : character, plot title.
#   background  : character vector of background Ensembl IDs with version
#                 to use as universe.
#
# Returns:
#   list of length 3:
#     [[1]] simplified enrichResult object (upSimGO)
#     [[2]] ggplot barplot
#     [[3]] original enrichResult object (ego)
############################################
GO_enrichment <- function(genes, title_plot, background){
  symbols       <- unlist(removeVersion(genes))
  background_id <- unlist(removeVersion(background))
  
  ego <- enrichGO(
    gene          = symbols,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff  = 1,
    qvalueCutoff  = 0.05,
    minGSSize     = 5,
    keyType       = "ENSEMBL",
    universe      = background_id
  )
  
  upSimGO <- clusterProfiler::simplify(
    ego,
    cutoff     = 0.8,
    by         = "p.adjust",
    select_fun = min,
    measure    = "Wang"
  )
  
  bplot <- barplot(
    upSimGO,
    showCategory = 10,
    title        = title_plot,
    color        = "qvalue"
  ) +
    theme(
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      axis.title  = element_text(size = 24)
    )
  
  list(upSimGO, bplot, ego)
}

############################################
# Function: KEGG_enrichment
#
# Description:
#   Perform KEGG pathway enrichment using enrichKEGG on a list
#   of Ensembl IDs with version, converting to Entrez IDs via
#   biomaRt first, and generate a barplot.
#
# Arguments:
#   genes       : character vector of Ensembl IDs with version.
#   title_plot  : character, plot title.
#   background  : character vector of background Ensembl IDs with version
#                 to use as universe.
#
# Returns:
#   list of length 2:
#     [[1]] enrichResult object (kegen)
#     [[2]] ggplot barplot
############################################
KEGG_enrichment <- function(genes, title_plot, background){
  gene.df     <- convertFromEnsembl(genes,      "entrezgene_id")
  universe.df <- convertFromEnsembl(background, "entrezgene_id")
  
  kegen <- enrichKEGG(
    gene          = as.character(gene.df$entrezgene_id),
    organism      = "hsa",
    pvalueCutoff  = 1,
    qvalueCutoff  = 0.05,
    minGSSize     = 5,
    maxGSSize     = 500,
    pAdjustMethod = "fdr",
    universe      = as.character(universe.df$entrezgene_id)
  )
  
  bplot <- barplot(
    kegen,
    showCategory = 10,
    title        = title_plot,
    color        = "qvalue"
  ) +
    theme(
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      axis.title  = element_text(size = 24)
    )
  
  list(kegen, bplot)
}

############################################
# Function: GO_module
#
# Description:
#   Perform GO BP enrichment for each WGCNA module and save
#   barplots to disk, using a common background universe.
#
# Arguments:
#   module_df : data frame with at least:
#                 - gene_id : Ensembl IDs with version
#                 - colors  : module colour per gene
#   path      : character, directory where plots will be saved.
#   background: character vector of background Ensembl IDs (no version)
#               to use as universe.
#
# Returns:
#   Invisibly returns list of enrichResult objects per module colour.
############################################
GO_module <- function(module_df, path, background){
  col_enrich_resp_inf_tot <- lapply(unique(module_df$colors), function(color){
    gene_list <- removeVersion(module_df[module_df$colors == color, "gene_id"])
    
    ego <- enrichGO(
      gene          = gene_list,
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      pAdjustMethod = "fdr",
      pvalueCutoff  = 1,
      qvalueCutoff  = 1,
      minGSSize     = 5,
      maxGSSize     = 500,
      keyType       = "ENSEMBL",
      universe      = background
    )
    
    upSimGO <- clusterProfiler::simplify(
      ego,
      cutoff     = 0.8,
      by         = "p.adjust",
      select_fun = min,
      measure    = "Wang"
    )
    
    upSimGO
  })
  
  names(col_enrich_resp_inf_tot) <- unique(module_df$colors)
  col_enrich_resp_inf_tot_filt   <- col_enrich_resp_inf_tot %>% keep(~ nrow(.) != 0)
  
  lapply(names(col_enrich_resp_inf_tot_filt), function(i){
    file_path <- paste0(path, i, "_enrich.tiff")
    barplot(col_enrich_resp_inf_tot_filt[[i]], showCategory = 10) +
      theme(
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title  = element_text(size = 14)
      )
    ggsave(file_path, width = 8, height = 8)
  })
  
  invisible(col_enrich_resp_inf_tot_filt)
}

############################################
# Function: KEGG_module
#
# Description:
#   Perform KEGG enrichment for each WGCNA module and save
#   barplots to disk, using a common background.
#
# Arguments:
#   module_df : data frame with at least:
#                 - gene_id : Ensembl IDs with version
#                 - colors  : module colour per gene
#   path      : character, directory where plots will be saved.
#   background: vector of background Entrez IDs (if used directly).
#
# Returns:
#   Invisibly returns list of enrichResult objects per module colour.
############################################
KEGG_module <- function(module_df, path, background){
  col_enrich_kegg <- lapply(unique(module_df$colors), function(color){
    gene_list <- strsplit(module_df[module_df$colors == color, "gene_id"], "\\.")
    gene_list <- lapply(seq_along(gene_list), function(x) gene_list[[x]][[1]])
    
    gene.df <- bitr(
      gene_list,
      fromType = "ENSEMBL",
      toType   = "ENTREZID",
      OrgDb    = org.Hs.eg.db
    )
    
    kegen <- enrichKEGG(
      gene         = gene.df$ENTREZID,
      organism     = "hsa",
      pvalueCutoff = 0.05,
      universe     = background
    )
    
    kegen
  })
  
  names(col_enrich_kegg) <- unique(module_df$colors)
  col_enrich_kegg        <- col_enrich_kegg %>% keep(~ nrow(.) != 0)
  
  lapply(names(col_enrich_kegg), function(i){
    file_path <- paste0(path, i, "_enrich.png")
    barplot(col_enrich_kegg[[i]], showCategory = 10)
    ggsave(file_path)
  })
  
  invisible(col_enrich_kegg)
}

############################################
# Function: GO_ARCHS4
#
# Description:
#   GO BP enrichment for T-test / ARCHS4-style gene lists
#   using gene symbols, returning a simplified enrichResult
#   and barplot.
#
# Arguments:
#   symbols : character vector of HGNC gene symbols.
#   title   : character, plot title.
#
# Returns:
#   enrichResult and barplot (invisibly via side-effect).
############################################
GO_ARCHS4 <- function(symbols, title){
  ego <- enrichGO(
    gene          = symbols,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    keyType       = "SYMBOL"
  )
  
  upSimGO <- clusterProfiler::simplify(
    ego,
    cutoff     = 0.8,
    by         = "p.adjust",
    select_fun = min,
    measure    = "Wang"
  )
  
  barplot(upSimGO, showCategory = 10) + labs(title = title)
}

############################################
# Function: KEGG_ARCHS4
#
# Description:
#   KEGG enrichment for T-test / ARCHS4-style gene lists
#   using gene symbols.
#
# Arguments:
#   symbols : character vector of HGNC gene symbols.
#   title   : character, used in filenames and plot title.
#
# Returns:
#   None (plots saved to disk; kegen and barplot not returned).
############################################
KEGG_ARCHS4 <- function(symbols, title){
  gene.df <- bitr(
    symbols,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )
  print(gene.df)
  
  kegen <- enrichKEGG(
    gene          = gene.df$ENTREZID,
    organism      = "hsa",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    pAdjustMethod = "fdr"
  )
  
  file_path_plot <- paste0(
    "C:\\Users\\gp487\\OneDrive - University of Exeter\\AD_infection\\Analyses\\RNA_sequencing\\output\\figures\\enrichment_plot\\KEGG\\respiratory_infection\\cell",
    title,
    ".png"
  )
  
  barplot(kegen, showCategory = 10) + labs(title = title)
  ggsave(
    "C:\\Users\\gp487\\OneDrive - University of Exeter\\AD_infection\\Analyses\\RNA_sequencing\\output\\figures\\enrichment_plot\\KEGG\\respiratory_infection\\cell"
  )
}

#########################################################################
################################## MAIN #################################
#########################################################################

setwd("C:\\Users\\gp487\\OneDrive - University of Exeter\\AD_infection\\Analyses\\RNA_sequencing")

##### Import background and hub genes ###################################

background_genes <- convertToSymbol(rownames(gene_expr_resp_inf))

hubgenesYellow <- read.table(
  "C:/Users/gp487/OneDrive - University of Exeter/AD_infection/Analyses/RNA_sequencing/output/data/hub_genes_yellow.txt",
  header = TRUE
)

##### Significant values from limma model ################################

genes_57 <- sign_AD_vs_Control_inf[
  !sign_AD_vs_Control_inf$ensembl_gene_id_version %in% sign_AD_vs_Control_no$ensembl_gene_id_version,
]

genes_57_GO <- GO_enrichment(
  genes_57$ensembl_gene_id_version,
  "AD vs. control with infection",
  background_genes$ensembl_gene_id_version
)

sign_AD_vs_Control_inf_ensembl <- sign_AD_vs_Control_inf$Row.names

AD_vs_Control_inf_GO <- GO_enrichment(
  sign_AD_vs_Control_inf$ensembl_gene_id_version,
  "AD vs. control with infection",
  background_genes$ensembl_gene_id_version
)
AD_vs_Control_inf_GO_down <- GO_enrichment(
  sign_AD_vs_Control_inf[sign_AD_vs_Control_inf$logFC < 0, "ensembl_gene_id_version"],
  "",
  background_genes$ensembl_gene_id_version
)
AD_vs_Control_inf_GO_up <- GO_enrichment(
  sign_AD_vs_Control_inf[sign_AD_vs_Control_inf$logFC > 0, "ensembl_gene_id_version"],
  "",
  background_genes$ensembl_gene_id_version
)

AD_vs_Control_no_GO <- GO_enrichment(
  sign_AD_vs_Control_no$ensembl_gene_id_version,
  "",
  background_genes$ensembl_gene_id_version
)
AD_vs_Control_no_GO_down <- GO_enrichment(
  sign_AD_vs_Control_no[sign_AD_vs_Control_no$logFC < 0, "ensembl_gene_id_version"],
  "",
  background_genes$ensembl_gene_id_version
)
AD_vs_Control_no_GO_up <- GO_enrichment(
  sign_AD_vs_Control_no[sign_AD_vs_Control_no$logFC > 0, "ensembl_gene_id_version"],
  "",
  background_genes$ensembl_gene_id_version
)

AD_vs_Control_inf_KEGG <- KEGG_enrichment(
  sign_AD_vs_Control_inf$ensembl_gene_id_version,
  "AD vs. control with infection",
  background_genes$ensembl_gene_id_version
)
genes_57_KEGG <- KEGG_enrichment(
  genes_57$ensembl_gene_id_version,
  "AD vs. control with infection",
  background_genes$ensembl_gene_id_version
)

AD_vs_Control_noinf_KEGG <- KEGG_enrichment(
  sign_AD_vs_Control_no$ensembl_gene_id_version,
  "",
  background_genes$ensembl_gene_id_version
)
AD_vs_Control_no_KEGG_down <- KEGG_enrichment(
  sign_AD_vs_Control_no[sign_AD_vs_Control_no$logFC < 0, "ensembl_gene_id_version"],
  "",
  background_genes$ensembl_gene_id_version
)
AD_vs_Control_no_KEGG_up <- KEGG_enrichment(
  sign_AD_vs_Control_no[sign_AD_vs_Control_no$logFC > 0, "ensembl_gene_id_version"],
  "",
  background_genes$ensembl_gene_id_version
)

tiff("output\\figures\\enrichment_plot\\GO\\respiratory_infection\\AD_vs_Control_no_GO.tiff",
     height = 700, width = 650)
AD_vs_Control_no_GO[[2]]
dev.off()

tiff("output\\figures\\enrichment_plot\\KEGG\\respiratory_infection\\AD_vs_Control_no_KEGG.tiff",
     height = 700, width = 650)
AD_vs_Control_noinf_KEGG[[2]]
dev.off()

##### MSigDB-based ORA for limma results #################################

msigdb_genesets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP")
msigdbr_t2g     <- msigdb_genesets %>%
  dplyr::distinct(gs_name, gene_symbol) %>%
  as.data.frame()

ora_results_c2_cp <- enricher(
  gene       = results_combined_annot[
    results_combined_annot$pvalue_AD_vs_Control_inf < 0.05 &
      abs(results_combined_annot$logFC_AD_vs_Control_inf) > log2(1.1),
    "ensembl_gene_id_version"
  ],
  TERM2GENE  = msigdbr_t2g,
  pvalueCutoff = 0.05
)
barplot(ora_results_c2_cp, showCategory = 10, title = "Top 10 Enriched Pathways in C2.CP")

##### WGCNA modules: GO + KEGG ##########################################

resp_inf_path_cell_GO   <- "output\\figures\\enrichment_plot\\GO\\respiratory_infection\\"
resp_inf_path_cell_KEGG <- "output\\figures\\enrichment_plot\\KEGG\\respiratory_infection\\"

GO_module(
  module_df_cell_no_fetal[[1]],
  resp_inf_path_cell_GO
)
KEGG_module(
  module_df_cell_no_fetal[[1]],
  resp_inf_path_cell_KEGG
)

GO_tan <- GO_enrichment(
  module_df_cell_no_fetal[[1]][module_df_cell_no_fetal[[1]]$colors == "tan", "gene_id"],
  "",
  background_genes$ensembl_gene_id_version
)
KEGG_tan <- KEGG_enrichment(
  module_df_cell_no_fetal[[1]][module_df_cell_no_fetal[[1]]$colors == "tan", "gene_id"],
  "",
  background_genes$ensembl_gene_id_version
)

GO_red <- GO_enrichment(
  module_df_cell_no_fetal[[1]][module_df_cell_no_fetal[[1]]$colors == "red", "gene_id"],
  "",
  background_genes$ensembl_gene_id_version
)
KEGG_red <- KEGG_enrichment(
  module_df_cell_no_fetal[[1]][module_df_cell_no_fetal[[1]]$colors == "red", "gene_id"],
  "",
  background_genes$ensembl_gene_id_version
)

tiff("output\\figures\\enrichment_plot\\GO\\respiratory_infection\\tan_GO.tiff",
     height = 700, width = 650)
GO_tan[[2]]
dev.off()

tiff("output\\figures\\enrichment_plot\\KEGG\\respiratory_infection\\red_KEGG.tiff",
     height = 700, width = 650)
KEGG_red[[2]]
dev.off()

##### WGCNA hub genes ###################################################

hubgenesYellow_plot      <- GO_enrichment(
  hubgenesYellow$ensembl_gene_id_version,
  "",
  background_genes$ensembl_gene_id_version
)
hubgenesYellow_plot_KEGG <- KEGG_enrichment(
  hubgenesYellow$ensembl_gene_id_version,
  "",
  background_genes$ensembl_gene_id_version
)

tiff("output\\figures\\enrichment_plot\\GO\\respiratory_infection\\hubgenesYellow_GO.tiff",
     height = 700, width = 650)
hubgenesYellow_plot[[2]]
dev.off()

tiff("output\\figures\\enrichment_plot\\KEGG\\respiratory_infection\\hubgenesYellow_KEGG.tiff",
     height = 700, width = 650)
hubgenesYellow_plot_KEGG[[2]]
dev.off()

##### Example: GO enrichment on green module symbols ####################

ego <- enrichGO(
  gene          = symbols_green$hgnc_symbol,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  keyType       = "SYMBOL"
)

KEGG_green[[2]] <- KEGG_green[[2]] +
  theme(
    axis.title = element_text(size = 22),
    axis.text  = element_text(size = 20)
  )

GO_green[[2]] <- GO_green[[2]] +
  theme(
    axis.title = element_text(size = 22),
    axis.text  = element_text(size = 20)
  )

combined_plot_GO <- plot_grid(
  KEGG_green[[2]],
  GO_green[[2]],
  ncol      = 2,
  labels    = LETTERS[1:2],
  label_size = 24
)

ggsave(
  "plots\\enrichment_plot\\WGCNA.tiff",
  combined_plot_GO, width = 26, height = 13
)

##### MSigDB C2 enrichment for green module (Entrez) ####################

msigdb_genesets      <- msigdbr(species = "Homo sapiens", category = "C2")
gene_list_green      <- strsplit(
  module_df_cell_no_fetal[[1]][module_df_cell_no_fetal[[1]]$colors == "green", "gene_id"],
  "\\."
)
gene_list_green      <- lapply(seq_along(gene_list_green), function(x) gene_list_green[[x]][[1]])
msigdbr_t2g_entrez   <- msigdb_genesets %>%
  dplyr::distinct(gs_name, entrez_gene) %>%
  as.data.frame()
gene.df_green        <- bitr(
  gene_list_green,
  fromType = "ENSEMBL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)
ora_results_c2_cp <- enricher(
  gene        = gene.df_green$ENTREZID,
  TERM2GENE   = msigdbr_t2g_entrez,
  pvalueCutoff = 0.05
)
barplot(ora_results_c2_cp, showCategory = 10, title = "Top 10 Enriched Pathways in C2.CP")

##### KEGG enrichment: ANOVA significant genes ##########################

gene.df_resp_inf <- bitr(
  p_values_AD_resp_inf[
    p_values_AD_resp_inf$`anova_p_resp_tot[, "AD:Infection"]` < 0.05,
    "genes"
  ],
  fromType = "ENSEMBL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

kegen_resp_inf <- enrichKEGG(
  gene         = gene.df_resp_inf$ENTREZID,
  organism     = "hsa",
  pvalueCutoff = 0.1
)
barplot(kegen_resp_inf, showCategory = 10)

##### Module-level KEGG enrichment (respiratory infection) ###############

resp_inf_path_cell_kegg <- "output\\figures\\enrichment_plot\\KEGG\\respiratory_infection\\cell\\"
KEGG_module(module_df_cell_no_fetal[[1]], resp_inf_path_cell_kegg)

##### Example: MSigDB Hallmark / C7 enrichment ##########################

m_t2g <- msigdbr(species = "Homo sapiens", category = "C7")

gene_list <- module_df_cell_no_fetal[[1]][
  module_df_cell_no_fetal[[1]]$colors == "green",
  "gene_id"
]
gene_list <- convertFromEnsembl(gene_list)
symbol    <- gene_list$hgnc_symbol
symbol    <- symbol[symbol != ""]

gene_df <- bitr(
  gene_list,
  fromType = "ENSEMBL",
  toType   = "SYMBOL",
  OrgDb    = org.Hs.eg.db
)
entrez_gene_list <- gene_df$ENTREZID

enrich_results <- enricher(entrez_gene_list, TERM2GENE = m_t2g)
head(enrich_results)
dotplot(enrich_results)
