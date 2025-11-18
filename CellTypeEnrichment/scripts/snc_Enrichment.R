#########################################################################################
################ ROSMAP snRNAseq Processing for EWCE Analysis ###########################
#########################################################################################
# @author: Giulia Pegoraro
#
# Cell-type enrichment analysis using EWCE on ROSMAP snRNA-seq data.
#
# This script:
#  - Loads ROSMAP cell-type reference (PFC snRNA-seq)
#  - Loads gene sets from WGCNA modules, limma DE results, ANOVA hits and miRNA modules
#  - Converts Ensembl IDs to HGNC symbols via Ensembl (biomaRt)
#  - Runs EWCE bootstrap enrichment tests at broad (annotation level 1) and fine
#    (annotation level 2) cell-type resolutions
#  - Performs analogous analyses with mouse reference (quality control / sanity check)
#  - Aggregates EWCE results across modules and builds cell-type × module heatmaps
#    for the full, tan, and red module sets
#########################################################################################

#########################################################################################
################################### LIBRARIES ###########################################
#########################################################################################

if (!requireNamespace("DropletUtils", quietly = TRUE))
  BiocManager::install("DropletUtils", dependencies = TRUE)
if (!requireNamespace("EWCE", quietly = TRUE))
  BiocManager::install("EWCE")
if (!requireNamespace("ewceData", quietly = TRUE))
  BiocManager::install("ewceData")
if (!requireNamespace("SingleCellExperiment", quietly = TRUE))
  BiocManager::install("SingleCellExperiment")
if (!requireNamespace("Seurat", quietly = TRUE))
  BiocManager::install("Seurat")
if (!requireNamespace("cowplot", quietly = TRUE))
  BiocManager::install("cowplot")
if (!requireNamespace("som", quietly = TRUE))
  BiocManager::install("som")
if (!requireNamespace("som", quietly = TRUE))
  BiocManager::install("som")
if (!requireNamespace("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")
if (!requireNamespace("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

library(DropletUtils)
library(EWCE)
library(ewceData)
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(som)
library(biomaRt)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(stringr)

#########################################################################################
################################### FUNCTIONS ###########################################
#########################################################################################

############################################
# Function: calculateAnnotation1
#
# Description:
#   Run EWCE bootstrap enrichment at annotation level 1 (broad cell types)
#   and build a barplot of SD from mean per broad cell class.
#
# Arguments:
#   ctd         : EWCE cell-type data list (human reference).
#   targetGenes : character vector of HGNC symbols to test.
#   title       : character, plot title.
#
# Returns:
#   list of length 2:
#     [[1]] EWCE bootstrap_enrichment_test result object (unconditional_results)
#     [[2]] ggplot barplot of SD from mean by broad cell type.
############################################
calculateAnnotation1 <- function(ctd, targetGenes, title){
  unconditional_results <- EWCE::bootstrap_enrichment_test(
    sct_data         = ctd,
    hits             = targetGenes,
    sctSpecies       = "human",
    genelistSpecies  = "human",
    sctSpecies_origin = "human",
    geneSizeControl  = TRUE,
    mtc_method       = "BH",
    reps             = 100000,
    annotLevel       = 1,
    no_cores         = 1
  )
  
  Res1 <- unconditional_results$results
  Res1 <- Res1[-(which(colnames(Res1) == "CellType"))]
  
  cell_df <- data.frame(
    row.names = c("Oli","Opc","Per","End","Mic","Ast","In","Ex"),
    CellType  = c("Oligodendrocytes","OPC","Pericytes","Endothelia",
                  "Microglia","Astrocytes","Inhibitory Neurons","Excitatory Neurons")
  )
  Res1 <- merge(Res1, cell_df, by = "row.names")
  
  # Set negative SD from mean to zero for plotting
  Res1[Res1$sd_from_mean < 0, "sd_from_mean"] <- 0
  
  # Add significance stars based on adjusted q
  Res1[Res1$q < 0.001,                       "sig"] <- "***"
  Res1[Res1$q < 0.01  & Res1$q > 0.001,      "sig"] <- "**"
  Res1[Res1$q < 0.05  & Res1$q > 0.01,       "sig"] <- "*"
  Res1[Res1$q > 0.05,                        "sig"] <- ""
  
  rosmapBroad <- ggplot(Res1, aes(x = as.factor(CellType), y = abs(sd_from_mean))) +
    geom_bar(stat = "identity", color = "black", lwd = 1) +
    geom_text(aes(label = sig, y = sd_from_mean + 0.05), size = 10) +
    scale_fill_distiller(direction = -1) +
    ylab("SD from mean") +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab(NULL) +
    ggtitle(title)
  
  return(list(unconditional_results, rosmapBroad))
}

############################################
# Function: calculateAnnotation2
#
# Description:
#   Run EWCE bootstrap enrichment at annotation level 2 (fine cell types)
#   and build a barplot of SD from mean with significance stars.
#
# Arguments:
#   ctd         : EWCE cell-type data list (human reference).
#   targetGenes : character vector of HGNC symbols to test.
#   title       : character, plot title.
#
# Returns:
#   list of length 2:
#     [[1]] data.frame of EWCE results (Res2)
#     [[2]] ggplot barplot of SD from mean by fine cell type.
############################################
calculateAnnotation2 <- function(ctd, targetGenes, title){
  unconditional_results2 <- EWCE::bootstrap_enrichment_test(
    sct_data         = ctd,
    hits             = targetGenes,
    sctSpecies       = "human",
    genelistSpecies  = "human",
    sctSpecies_origin = "human",
    geneSizeControl  = TRUE,
    reps             = 100000,
    mtc_method       = "BH",
    annotLevel       = 2,
    no_cores         = 1
  )
  
  Res2 <- unconditional_results2$results
  
  # Add significance stars based on adjusted q
  Res2[Res2$q < 0.001,                       "sig"] <- "***"
  Res2[Res2$q < 0.01  & Res2$q > 0.001,      "sig"] <- "**"
  Res2[Res2$q < 0.05  & Res2$q > 0.01,       "sig"] <- "*"
  Res2[Res2$q > 0.05,                        "sig"] <- ""
  
  rosmapBroad <- ggplot(Res2, aes(x = as.factor(CellType), y = sd_from_mean)) +
    geom_bar(stat = "identity", color = "black", lwd = 1, fill = "gray") +
    geom_text(
      aes(label = sig,
          y = ifelse(sd_from_mean > 0, sd_from_mean + 0.05, sd_from_mean - 0.2)),
      size = 12
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(Res2$sd_from_mean) + 1)) +
    ylab("SD from mean") +
    theme_cowplot() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 18),
      axis.text.y  = element_text(size = 16),
      axis.title.y = element_text(size = 18),
      plot.title   = element_text(size = 18, face = "bold")
    ) +
    xlab(NULL) +
    ggtitle(title)
  
  return(list(Res2, rosmapBroad))
}

############################################
# Function: convertToSymbol
#
# Description:
#   Convert Ensembl IDs with version to HGNC symbols using Ensembl (biomaRt).
#
# Arguments:
#   targetList : character vector of Ensembl IDs with version.
#
# Returns:
#   data frame with columns:
#     ensembl_gene_id_version
#     hgnc_symbol
############################################
convertToSymbol <- function(targetList){
  ensembl <- useEnsembl(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl",
    version = 104
  )
  
  targetList <- getBM(
    attributes = c("ensembl_gene_id_version", "hgnc_symbol"),
    filters    = "ensembl_gene_id_version",
    values     = targetList,
    mart       = ensembl
  )
  
  return(targetList)
}

############################################
# Function: MouseDatasetControl
#
# Description:
#   Run EWCE enrichment using the default mouse cortical reference (ewceData::ctd)
#   as a cross-species control for human gene lists.
#
# Arguments:
#   targetGenes : character vector of HGNC symbols to test.
#
# Returns:
#   list of length 2:
#     [[1]] EWCE results data.frame (unconditional_results_Mouse$results)
#     [[2]] EWCE default plot (ewce_plot)
############################################
MouseDatasetControl <- function(targetGenes){
  MouseCTD <- ewceData::ctd()
  
  unconditional_results_Mouse <- EWCE::bootstrap_enrichment_test(
    sct_data        = MouseCTD,
    hits            = targetGenes,
    sctSpecies      = "mouse",
    genelistSpecies = "human",
    geneSizeControl = TRUE,
    reps            = 100000,
    annotLevel      = 1,
    no_cores        = 1
  )
  
  ResMice <- unconditional_results_Mouse$results
  ResMice[ResMice$sd_from_mean < 0, "sd_from_mean"] <- 0
  
  miceBroad <- EWCE::ewce_plot(unconditional_results_Mouse$results, mtc_method = "BH")
  
  return(list(unconditional_results_Mouse$results, miceBroad))
}

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
# Function: rename_mod
#
# Description:
#   Re-label module names of the form "M<number>" into a modified form.
#   For M >= 10, subtract 1 and append "A"; for others, simply append "A".
#   (Used to re-index and label modules in the heatmap.)
#
# Arguments:
#   x : character vector of module IDs (e.g. "M1", "M10").
#
# Returns:
#   character vector of renamed modules (e.g. "M1A", "M9A").
############################################
rename_mod <- function(x) {
  nums <- as.integer(str_match(x, "^M(\\d+)$")[,2])  # NA if not M<digits>
  out  <- x
  is_m <- !is.na(nums)
  
  nums_new <- ifelse(is_m & nums >= 10, nums - 1L, nums)
  out[is_m] <- paste0("M", nums_new[is_m], "A")
  out
}

#########################################################################################
####################################### MAIN ############################################
#########################################################################################

##############################
# Load EWCE reference objects
# and gene set inputs
##############################

setwd("C:\\Users\\gp487\\OneDrive - University of Exeter\\AD_infection\\Analyses\\sncEnrichment\\")

load("input\\objects\\ctd_scROS.rda")

targetListWGCNA_yellow <- read.csv("input\\data\\yellow_module.csv", row.names = 1)
targetListWGCNA_tan    <- read.csv("input\\data\\tan_module.csv",    row.names = 1)
targetListWGCNA_red    <- read.csv("input\\data\\red_module.csv",    row.names = 1)

hubgenesYellow <- read.table(
  "C:/Users/gp487/OneDrive - University of Exeter/AD_infection/Analyses/RNA_sequencing/output/data/hub_genes_yellow.txt",
  header = TRUE
)

limma_AD_vs_Control_no  <- read.table("input\\data\\sign_AD_vs_Control_no.txt",  header = TRUE)
limma_AD_vs_Control_inf <- read.table("input\\data\\sign_AD_vs_Control_inf.txt", header = TRUE)

yellow_module <- read.delim(
  "RNA_sequencing\\output\\data\\module_discovery_yellow_module.tsv",
  stringsAsFactors = FALSE
)
yellow_module <- yellow_module[, 1:2]

tan_module <- read.delim(
  "RNA_sequencing\\output\\data\\module_discovery_tan_module.tsv",
  stringsAsFactors = FALSE
)
tan_module <- tan_module[, 1:2]

red_module <- read.delim(
  "RNA_sequencing\\output\\data\\module_discovery_red_module.tsv",
  stringsAsFactors = FALSE
)
red_module <- red_module[, 1:2]

############################################
# EWCE on WGCNA modules (yellow, tan, red)
############################################

targetListWGCNA_yellow <- convertToSymbol(targetListWGCNA_yellow$x)
targetListWGCNA_tan    <- convertToSymbol(targetListWGCNA_tan$x)

results_ann1WGCNA_tan <- calculateAnnotation1(ctd, targetListWGCNA_tan$hgnc_symbol, "WGCNA tan module")
results_ann2WGCNA_tan <- calculateAnnotation2(ctd, targetListWGCNA_pink$hgnc_symbol, "WGCNA tan module")

results_ann2WGCNA_tan_IMMUNE <- calculateAnnotation1(
  scPFC427_Immune_ctd[["ctd"]],
  hubgenesYellow$hgnc_symbol,
  "WGCNA yellow module hub genes on immune cells"
)

targetListWGCNA_red <- convertToSymbol(targetListWGCNA_red$x)

results_ann1WGCNA_red <- calculateAnnotation1(ctd, targetListWGCNA_red$hgnc_symbol, "WGCNA red module")
results_ann2WGCNA_red <- calculateAnnotation2(ctd, targetListWGCNA_red$hgnc_symbol, "WGCNA red module")

tiff("C:\\Users\\gp487\\OneDrive - University of Exeter\\AD_infection\\Analyses\\RNA_sequencing\\output\\figures\\EWCE\\results_ann1WGCNA_tan.tiff",
     height = 400, width = 700)
results_ann1WGCNA_tan
dev.off()

tiff("C:\\Users\\gp487\\OneDrive - University of Exeter\\AD_infection\\Analyses\\RNA_sequencing\\output\\figures\\EWCE\\results_ann1WGCNA_red.tiff",
     height = 400, width = 700)
results_ann1WGCNA_red
dev.off()

############################################
# EWCE on WGCNA hub genes (yellow module)
############################################

hubgenesYellow <- convertToSymbol(hubgenesYellow$x)

results_ann1hubWGCNA_yellow <- calculateAnnotation1(ctd, hubgenesYellow$hgnc_symbol,
                                                    "WGCNA yellow module hub genes")
results_ann2hubWGCNA_yellow <- calculateAnnotation2(ctd, hubgenesYellow$hgnc_symbol,
                                                    "WGCNA yellow module hub genes")
results_ann2hubWGCNA_yellow_IMMUNE <- calculateAnnotation1(
  scPFC427_Immune_ctd[["ctd"]],
  hubgenesYellow$hgnc_symbol,
  "WGCNA yellow module hub genes on immune cells"
)

############################################
# EWCE on limma DE gene sets
############################################

limma_AD_vs_Control_no <- convertToSymbol(limma_AD_vs_Control_no$x)
results_ann1_limma_AD_vs_Control_no <- calculateAnnotation1(
  ctd,
  limma_AD_vs_Control_no$hgnc_symbol,
  "DEGs AD vs Control without infection"
)
results_ann2_limma_AD_vs_Control_no <- calculateAnnotation2(
  ctd,
  limma_AD_vs_Control_no$hgnc_symbol,
  "DEGs AD vs Control without infection"
)

limma_AD_vs_Control_inf <- convertToSymbol(limma_AD_vs_Control_inf$x)
results_ann1_limma_AD_vs_Control_inf <- calculateAnnotation1(
  ctd,
  limma_AD_vs_Control_inf$hgnc_symbol,
  "DEGs AD vs Control with infection"
)
results_ann2_limma_AD_vs_Control_inf <- calculateAnnotation2(
  ctd,
  limma_AD_vs_Control_inf$hgnc_symbol,
  "DEGs AD vs Control with infection"
)

ggarrange(
  results_ann1Comparison1[[2]],
  results_ann1Comparison2[[2]],
  results_ann1Comparison5[[2]],
  results_ann1Comparison6[[2]]
)


############################################
# EWCE on miNA module discovery (yellow)
############################################

m1  <- yellow_module[yellow_module$CLUSTER_NAME == "M1", "CLUSTER_GENES"];  m1  <- list(unlist(str_split(m1,  ",")))
results_ann1_m1  <- calculateAnnotation1(ctd, m1[[1]],  "Module M1")
results_ann2_m1  <- calculateAnnotation2(ctd, m1[[1]],  "Module M1")

m2  <- yellow_module[yellow_module$CLUSTER_NAME == "M2", "CLUSTER_GENES"];  m2  <- list(unlist(str_split(m2,  ",")))
results_ann1_m2  <- calculateAnnotation1(ctd, m2[[1]],  "Module M2")
results_ann2_m2  <- calculateAnnotation2(ctd, m2[[1]],  "Module M2")

m3  <- yellow_module[yellow_module$CLUSTER_NAME == "M3", "CLUSTER_GENES"];  m3  <- list(unlist(str_split(m3,  ",")))
results_ann1_m3  <- calculateAnnotation1(ctd, m3[[1]],  "Module M3")
results_ann2_m3  <- calculateAnnotation2(ctd, m3[[1]],  "Module M3")

m4  <- yellow_module[yellow_module$CLUSTER_NAME == "M4", "CLUSTER_GENES"];  m4  <- list(unlist(str_split(m4,  ",")))
results_ann1_m4  <- calculateAnnotation1(ctd, m4[[1]],  "Module M4")
results_ann2_m4  <- calculateAnnotation2(ctd, m4[[1]],  "Module M4")

m5  <- yellow_module[yellow_module$CLUSTER_NAME == "M5", "CLUSTER_GENES"];  m5  <- list(unlist(str_split(m5,  ",")))
results_ann1_m5  <- calculateAnnotation1(ctd, m5[[1]],  "Module M5")
results_ann2_m5  <- calculateAnnotation2(ctd, m5[[1]],  "Module M5")

m6  <- yellow_module[yellow_module$CLUSTER_NAME == "M6", "CLUSTER_GENES"];  m6  <- list(unlist(str_split(m6,  ",")))
results_ann1_m6  <- calculateAnnotation1(ctd, m6[[1]],  "Module M6")
results_ann2_m6  <- calculateAnnotation2(ctd, m6[[1]],  "Module M6")

m7  <- yellow_module[yellow_module$CLUSTER_NAME == "M7", "CLUSTER_GENES"];  m7  <- list(unlist(str_split(m7,  ",")))
results_ann1_m7  <- calculateAnnotation1(ctd, m7[[1]],  "Module M7")
results_ann2_m7  <- calculateAnnotation2(ctd, m7[[1]],  "Module M7")

m8  <- yellow_module[yellow_module$CLUSTER_NAME == "M8", "CLUSTER_GENES"];  m8  <- list(unlist(str_split(m8,  ",")))
results_ann1_m8  <- calculateAnnotation1(ctd, m8[[1]],  "Module M8")
results_ann2_m8  <- calculateAnnotation2(ctd, m8[[1]],  "Module M8")

m10 <- yellow_module[yellow_module$CLUSTER_NAME == "M10","CLUSTER_GENES"]; m10 <- list(unlist(str_split(m10, ",")))
results_ann1_m10 <- calculateAnnotation1(ctd, m10[[1]], "Module M10")
results_ann2_m10 <- calculateAnnotation2(ctd, m10[[1]], "Module M10")

m11 <- yellow_module[yellow_module$CLUSTER_NAME == "M11","CLUSTER_GENES"]; m11 <- list(unlist(str_split(m11, ",")))
results_ann1_m11 <- calculateAnnotation1(ctd, m11[[1]], "Module M11")
results_ann2_m11 <- calculateAnnotation2(ctd, m11[[1]], "Module M11")

############################################
# EWCE on tan module subclusters
############################################

m1_tan <- tan_module[tan_module$CLUSTER_NAME == "M1", "CLUSTER_GENES"]; m1_tan <- list(unlist(str_split(m1_tan, ",")))
results_ann1_m1_tan <- calculateAnnotation1(ctd, m1_tan[[1]], "Module M1")
results_ann2_m1_tan <- calculateAnnotation2(ctd, m1_tan[[1]], "Module M1")

m2_tan <- tan_module[tan_module$CLUSTER_NAME == "M5", "CLUSTER_GENES"]; m2_tan <- list(unlist(str_split(m2_tan, ",")))
results_ann1_m2_tan <- calculateAnnotation1(ctd, m2_tan[[1]], "Module M2")
results_ann2_m2_tan <- calculateAnnotation2(ctd, m2_tan[[1]], "Module M2")

m3_tan <- tan_module[tan_module$CLUSTER_NAME == "M6", "CLUSTER_GENES"]; m3_tan <- list(unlist(str_split(m3_tan, ",")))
results_ann1_m3_tan <- calculateAnnotation1(ctd, m3_tan[[1]], "Module M3")
results_ann2_m3_tan <- calculateAnnotation2(ctd, m3_tan[[1]], "Module M3")

############################################
# EWCE on red module subclusters
############################################

m1_red <- red_module[red_module$CLUSTER_NAME == "M1", "CLUSTER_GENES"]; m1_red <- list(unlist(str_split(m1_red, ",")))
results_ann1_m1_red <- calculateAnnotation1(ctd, m1_red[[1]], "Module M1")
results_ann2_m1_red <- calculateAnnotation2(ctd, m1_red[[1]], "Module M1")

m2_red <- red_module[red_module$CLUSTER_NAME == "M2", "CLUSTER_GENES"]; m2_red <- list(unlist(str_split(m2_red, ",")))
results_ann1_m2_red <- calculateAnnotation1(ctd, m2_red[[1]], "Module M2")
results_ann2_m2_red <- calculateAnnotation2(ctd, m2_red[[1]], "Module M2")

m3_red <- red_module[red_module$CLUSTER_NAME == "M3", "CLUSTER_GENES"]; m3_red <- list(unlist(str_split(m3_red, ",")))
results_ann1_m3_red <- calculateAnnotation1(ctd, m3_red[[1]], "Module M3")
results_ann2_m3_red <- calculateAnnotation2(ctd, m3_red[[1]], "Module M3")

m4_red <- red_module[red_module$CLUSTER_NAME == "M4", "CLUSTER_GENES"]; m4_red <- list(unlist(str_split(m4_red, ",")))
results_ann1_m4_red <- calculateAnnotation1(ctd, m4_red[[1]], "Module M4")
results_ann2_m4_red <- calculateAnnotation2(ctd, m4_red[[1]], "Module M4")

############################################
# Aggregate annotation-level 1 results
# across modules and build heatmap
############################################

# Optional mapping to pretty cell-type names
cell_map <- tibble(
  CellType         = c("Oli","Opc","Per","End","Mic","Ast","In","Ex"),
  CellType_pretty  = c("Oligodendrocytes","OPC","Pericytes","Endothelia",
                       "Microglia","Astrocytes","Inhibitory Neurons","Excitatory Neurons")
)

# Collect results_ann1_m* objects from environment
obj_names <- ls(pattern = "^results_ann1_m[0-9]+$")

ann1_all <- map_dfr(obj_names, function(nm) {
  x   <- get(nm, inherits = TRUE)
  tab <- tryCatch(x[[1]]$results, error = function(e) NULL)
  if (is.null(tab)) return(tibble())
  
  if (!"CellType" %in% names(tab)) {
    tab <- tibble(CellType = rownames(tab)) %>% bind_cols(as_tibble(tab))
  }
  
  module_lab <- nm %>%
    sub("^results_ann1_?", "", .) %>%
    str_replace_all("_", " ") %>%
    str_replace_all("\\.", " ") %>%
    toupper()
  
  as_tibble(tab) %>%
    mutate(
      Module    = module_lab,
      neglog10q = -log10(pmax(q, .Machine$double.xmin))
    ) %>%
    select(Module, CellType, q, neglog10q, sd_from_mean, everything()) %>%
    left_join(cell_map, by = "CellType") %>%
    mutate(CellType = coalesce(CellType_pretty, CellType)) %>%
    select(-CellType_pretty)
})

ann1_all <- ann1_all %>%
  mutate(
    Module = rename_mod(as.character(Module)),
    .ord   = as.integer(str_extract(Module, "\\d+"))
  ) %>%
  arrange(.ord) %>%
  mutate(Module = factor(Module, levels = unique(Module))) %>%
  select(-.ord)

L <- ann1_all$sd_from_mean %>% abs() %>% quantile(0.99, na.rm = TRUE) %>% as.numeric()

ann1_hm <- ann1_all %>%
  mutate(
    sd_cap = pmax(pmin(sd_from_mean, L), -L),
    sig    = case_when(
      q < 0.001 ~ "***",
      q < 0.01  ~ "**",
      q < 0.05  ~ "*",
      TRUE      ~ ""
    )
  )

p_heat <- ggplot(ann1_hm, aes(x = Module, y = CellType, fill = sd_cap)) +
  geom_tile(color = "grey92", linewidth = 0.25) +
  scale_fill_gradient2(
    name     = "sd from mean",
    low      = "#1f4e79",
    mid      = "white",
    high     = "#a91515",
    midpoint = 0,
    na.value = "grey90",
    limits   = c(-L, L)
  ) +
  geom_text(
    data  = subset(ann1_hm, q < 0.05),
    aes(label = sig),
    size     = 6,
    fontface = "bold",
    color    = "black"
  ) +
  labs(
    title    = "Cell-type enrichment",
    subtitle = "",
    x        = NULL,
    y        = NULL
  ) +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(face = "bold")
  )

p_heat +
  scale_x_discrete(labels = function(x) {
    nums   <- as.integer(stringr::str_match(x, "^M(\\d+)$")[,2])
    is_m   <- !is.na(nums)
    nums_new <- ifelse(is_m & nums >= 10, nums - 1L, nums)
    out    <- x
    out[is_m] <- paste0("M", nums_new[is_m], "A")
    out
  })

print(p_heat)

############################################
# Heatmap for tan submodules
############################################

obj_names_tan <- ls(pattern = "^results_ann1_m[0-9]+_tan$")

ann1_all_tan <- map_dfr(obj_names_tan, function(nm) {
  x   <- get(nm, inherits = TRUE)
  tab <- tryCatch(x[[1]]$results, error = function(e) NULL)
  if (is.null(tab)) return(tibble())
  
  tab       <- as_tibble(tab)
  names(tab) <- gsub("\\.", "_", names(tab))
  
  if (!"CellType" %in% names(tab)) {
    tab <- tibble(CellType = rownames(x[[1]]$results)) %>% bind_cols(tab)
  }
  
  module_lab <- nm %>%
    sub("^results_ann1_?", "", .) %>%
    str_replace_all("_|\\.", " ") %>%
    toupper() %>%
    str_squish() %>%
    sub(" TAN$", "B", .)
  
  tab %>%
    mutate(
      Module    = module_lab,
      neglog10q = -log10(pmax(q, .Machine$double.xmin))
    ) %>%
    dplyr::select(Module, CellType, q, neglog10q, sd_from_mean, dplyr::everything())
})

ann1_all_tan <- ann1_all_tan %>%
  mutate(CellType_key = tolower(str_squish(CellType))) %>%
  left_join(cell_map, by = "CellType") %>%
  mutate(CellType = coalesce(CellType_pretty, CellType)) %>%
  dplyr::select(-CellType_pretty)

L <- ann1_all_tan$sd_from_mean %>%
  abs() %>% quantile(0.99, na.rm = TRUE) %>% as.numeric()

ann1_hm_tan <- ann1_all_tan %>%
  mutate(
    sd_cap = pmax(pmin(sd_from_mean, L), -L),
    sig    = case_when(
      q < 0.001 ~ "***",
      q < 0.01  ~ "**",
      q < 0.05  ~ "*",
      TRUE      ~ ""
    )
  )

p_heat_tan <- ggplot(ann1_hm_tan, aes(x = Module, y = CellType, fill = sd_cap)) +
  geom_tile(color = "grey92", linewidth = 0.25) +
  scale_fill_gradient2(
    name     = "sd from mean",
    low      = "#1f4e79",
    mid      = "white",
    high     = "#a91515",
    midpoint = 0,
    na.value = "grey90",
    limits   = c(-L, L)
  ) +
  geom_text(
    data  = subset(ann1_hm_tan, q < 0.05),
    aes(label = sig),
    size     = 6,
    fontface = "bold",
    color    = "black"
  ) +
  labs(title = "Cell-type enrichment", x = NULL, y = NULL) +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(face = "bold")
  )

print(p_heat_tan)

############################################
# Heatmap for red submodules
############################################

obj_names_red <- ls(pattern = "^results_ann1_m[0-9]+_red$")

ann1_all_red <- map_dfr(obj_names_red, function(nm) {
  x   <- get(nm, inherits = TRUE)
  tab <- tryCatch(x[[1]]$results, error = function(e) NULL)
  if (is.null(tab)) return(tibble())
  
  tab        <- as_tibble(tab)
  names(tab) <- gsub("\\.", "_", names(tab))
  
  if (!"CellType" %in% names(tab)) {
    tab <- tibble(CellType = rownames(x[[1]]$results)) %>% bind_cols(tab)
  }
  
  module_lab <- nm %>%
    sub("^results_ann1_?", "", .) %>%
    str_replace_all("_|\\.", " ") %>%
    toupper() %>%
    str_squish() %>%
    sub(" RED$", "C", .)
  
  tab %>%
    mutate(
      Module    = module_lab,
      neglog10q = -log10(pmax(q, .Machine$double.xmin))
    ) %>%
    dplyr::select(Module, CellType, q, neglog10q, sd_from_mean, dplyr::everything())
})

ann1_all_red <- ann1_all_red %>%
  mutate(CellType_key = tolower(str_squish(CellType))) %>%
  left_join(cell_map, by = "CellType") %>%
  mutate(CellType = coalesce(CellType_pretty, CellType)) %>%
  dplyr::select(-CellType_pretty)

L <- ann1_all_red$sd_from_mean %>%
  abs() %>% quantile(0.99, na.rm = TRUE) %>% as.numeric()

ann1_hm_red <- ann1_all_red %>%
  mutate(
    sd_cap = pmax(pmin(sd_from_mean, L), -L),
    sig    = case_when(
      q < 0.001 ~ "***",
      q < 0.01  ~ "**",
      q < 0.05  ~ "*",
      TRUE      ~ ""
    )
  )

p_heat_red <- ggplot(ann1_hm_red, aes(x = Module, y = CellType, fill = sd_cap)) +
  geom_tile(color = "grey92", linewidth = 0.25) +
  scale_fill_gradient2(
    name     = "sd from mean",
    low      = "#1f4e79",
    mid      = "white",
    high     = "#a91515",
    midpoint = 0,
    na.value = "grey90",
    limits   = c(-L, L)
  ) +
  geom_text(
    data  = subset(ann1_hm_red, q < 0.05),
    aes(label = sig),
    size     = 6,
    fontface = "bold",
    color    = "black"
  ) +
  labs(title = "Cell-type enrichment", x = NULL, y = NULL) +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(face = "bold")
  )

print(p_heat_red)

tiff("RNA_sequencing\\output\\figures\\EWCE\\p_heat_red.tiff",
     height = 700, width = 700)
p_heat_red
dev.off()
