########################
# @author: Giulia Pegoraro
#
# GSVA + limma pathway-level analysis using Gene Ontology (GO)
#
# This script:
#  - Loads log2-expression data and matched phenotype for the respiratory infection cohort
#  - Converts Ensembl gene IDs (with version) to HGNC symbols
#  - Builds GO-based pathway gene sets (Biological Process, from org.Hs.eg.db + GO.db)
#  - Runs GSVA to obtain pathway activity scores (ES matrix)
#  - Fits limma models on GSVA scores with AD, infection, and covariates
#  - Tests pathway-level contrasts (AD, infection, and interaction)
#  - Exports per-contrast results and a combined table across all contrasts
#  - Optionally visualises significant pathways and performs pathway clustering
########################

##########################################################
######################## LIBRARIES #######################
##########################################################
suppressPackageStartupMessages({
  library(GSVA)
  library(limma)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(GO.db)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(PanomiR)
  library(ensembldb)
})

##########################################################
######################## FUNCTIONS #######################
##########################################################

############################################
# Function: convertFromEnsembl
#
# Description:
#   Convert Ensembl gene IDs with version (e.g. ENSG00000141510.12)
#   to another identifier (e.g. HGNC symbol) using a fixed Ensembl
#   archive for reproducibility.
#
# Arguments:
#   targetList          : character vector of Ensembl IDs with version.
#   conversionAttribute : character, Ensembl attribute to retrieve
#                         (default: "hgnc_symbol").
#
# Returns:
#   data frame with columns:
#     ensembl_gene_id_version
#     <conversionAttribute>   (e.g. hgnc_symbol)
#   Rows with empty conversionAttribute are removed.
############################################
convertFromEnsembl <- function(targetList, conversionAttribute = "hgnc_symbol") {
  # Use the Ensembl 104 archive for reproducibility
  ensembl <- biomaRt::useEnsembl(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl",
    version = 104
  )
  
  # Query Ensembl: input Ensembl IDs WITH version (e.g., ENSG00000141510.12)
  result <- biomaRt::getBM(
    attributes = c("ensembl_gene_id_version", conversionAttribute),
    filters    = "ensembl_gene_id_version",
    values     = targetList,
    mart       = ensembl
  )
  
  # Optionally remove empty results
  result <- result[result[[conversionAttribute]] != "", , drop = FALSE]
  
  return(result)
}

############################################
# Function: write_tt_go
#
# Description:
#   Extract a limma topTable for a given GSVA contrast, reorder columns,
#   and write to CSV.
#
# Arguments:
#   fit2       : MArrayLM object after contrasts.fit + eBayes.
#   coef_name  : character, name of coefficient/contrast.
#   out_prefix : character, prefix for output filename.
#
# Returns:
#   data frame with Pathway as first column and limma statistics.
############################################
write_tt_go <- function(fit2, coef_name, out_prefix) {
  tt <- limma::topTable(
    fit2,
    coef          = coef_name,
    number        = Inf,
    sort.by       = "P",
    adjust.method = "fdr"
  )
  tt$Pathway <- rownames(tt)
  tt <- dplyr::relocate(tt, Pathway)
  
  out_csv <- file.path(
    out_dir_data,
    sprintf("%s_%s_GO_%s_%s.csv", out_prefix, coef_name, go_onto, go_subset_tag)
  )
  readr::write_csv(tt, out_csv)
  tt
}

############################################
# Function: combine_results
#
# Description:
#   Merge multiple per-contrast GSVA limma results tables into a single
#   wide table, suffixing statistic columns by contrast name.
#
# Arguments:
#   ... : data frames as returned by write_tt_go().
#
# Returns:
#   data frame with one row per pathway and all contrast-specific
#   statistics in separate columns.
############################################
combine_results <- function(...) {
  dfs <- list(...)
  nm  <- c(
    "AD_with_vs_without_infection",
    "Control_with_vs_without_infection",
    "AD_vs_Control_without_infection",
    "AD_vs_Control_with_infection"
  )
  names(dfs) <- nm[seq_along(dfs)]
  
  dfs <- lapply(dfs, function(d) {
    d %>% dplyr::select(Pathway, logFC, AveExpr, t, P.Value, adj.P.Val, B)
  })
  
  for (i in seq_along(dfs)) {
    names(dfs[[i]])[-1] <- paste0(names(dfs[[i]])[-1], "__", names(dfs)[i])
  }
  Reduce(function(x, y) dplyr::full_join(x, y, by = "Pathway"), dfs)
}

############################################
# Function: plot_sig_pathways
#
# Description:
#   Create a horizontal barplot of logFC for significantly enriched
#   pathways (FDR < 0.05), coloured by direction (up/down), and save
#   as PNG.
#
# Arguments:
#   tt             : data frame, limma topTable for a given contrast.
#   contrast_label : character, label used in plot title / filename.
#
# Returns:
#   None (PNG file is written to disk; returns invisibly).
############################################
plot_sig_pathways <- function(tt, contrast_label) {
  sig <- tt %>%
    dplyr::filter(!is.na(adj.P.Val), adj.P.Val < 0.05)
  
  if (nrow(sig) == 0) {
    message(sprintf("[PLOT] No significant pathways (FDR < 0.05) for %s", contrast_label))
    return(invisible(NULL))
  }
  
  sig <- sig %>% dplyr::arrange(logFC, adj.P.Val)
  cols <- ifelse(sig$logFC < 0, "blue", "red")
  safe_name <- gsub("[^A-Za-z0-9_]+", "_", contrast_label)
  
  png_h <- max(1200, 25 * nrow(sig))
  png(
    file.path(out_dir_figs, sprintf("GSVA_GO_%s_sig_%s_%s.png", go_onto, safe_name, go_subset_tag)),
    width  = 2500,
    height = png_h
  )
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mar = c(6, 75, 6, 2) + 0.1)
  
  barplot(
    sig$logFC,
    names.arg = sig$Pathway,
    horiz     = TRUE,
    las       = 1,
    col       = cols,
    xlab      = "logFC (GSVA score)",
    main      = sprintf("Significant pathways (Padj < 0.05): %s — GO %s",
                        contrast_label, go_onto),
    cex.names = 2,
    cex.axis  = 2,
    cex.lab   = 2,
    cex.main  = 2.5
  )
  abline(v = 0, lty = 2)
  legend(
    "bottomright", inset = 0.01, bty = "n",
    fill   = c("red", "blue"),
    legend = c("Up (logFC > 0)", "Down (logFC < 0)"),
    cex    = 1.3
  )
  dev.off()
}

##########################################################
######################## SETTINGS ########################
##########################################################

# Base project directory (edit if needed)
base_dir <- "C:/Users/gp487/OneDrive - University of Exeter/AD_infection/Analyses"
setwd(base_dir)

# Input files produced by upstream RNA-seq preprocessing
pheno_file <- file.path(base_dir, "RNA_sequencing/input/data/Metadata/pheno_resp_infection_no_out.csv")
expr_file  <- file.path(base_dir, "RNA_sequencing/input/data/Expression_data/gene_expr_resp_inf_no_out.csv")

# Output locations
out_dir_data <- file.path(base_dir, "RNA_sequencing/output/data/GSVA_GO")
out_dir_figs <- file.path(base_dir, "RNA_sequencing/output/figures/GSVA_GO")
if (!dir.exists(out_dir_data)) dir.create(out_dir_data, recursive = TRUE)
if (!dir.exists(out_dir_figs)) dir.create(out_dir_figs, recursive = TRUE)

# GO ontology choice
go_onto <- "BP"   # "BP" (Biological Process), "MF", or "CC"

# Tag for filenames (currently using all BP terms)
go_subset_tag <- "ALL_BP"

# GSVA options
min_gs_size <- 10
max_gs_size <- 500
method_gsva <- "gsva"
kcdf        <- "Gaussian"

##########################################################
######################## LOAD DATA #######################
##########################################################

# Expression matrix (genes x samples); Ensembl IDs with version as rownames
expr <- as.matrix(read.csv(expr_file, row.names = 1, check.names = FALSE))

# Phenotype
pheno <- read.csv(pheno_file, row.names = 1, check.names = FALSE)
if (!"BBNId" %in% colnames(pheno)) pheno$BBNId <- rownames(pheno)

# Ensure sample order matches between expression and phenotype
if (!identical(colnames(expr), pheno$BBNId)) {
  pheno <- pheno[match(colnames(expr), pheno$BBNId), ]
}
stopifnot(identical(colnames(expr), pheno$BBNId))

##########################################################
########### MAP ENSEMBL IDs (with version) → SYMBOL ######
##########################################################

# Extract Ensembl IDs (with version) and convert to HGNC symbols
ens_with_version <- rownames(expr)
conv <- convertFromEnsembl(ens_with_version, conversionAttribute = "hgnc_symbol")

# Lookup: Ensembl-with-version → SYMBOL
sym_lookup <- setNames(conv$hgnc_symbol, conv$ensembl_gene_id_version)
symvec     <- sym_lookup[rownames(expr)]

# Keep only successfully mapped genes
keep <- !is.na(symvec) & nzchar(symvec)
expr_sym <- expr[keep, , drop = FALSE]
rownames(expr_sym) <- symvec[keep]

# Collapse duplicates (multiple Ensembl IDs mapping to the same gene symbol)
expr_sym <- limma::avereps(expr_sym)

##########################################################
################# BUILD GO GENE SETS #####################
##########################################################
message("Building GO gene sets (ALL BP) from org.Hs.eg.db / GO.db ...")

# Use SYMBOLs present in the expression matrix
common_genes <- rownames(expr_sym)

# Map SYMBOL → GO; keep only selected ontology
map_go <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys     = common_genes,
  keytype  = "SYMBOL",
  columns  = c("GO", "ONTOLOGY")
)
map_go <- map_go[!is.na(map_go$GO) & map_go$ONTOLOGY == go_onto, c("GO", "SYMBOL")]

# Fetch GO term names for readability
term_names <- AnnotationDbi::select(
  GO.db,
  keys    = unique(map_go$GO),
  keytype = "GOID",
  columns = c("TERM")
)
colnames(term_names) <- c("GO", "TERM")

# Attach term names and build display names "TERM [GO:ID]"
map_go <- dplyr::left_join(map_go, term_names, by = "GO")
map_go$TERM <- ifelse(is.na(map_go$TERM), map_go$GO, map_go$TERM)
map_go$set  <- paste0(map_go$TERM, " [", map_go$GO, "]")

# Save the exact GO terms used
terms_used <- map_go %>%
  distinct(GO, TERM) %>%
  arrange(TERM)
readr::write_csv(
  terms_used,
  file.path(out_dir_data, sprintf("GO_terms_used_%s_%s.csv", go_onto, go_subset_tag))
)

# Construct list of gene sets and apply size filter
genesets <- split(map_go$SYMBOL, map_go$set)
# Deduplicate genes per set and intersect with available genes
genesets <- lapply(genesets, function(x) intersect(unique(x), common_genes))
# Apply size bounds
genesets <- genesets[lengths(genesets) >= min_gs_size & lengths(genesets) <= max_gs_size]

message(sprintf("Retained %d GO %s gene sets (size %d..%d).",
                length(genesets), go_onto, min_gs_size, max_gs_size))

##########################################################
######################## RUN GSVA ########################
##########################################################
message("Running GSVA on GO BP gene sets ...")

param <- GSVA::gsvaParam(
  expr_sym,
  genesets,
  minSize = min_gs_size,
  maxSize = max_gs_size,
  kcdf    = kcdf
)

# Pathway x sample ES matrix
es <- GSVA::gsva(param, verbose = FALSE)

# Save ES matrix
write.csv(
  es,
  file.path(out_dir_data, sprintf("GSVA_ES_GO_%s_%s.csv", go_onto, go_subset_tag))
)

##########################################################
################ DESIGN & CONTRASTS ######################
##########################################################

# Encode phenotype covariates as factors / numeric
factor_list <- c(
  "AD","Infection","Gender",
  "Brain.Bank_Bristol","Brain.Bank_Edinburgh",
  "Brain.Bank_Manchester","Brain.Bank_Queen_Square"
)
for (i in factor_list) pheno[[i]] <- as.factor(pheno[[i]])

num_list <- c(
  "Age","fpkm_astrocytes","fpkm_endothelial",
  "fpkm_microglia","fpkm_neurons","fpkm_OPC",
  "fpkm_oligodendrocytes"
)
for (i in num_list) pheno[[i]] <- suppressWarnings(as.numeric(pheno[[i]]))

# Design matrix (same structure as gene-level limma)
design <- model.matrix(
  ~ AD * Infection +
    Age + Gender +
    Brain.Bank_Bristol + 
    Brain.Bank_Edinburgh + 
    Brain.Bank_Manchester + 
    Brain.Bank_Queen_Square +
    Brain.Bank_Newcastle +
    fpkm_astrocytes + 
    fpkm_endothelial + 
    fpkm_microglia + 
    fpkm_neurons + 
    fpkm_OPC + 
    fpkm_oligodendrocytes,
  data = pheno
)
colnames(design) <- make.names(colnames(design))

# Align design with GSVA ES columns
rownames(design) <- pheno$BBNId
stopifnot(identical(colnames(es), rownames(design)))

# Fit limma model at pathway level
fit <- limma::lmFit(es, design)
fit <- limma::eBayes(fit)

# Contrasts: same as gene-level analysis
contrast.matrix <- makeContrasts(
  AD_with_vs_without_infection      = Infection1 + AD1.Infection1,
  Control_with_vs_without_infection = Infection1,
  AD_vs_Control_without_infection   = AD1,
  AD_vs_Control_with_infection      = AD1 + AD1.Infection1,
  levels = design
)

fit2 <- limma::contrasts.fit(fit, contrast.matrix)
fit2 <- limma::eBayes(fit2)

##########################################################
######################## RESULTS #########################
##########################################################

# Per-contrast tables
res_AD_with_vs_without_infection      <- write_tt_go(fit2, "AD_with_vs_without_infection",      "GSVA_limma")
res_Control_with_vs_without_infection <- write_tt_go(fit2, "Control_with_vs_without_infection", "GSVA_limma")
res_AD_vs_Control_no                  <- write_tt_go(fit2, "AD_vs_Control_without_infection",   "GSVA_limma")
res_AD_vs_Control_inf                 <- write_tt_go(fit2, "AD_vs_Control_with_infection",      "GSVA_limma")

# Combined table across all contrasts
res_combined <- combine_results(
  res_AD_with_vs_without_infection,
  res_Control_with_vs_without_infection,
  res_AD_vs_Control_no,
  res_AD_vs_Control_inf
)
readr::write_csv(
  res_combined,
  file.path(out_dir_data, sprintf("GSVA_limma_ALL_CONTRASTS_GO_%s_%s.csv", go_onto, go_subset_tag))
)

# Optional visualisation of significant pathways
plot_sig_pathways(res_AD_with_vs_without_infection,      "AD: with vs without infection")
plot_sig_pathways(res_Control_with_vs_without_infection, "Control: with vs without infection")
plot_sig_pathways(res_AD_vs_Control_no,                  "AD vs Control: no infection")
plot_sig_pathways(res_AD_vs_Control_inf,                 "AD vs Control: with infection")


