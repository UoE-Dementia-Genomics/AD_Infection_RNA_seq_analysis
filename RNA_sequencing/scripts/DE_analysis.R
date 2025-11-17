########################
# @author: Giulia Pegoraro
#
# RNA-sequencing data analysis
# Limma RNA-sequencing pipeline
#
# This script:
#  - Loads normalized RNA-seq expression data and phenotype metadata
#  - Prepares covariates (factors, numeric variables, infection types)
#  - Selects a respiratory infection cohort and removes outliers
#  - Fits a limma linear model with AD, infection, and covariates
#  - Computes differential expression contrasts of interest
#  - Maps Ensembl IDs to HGNC symbols and exports DE results
#  - Generates Venn diagrams and volcano plots for key comparisons
########################

##########################################################
######################## LIBRARIES #######################
##########################################################

# Load required libraries for data processing, statistics, and visualization
library(readxl)
library(limma)
library(edgeR)
library(WGCNA)
library(biomaRt)
library(org.Hs.eg.db)
library(magrittr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)
library(corrplot)
library(ggvenn)
library(VennDiagram)
library(ggpubr)
library(EnhancedVolcano)
library(fastDummies)
library(dplyr)
library(stringr)

##########################################################
######################## FUNCTIONS #######################
##########################################################

############################################
# Function: convertToSymbol
#
# Description:
#   Map Ensembl gene IDs including version (e.g. ENSG000001.14) to HGNC
#   gene symbols using Ensembl via biomaRt.
#
# Arguments:
#   targetList : character vector of Ensembl IDs with version.
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
# Function: pca.plot
#
# Description:
#   Generate a PCA scatter plot of the first two PCs, colouring samples
#   by a selected variable, and label points. The proportion of variance
#   explained by PC1 and PC2 is reported in the axis labels.
#
# Arguments:
#   pca_res     : prcomp object containing PCA results.
#   sample_info : data frame with sample metadata, including BBNId.
#   title       : character, plot title and output file prefix.
#   variable    : column name in 'sample_info' used to colour/label points.
#
# Returns:
#   None (PNG file is written to disk).
############################################
pca.plot <- function(pca_res, sample_info, title, variable) {
  pc_sum <- summary(pca_res)$importance
  PC1_var <- signif(as.numeric(pc_sum["Proportion of Variance","PC1"])*100,2)
  PC2_var <- signif(as.numeric(pc_sum["Proportion of Variance","PC2"])*100,2)
  
  sample_data <- pca_res$x %>%
    as_tibble(rownames = "BBNId") %>%
    full_join(sample_info, by = "BBNId")
  
  ggplot(sample_data, aes(x = PC1, y = PC2, colour = variable)) +
    geom_point() +
    geom_text(aes(label = {{ variable }}), nudge_x = 0.02, nudge_y = 0.02, check_overlap = FALSE) +
    theme(text = element_text(size = 0.02)) +
    xlab(paste0("PC1 (", PC1_var, "%)")) +
    ylab(paste0("PC2 (", PC2_var, "%)")) +
    labs(colour = title)
  
  ggsave(file = paste0(title, "_pca_plot", ".png"))
}

############################################
# Function: pca.correlation
#
# Description:
#   Correlate phenotype variables with a given number of principal
#   components (PC1..PCn) and visualise the correlations (with p-values)
#   using a correlation heatmap.
#
# Arguments:
#   pca_score : numeric matrix/data frame of PCA scores (samples x PCs).
#   pheno_df  : data frame with phenotype variables (samples in rows).
#   n_PC      : integer, number of PCs to include (from PC1 to PC_n_PC).
#   file_name : character, output PNG file path.
#
# Returns:
#   None (PNG file is written to disk).
############################################
pca.correlation <- function(pca_score, pheno_df, n_PC, file_name){
  princcomps <- as.data.frame(pca_score[,1:n_PC])
  princcomps_PCA <- merge(pheno_df, princcomps,  by ="row.names")
  princcomps_PCA <- princcomps_PCA[,-1]
  
  # p-values for correlations
  dat_mat <- data.matrix(princcomps_PCA)
  cor_p <- corrplot::cor.mtest(princcomps_PCA, conf.level=0.95)
  p.val_df <- cor_p$p
  p.val_df <- p.val_df[1:ncol(pheno_df),(ncol(pheno_df)+1):ncol(p.val_df)]
  
  # Raw correlation matrix (PCs vs phenotype variables)
  cor_PCA <- stats::cor(
    princcomps_PCA[,(ncol(pheno_df)+1):ncol(princcomps_PCA)], 
    princcomps_PCA[,1:ncol(pheno_df)],
    use = "complete.obs"
  )
  colnames(cor_PCA) <- colnames(pheno_df)
  
  png(file_name, width=2400, height=1500, res=300)
  corrplot::corrplot(
    cor_PCA,
    p.mat = t(p.val_df),
    addrect = 2,
    method = "color",
    insig = "label_sig",
    sig.level = c(0.001, 0.01, 0.05), 
    tl.cex = 0.6,
    cl.cex = 0.4,
    cl.offset = 0.2,
    cl.align.text = "l",
    pch.cex = 0.6,
    pch.col = "grey20"
  )
  dev.off()
}

############################################
# Function: removeVersion
#
# Description:
#   Strip version numbers from Ensembl IDs, e.g.
#   "ENSG00000000003.14" -> "ENSG00000000003".
#
# Arguments:
#   genes : character vector of Ensembl IDs with version.
#
# Returns:
#   list of Ensembl IDs without version (one element per input).
############################################
removeVersion <- function(genes){
  gene_list <- strsplit((genes),"\\.")
  gene_list <- lapply(1:length(gene_list), function(x){ gene_list[[x]][[1]]})
}

##########################################################
########################## MAIN ##########################
##########################################################

############################ Working directory ############################
setwd("C:\\Users\\gp487\\OneDrive - University of Exeter\\AD_infection\\Analyses")

############################ Input expression and phenotype data ############################

# Normalized count matrix (log2 CPM, batch effects removed)
count <- read.csv(
  "input\\data\\Expression_data\\log2.cpm_whitout_be.csv",
  row.names = 1
)

# Phenotype data with batch covariates and infection information
pheno_with_batch <- read.csv("input\\data\\Metadata\\pheno_with_batch_data.csv", row.names = 1)
pheno_infection  <- read.csv("input\\data\\Metadata\\Pheno_Final_with_infection_type.csv")
pheno            <- read.csv("input\\data\\Metadata\\Filtered_pheno_RNA.csv", row.names = 1)

############################ Convert phenotype columns to factors / numeric ############################

factor_list <- c(
  "AD","Infection","Gender",
  "Batch"
)
numeric_list <- c(
  "Age","PMD","TIN.median.",
  "fpkm_astrocytes","fpkm_endothelial","fpkm_microglia","fpkm_neurons",
  "fpkm_OPC","fpkm_oligodendrocytes","fpkm_fetal_quiescent",
  "fpkm_fetal_replicating"
)

for(i in factor_list){
  print(i)
  pheno[,i] <- as.factor(pheno[,i])
}

for(i in numeric_list){
  print(i)
  pheno[,i] <- as.numeric(pheno[,i])
}

# Create dummy variables for Brain.Bank and Batch
pheno <- dummy_cols(
  pheno,
  select_columns = c("Brain.Bank","Batch")
)

pheno[, grep("Brain.Bank_",colnames(pheno), ignore.case = TRUE)] <-
  apply(pheno[, grep("Brain.Bank_",colnames(pheno), ignore.case = TRUE)], 2, as.factor)
pheno[, grep("Batch_",colnames(pheno), ignore.case = TRUE)] <-
  apply(pheno[, grep("Batch_",colnames(pheno), ignore.case = TRUE)], 2, as.factor)

############################ Impute missing values ############################

# Impute missing PMD with median
pheno[is.na(pheno$PMD),"PMD"] <- median(pheno$PMD, na.rm = TRUE)

# Check and align IDs between expression and phenotype matrices
if(!identical(colnames(count), pheno$BBNId)){
  errorCondition("ID not identical")
  pheno <- pheno_data[match(rownames(count),rownames(pheno)),]
}

############################ Clean and merge infection metadata ############################

# Harmonise BBNId in infection metadata
pheno_infection$BBNId <- gsub("\\.", "", pheno_infection$BBNId)
pheno_infection$BBNId <- gsub("_", "", pheno_infection$BBNId)
pheno_infection$BBNId <- gsub("-", "", pheno_infection$BBNId)

# Merge phenotype with infection information
pheno_data_infection <- merge(
  pheno,
  pheno_infection[,c(3,13:23)],
  by.x = "BBNId",
  by.y = "BBNId"
)

if(!identical(colnames(count), pheno_data_infection$BBNId)){
  errorCondition("Not matching Ids")
  pheno_data_infection <- pheno_data_infection[match(colnames(count),pheno_data_infection$BBNId),]
}
write.csv(
  pheno_data_infection,
  "RNA_sequencing\\input\\data\\Expression_data\\pheno_data_infection.csv"
)

############################ Restrict to respiratory infection cohort ############################

selected_col <- c(
  "Cardiac_Infection", "Sepsis", "Urinary_infection",
  "Bedsores", "Parotide_abscess", "Diverculitis", "Immunised",
  "Pyelonephritis", "Influenza", "Peritonitis", "not_specified_inf"
)

# Exclude samples with any of selected infections (keep "respiratory infection" cohort)
pheno_resp_infection <- na.omit(
  pheno_data_infection[-which(rowSums(pheno_data_infection[, 54:63] == 1) > 0),]
)

# Expression matrix for respiratory infection subset
gene_expr_resp_inf <- as.matrix(
  count[,which(colnames(count) %in% pheno_resp_infection$BBNId)]
)

write.csv(
  pheno_resp_infection,
  "input\\data\\Metadata\\pheno_resp_infection.csv"
)
write.csv(
  gene_expr_resp_inf,
  "input\\data\\Expression_data\\log2.cpm_resp_inf.csv"
)

############################ Outlier detection and export cleaned cohort ############################

# PCA before outlier removal
pca_res <- prcomp(t(gene_expr_resp_inf))
rownames(pheno_resp_infection) <- pheno_resp_infection$BBNId
autoplot(pca_res, data = pheno_resp_infection, colour = ("Infection"), label = TRUE, label.size = 3)

# Mahalanobis outlier detection (function defined elsewhere in pipeline)
outliers_results <- mahalanobis.outlier(gene_expr_resp_inf)

# Re-align phenotype rows to expression columns
if(!identical(colnames(gene_expr_resp_inf), pheno_resp_infection$BBNId)){
  errorCondition("Not matching IDs")
  pheno_resp_infection <- pheno_resp_infection[
    match(colnames(gene_expr_resp_inf),pheno_resp_infection$BBNId),
  ]
}

# Fix Queen Square bank naming
colnames(pheno_resp_infection)[grep("Brain.Bank_Queen",colnames(pheno_resp_infection),ignore.case = TRUE)] <- "Brain.Bank_Queen_Square"

# Export phenotype and expression matrices (with outliers retained; removal may be elsewhere)
write.csv(
  pheno_resp_infection,
  "input\\data\\Metadata\\pheno_resp_infection_no_out.csv"
)
write.csv(
  gene_expr_resp_inf,
  "input\\data\\Expression_data\\gene_expr_resp_inf_no_out.csv"
)

############################ Re-encode covariates for limma ############################

factor_list <- c(
  "AD","Infection","Gender",
  "Brain.Bank_Bristol",
  "Brain.Bank_Edinburgh",
  "Brain.Bank_Manchester",
  "Brain.Bank_Queen_Square",
  "Brain.Bank_Newcastle"
)
for(i in factor_list){
  print(i)
  pheno_resp_infection[,i] <- as.factor(pheno_resp_infection[,i])
}

# PCA after outlier handling
pca_res_resp_inf <- prcomp(t(gene_expr_resp_inf))
autoplot(pca_res_resp_inf, data = pheno_resp_infection, colour = ("Infection"), label = TRUE, label.size = 3)

############################ Differential expression analysis (limma) ############################

# Design matrix:
# AD, Infection, AD*Infection and covariates (age, gender, brain bank, cell-type proportions)
design <- model.matrix(
  ~ AD * Infection +
    Age +
    Gender + 
    Brain.Bank_Bristol +
    Brain.Bank_Edinburgh +
    Brain.Bank_Manchester +
    Brain.Bank_Queen_Square +
    #Brain.Bank_Newcastle +
    fpkm_astrocytes +
    fpkm_endothelial +
    fpkm_microglia +
    fpkm_neurons +
    fpkm_OPC +
    fpkm_oligodendrocytes,
  data = pheno_resp_infection
)
colnames(design) <- make.names(colnames(design))

# Fit linear model and empirical Bayes
fit <- lmFit(gene_expr_resp_inf, design)
fit <- eBayes(fit)

# Contrasts:

contrast.matrix <- makeContrasts(
  AD_with_vs_without_infection      = Infection1 + AD1.Infection1,
  Control_with_vs_without_infection = Infection1,
  AD_vs_Control_without_infection   = AD1,
  AD_vs_Control_with_infection      = AD1 + AD1.Infection1,
  levels = design
)

# Apply contrasts and recompute statistics
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract results
tt_AD_with_vs_without_infection      <- topTable(
  fit2,
  coef = "AD_with_vs_without_infection",
  number = Inf,
  sort.by = "p",
  adjust.method = "fdr"
)
tt_Control_with_vs_without_infection <- topTable(
  fit2,
  coef = "Control_with_vs_without_infection",
  number = Inf,
  sort.by = "none",
  adjust.method = "fdr"
)
tt_AD_vs_Control_no  <- topTable(
  fit2,
  coef = "AD_vs_Control_without_infection",
  number = Inf,
  sort.by = "none",
  adjust.method = "fdr"
)
tt_AD_vs_Control_inf <- topTable(
  fit2,
  coef = "AD_vs_Control_with_infection",
  number = Inf,
  sort.by = "none",
  adjust.method = "fdr"
)

############################ Map Ensembl IDs to HGNC symbols ############################

tt_AD_vs_Control_inf_symbol <- convertToSymbol(rownames(tt_AD_vs_Control_inf))
tt_AD_vs_Control_inf_symbol <- merge(
  tt_AD_vs_Control_inf,
  tt_AD_vs_Control_inf_symbol,
  by.x = "row.names",
  by.y = "ensembl_gene_id_version"
)
names(tt_AD_vs_Control_inf_symbol)[names(tt_AD_vs_Control_inf_symbol) == "Row.names"] <- "ensembl_gene_id_version"

tt_AD_vs_Control_no_symbol <- convertToSymbol(rownames(tt_AD_vs_Control_no))
tt_AD_vs_Control_no_symbol <- merge(
  tt_AD_vs_Control_no,
  tt_AD_vs_Control_no_symbol,
  by.x = "row.names",
  by.y = "ensembl_gene_id_version"
)
names(tt_AD_vs_Control_no_symbol)[names(tt_AD_vs_Control_no_symbol) == "Row.names"] <- "ensembl_gene_id_version"

tt_AD_with_vs_without_infection_symbol <- convertToSymbol(rownames(tt_AD_with_vs_without_infection))
tt_AD_with_vs_without_infection_symbol <- merge(
  tt_AD_with_vs_without_infection,
  tt_AD_with_vs_without_infection_symbol,
  by.x = "row.names",
  by.y = "ensembl_gene_id_version"
)
names(tt_AD_with_vs_without_infection_symbol)[names(tt_AD_with_vs_without_infection_symbol) == "Row.names"] <- "ensembl_gene_id_version"

tt_Control_with_vs_without_infection_symbol <- convertToSymbol(rownames(tt_Control_with_vs_without_infection))
tt_Control_with_vs_without_infection_symbol <- merge(
  tt_Control_with_vs_without_infection,
  tt_Control_with_vs_without_infection_symbol,
  by.x = "row.names",
  by.y = "ensembl_gene_id_version"
)
names(tt_Control_with_vs_without_infection_symbol)[names(tt_Control_with_vs_without_infection_symbol) == "Row.names"] <- "ensembl_gene_id_version"

############################ Extract significant genes ############################

sign_AD_with_vs_without_infection      <- tt_AD_with_vs_without_infection_symbol[tt_AD_with_vs_without_infection_symbol$adj.P.Val < 0.1 ,]
sign_Control_with_vs_without_infection <- tt_Control_with_vs_without_infection_symbol[tt_Control_with_vs_without_infection_symbol$adj.P.Val < 0.1,]
sign_AD_vs_Control_no  <- tt_AD_vs_Control_no_symbol[tt_AD_vs_Control_no_symbol$adj.P.Val < 0.05, ]
sign_AD_vs_Control_inf <- tt_AD_vs_Control_inf_symbol[tt_AD_vs_Control_inf_symbol$adj.P.Val < 0.05,]

# Add Ensembl IDs without version
tt_AD_with_vs_without_infection_symbol$ensembl_id <- unlist(removeVersion(tt_AD_with_vs_without_infection_symbol$ensembl_gene_id_version))
tt_Control_with_vs_without_infection_symbol$ensembl_id <- unlist(removeVersion(tt_Control_with_vs_without_infection_symbol$ensembl_gene_id_version))
sign_AD_vs_Control_no$ensembl_id  <- unlist(removeVersion(sign_AD_vs_Control_no$ensembl_gene_id_version))
sign_AD_vs_Control_inf$ensembl_id <- unlist(removeVersion(sign_AD_vs_Control_inf$ensembl_gene_id_version))

# Export significant gene lists for downstream analysis
write.table(sign_AD_vs_Control_no, "sncEnrichment\\input\\data\\sign_AD_vs_Control_no.txt",  row.names = FALSE, quote = FALSE)
write.table(sign_AD_vs_Control_inf,"sncEnrichment\\input\\data\\sign_AD_vs_Control_inf.txt", row.names = FALSE, quote = FALSE)

############################ Combine results across contrasts ############################

# Integrate effect sizes and p-values into a single table
results_combined <- data.frame(
  Gene                                 = row.names(tt_AD_vs_Control_inf),
  logFC_AD_vs_Control_inf              = tt_AD_vs_Control_inf$logFC,
  pvalue_AD_vs_Control_inf             = tt_AD_vs_Control_inf$adj.P.Val,
  logFC_AD_with_vs_without_infection   = tt_AD_with_vs_without_infection$logFC,
  pvalue_AD_with_vs_without_infection  = tt_AD_with_vs_without_infection$adj.P.Val,
  logFC_Control_with_vs_without_infection  = tt_Control_with_vs_without_infection$logFC,
  pvalue_Control_with_vs_without_infection = tt_Control_with_vs_without_infection$adj.P.Val,
  logFC_AD_vs_Control_no               = tt_AD_vs_Control_no$logFC,
  pvalue_AD_vs_Control_no              = tt_AD_vs_Control_no$adj.P.Val
)

results_sorted <- results_combined[order(results_combined$pvalue_AD_vs_Control_inf), ]

results_combined_symbol <- convertToSymbol(results_combined$Gene)

# Merge with annotation (Ensembl/HGNC)
results_combined_annot <- merge(
  results_combined_symbol,
  results_combined,
  by.x = "ensembl_gene_id_version",
  by.y = "Gene"
)

ensembl_gene_ids <- removeVersion(results_combined_annot$ensembl_gene_id_version) 
results_combined_annot$ensembl_gene_id <- unlist(ensembl_gene_ids)
write.csv(results_combined_annot, "output\\data\\results_limma.csv")

# Top 100 genes for infection contrasts (by P.Value)
top100_AD_with_vs_without_infection <- tt_AD_with_vs_without_infection_symbol[order(tt_AD_with_vs_without_infection_symbol$P.Value), ][1:100, ]
top100_Control_with_vs_without_infection <- tt_Control_with_vs_without_infection_symbol[order(tt_Control_with_vs_without_infection_symbol$P.Value), ][1:100, ]

############################ Save key result tables ############################

write.csv(sign_AD_vs_Control_inf, 
          "output\\data\\AD_vs_Control_inf_sign_results.csv", quote = FALSE)
write.csv(top100_AD_with_vs_without_infection, 
          "output\\data\\AD_with_vs_without_infection_top100_results.csv", quote = FALSE, row.names = FALSE)
write.csv(top100_Control_with_vs_without_infection, 
          "output\\data\\Control_with_vs_without_infection_top100_results.csv", quote = FALSE, row.names = FALSE)
write.csv(sign_AD_vs_Control_no, 
          "output\\data\\AD_vs_Control_no_sign_results.csv", quote = FALSE, row.names = FALSE)

############################ Venn diagram: AD vs Control (with/without infection) ############################

tiff("output\\figures\\Venn_diagram\\venn_limma_results.tif",
     units = "in", width = 12, height = 8, res = 600, compression = "lzw")

grid.newpage()
# Space for labels
pushViewport(viewport(
  x = 0.5, y = 0.46,
  width = 0.8,
  height = 0.96,
  clip = "off"
))

venn.plot <- draw.pairwise.venn(
  area1      = length(sign_AD_vs_Control_no$ensembl_gene_id_version),
  area2      = length(sign_AD_vs_Control_inf$ensembl_gene_id_version),
  cross.area = length(intersect(
    sign_AD_vs_Control_no$ensembl_gene_id_version,
    sign_AD_vs_Control_inf$ensembl_gene_id_version
  )),
  category = c("AD vs. Control no infection", "AD vs. Control with infection"),
  fill     = c("#0072B2", "#E69F00"),
  alpha    = 0.5,
  lty      = "blank",
  cex      = 1.5,
  cat.cex  = 1.2,
  cat.pos  = c(-20, 20)
)
popViewport()
dev.off()

############################ Venn diagram: up/down genes across infection status ############################

# Count genes passing logFC and FDR thresholds
nrow(tt_AD_vs_Control_inf_symbol[abs(tt_AD_vs_Control_inf_symbol$logFC) > log2(1.1) &
                                   tt_AD_vs_Control_inf_symbol$adj.P.Val < 0.05,])

up_genes_AD_vs_Control_inf <- tt_AD_vs_Control_inf_symbol[
  tt_AD_vs_Control_inf_symbol$logFC > log2(1.1) & 
    tt_AD_vs_Control_inf_symbol$adj.P.Val < 0.05,
  "hgnc_symbol"
]

down_genes_AD_vs_Control_inf <- tt_AD_vs_Control_inf_symbol[
  tt_AD_vs_Control_inf_symbol$logFC < -log2(1.1) & 
    tt_AD_vs_Control_inf_symbol$adj.P.Val < 0.05,
  "hgnc_symbol"
]

up_genes_AD_vs_Control_no <- tt_AD_vs_Control_no_symbol[
  tt_AD_vs_Control_no_symbol$logFC > log2(1.1) & 
    tt_AD_vs_Control_no_symbol$adj.P.Val < 0.05,
  "hgnc_symbol"
]
down_genes_AD_vs_Control_no <- tt_AD_vs_Control_no_symbol[
  tt_AD_vs_Control_no_symbol$logFC < -log2(1.1) & 
    tt_AD_vs_Control_no_symbol$adj.P.Val < 0.05,
  "hgnc_symbol"
]

# Remove empty symbols
up_genes_AD_vs_Control_inf    <- up_genes_AD_vs_Control_inf[up_genes_AD_vs_Control_inf    != ""]
down_genes_AD_vs_Control_inf  <- down_genes_AD_vs_Control_inf[down_genes_AD_vs_Control_inf != ""]
up_genes_AD_vs_Control_no     <- up_genes_AD_vs_Control_no[up_genes_AD_vs_Control_no     != ""]
down_genes_AD_vs_Control_no   <- down_genes_AD_vs_Control_no[down_genes_AD_vs_Control_no != ""]

venn.plot <- draw.quad.venn(
  area1 = length(up_genes_AD_vs_Control_inf),
  area2 = length(down_genes_AD_vs_Control_inf),
  area3 = length(up_genes_AD_vs_Control_no),
  area4 = length(down_genes_AD_vs_Control_no),
  n12   = length(intersect(up_genes_AD_vs_Control_inf, down_genes_AD_vs_Control_inf)),
  n13   = length(intersect(up_genes_AD_vs_Control_inf, up_genes_AD_vs_Control_no)),
  n14   = length(intersect(up_genes_AD_vs_Control_inf, down_genes_AD_vs_Control_no)),
  n23   = length(intersect(down_genes_AD_vs_Control_inf, up_genes_AD_vs_Control_no)),
  n24   = length(intersect(down_genes_AD_vs_Control_inf, down_genes_AD_vs_Control_no)),
  n34   = length(intersect(up_genes_AD_vs_Control_no, down_genes_AD_vs_Control_no)),
  n123  = length(Reduce(intersect, list(
    up_genes_AD_vs_Control_inf, down_genes_AD_vs_Control_inf, up_genes_AD_vs_Control_no
  ))),
  n124  = length(Reduce(intersect, list(
    up_genes_AD_vs_Control_inf, down_genes_AD_vs_Control_inf, down_genes_AD_vs_Control_no
  ))),
  n134  = length(Reduce(intersect, list(
    up_genes_AD_vs_Control_inf, up_genes_AD_vs_Control_no, down_genes_AD_vs_Control_no
  ))),
  n234  = length(Reduce(intersect, list(
    down_genes_AD_vs_Control_inf, up_genes_AD_vs_Control_no, down_genes_AD_vs_Control_no
  ))),
  n1234 = length(Reduce(intersect, list(
    up_genes_AD_vs_Control_inf, down_genes_AD_vs_Control_inf,
    up_genes_AD_vs_Control_no,  down_genes_AD_vs_Control_no
  ))),
  category = c(
    "Up AD vs Ctrl (inf)", "Down AD vs Ctrl (inf)", 
    "Up AD vs Ctrl (no inf)", "Down AD vs Ctrl (no inf)"
  ),
  fill = c("#D55E00", "#009E73", "#0072B2", "#CC79A7"),
  alpha = 0.5,
  cat.pos = c(-10, 10, 0, 0),
  cat.dist = c(0.2, 0.2, 0.15, 0.15),
  cat.just = list(c(0.5, 0.5), c(0.5, 0.5), c(0.5, 0.5), c(0.5, 0.5)),
  cex = 2,
  cat.cex = 1.5,
  scaled = FALSE
)

############################ Volcano plot: AD vs Control with infection ############################

# Custom colour key based on FC and FDR (fc_cutoff must be defined upstream)
keyvals <- ifelse(
  (tt_AD_vs_Control_inf_symbol$adj.P.Val > 0.05 | abs(tt_AD_vs_Control_inf_symbol$logFC) < fc_cutoff),
  "grey",
  ifelse(tt_AD_vs_Control_inf_symbol$logFC < -fc_cutoff, "blue", "red")
)
keyvals[is.na(keyvals)] <- "grey"

names(keyvals)[keyvals == "red"]  <- "Upregulated"
names(keyvals)[keyvals == "blue"] <- "Downregulated"
names(keyvals)[keyvals == "grey"] <- "Not Significant"

tt_AD_vs_Control_inf_symbol <- tt_AD_vs_Control_inf_symbol[order(tt_AD_vs_Control_inf_symbol$adj.P.Val),]
tt_AD_vs_Control_no_symbol  <- tt_AD_vs_Control_no_symbol[order(tt_AD_vs_Control_no_symbol$adj.P.Val),]

top10 <- head(tt_AD_vs_Control_inf_symbol$hgnc_symbol, 10)

# EnhancedVolcano for AD vs Ctrl (with infection)
p <- EnhancedVolcano(
  tt_AD_vs_Control_inf_symbol,
  lab        = tt_AD_vs_Control_inf_symbol$hgnc_symbol,
  selectLab  = top10,
  x          = "logFC",
  y          = "adj.P.Val",
  pCutoff    = 0.05,
  FCcutoff   = 0,
  colCustom  = keyvals,
  legendLabels = c("Not Significant", "Downregulated", "Upregulated"),
  title      = NULL,
  subtitle   = "",
  pointSize  = 4.0,
  labSize    = 10,
  labCol     = "black",
  labFace    = "plain",
  drawConnectors = TRUE, arrowheads = FALSE,
  widthConnectors = 1,
  xlim       = c(-1.4, 1.4), 
  ylim       = c(0, 3.5),
  axisLabSize = 24,
  legendLabSize = 24,
  legendIconSize = 8
)

tiff("RNA_sequencing\\output\\figures\\volcano_plot\\volcano_plot_ADvsCTR_infection.tiff", 
     height = 1200, width = 1200)
p + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")
dev.off()

############################ Volcano plot: AD vs Control without infection ############################

keyvals1 <- ifelse(
  (tt_AD_vs_Control_no_symbol$adj.P.Val > 0.05 | abs(tt_AD_vs_Control_no_symbol$logFC) < fc_cutoff),
  "grey",
  ifelse(tt_AD_vs_Control_no_symbol$logFC < -fc_cutoff, "blue", "red")
)
keyvals1[is.na(keyvals1)] <- "grey"

names(keyvals1)[keyvals1 == "red"]  <- "Upregulated"
names(keyvals1)[keyvals1 == "blue"] <- "Downregulated"
names(keyvals1)[keyvals1 == "grey"] <- "Not Significant"

# Build labels prioritising HGNC symbol; fallback to Ensembl without version
tt_AD_vs_Control_no_symbol$labels <- tt_AD_vs_Control_no_symbol %>%
  arrange(adj.P.Val) %>%
  mutate(
    gene_final = ifelse(
      hgnc_symbol == "" | is.na(hgnc_symbol),
      str_remove(ensembl_gene_id_version, "\\..*"),
      hgnc_symbol
    )
  ) %>%
  pull(gene_final)

top10_no_inf <- head(tt_AD_vs_Control_no_symbol$labels, 10)

p1 <- EnhancedVolcano(
  tt_AD_vs_Control_no_symbol,
  lab        = tt_AD_vs_Control_no_symbol$labels,
  x          = "logFC",
  y          = "adj.P.Val",
  selectLab  = top10_no_inf,
  pCutoff    = 0.05,
  FCcutoff   = 0,
  colCustom  = keyvals1,
  legendLabels = c("Not Significant", "Downregulated", "Upregulated"),
  title      = NULL,
  subtitle   = "",
  pointSize  = 4.0,
  labSize    = 10,
  labCol     = "black",
  labFace    = "plain",
  drawConnectors = TRUE, arrowheads = FALSE,
  widthConnectors = 1,
  xlim       = c(-1.5, 1.5), 
  ylim       = c(0, 6),
  axisLabSize = 24,
  legendLabSize = 24,
  legendIconSize = 8
)

tiff("output\\figures\\volcano_plot\\volcano_plot_ADvsCTR_no_infection.tiff", 
     height = 1200, width = 1200)
p1 + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")
dev.off()
