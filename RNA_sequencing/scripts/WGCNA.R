###########################################
# @author: Giulia Pegoraro
#
# Weighted Gene Co-expression Network Analysis (WGCNA)
#
###########################################

################################################################################
################################### LIBRARIES ##################################
################################################################################

if (!requireNamespace("WGCNA", quietly = TRUE))
  BiocManager::install("WGCNA")

library(WGCNA)
library(magrittr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggvenn)
library(VennDiagram)
library(grid)
library(fastDummies)

################################################################################
################################### FUNCTIONS ##################################
################################################################################

############################################
# Function: calculate_power
#
# Description:
#   Estimate the soft-thresholding power for WGCNA by evaluating scale-free
#   topology and mean connectivity over a range of candidate powers.
#
# Arguments:
#   expr_matrix : numeric matrix/data frame of expression values,
#                 genes in rows, samples in columns.
#
# Returns:
#   sft         : list returned by pickSoftThreshold containing the
#                 fit indices for each power.
############################################
calculate_power <- function(expr_matrix){
  allowWGCNAThreads()          # allow multi-threading (optional)
  
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
  
  # Network topology analysis
  sft = pickSoftThreshold(
    t(expr_matrix), 
    powerVector = powers,
    verbose = 5
  )
  
  par(mfrow = c(1,2))
  cex1 = 0.9
  
  # Scale-free topology fit
  plot(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    xlab = "Soft Threshold (power)",
    ylab = "Scale Free Topology Model Fit, signed R^2",
    main = "Scale independence"
  )
  text(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    labels = powers, cex = cex1, col = "red"
  )
  abline(h = 0.90, col = "red")
  
  # Mean connectivity
  plot(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    xlab = "Soft Threshold (power)",
    ylab = "Mean Connectivity",
    type = "n",
    main = "Mean connectivity"
  )
  text(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    labels = powers,
    cex = cex1, col = "red"
  )
  
  return(sft)
}

############################################
# Function: perform_WCGNA
#
# Description:
#   Run blockwise WGCNA on an expression matrix and return the network object
#   containing modules and associated information.
#
# Arguments:
#   expr_matrix     : numeric matrix/data frame of expression values,
#                     genes in rows, samples in columns.
#   power           : integer, soft-thresholding power.
#   max_Block_Size  : integer, maximum block size for blockwiseModules.
#
# Returns:
#   netwk           : list containing the WGCNA network and module information.
############################################
perform_WCGNA <- function(expr_matrix, power, max_Block_Size){
  picked_power = power
  temp_cor <- cor       
  cor <- WGCNA::cor   # Use WGCNA cor function (fix namespace conflict)
  
  netwk <- blockwiseModules(
    t(expr_matrix),
    
    # Adjacency / topology
    power = picked_power,
    networkType = "signed",
    
    # Tree and block options
    deepSplit = 2,
    pamRespectsDendro = FALSE,
    minModuleSize = 50,
    maxBlockSize = max_Block_Size,
    
    # Module adjustments
    reassignThreshold = 0,
    mergeCutHeight = 0.15,
    
    # TOM options (archive results)
    saveTOMs = TRUE,
    saveTOMFileBase = "ER",
    
    # Output options
    numericLabels = TRUE,
    verbose = 3
  )
  
  cor <- temp_cor
  return(netwk)
}

############################################
# Function: plot_dendogram
#
# Description:
#   Plot gene dendrogram with module colors annotated underneath.
#
# Arguments:
#   netwk        : WGCNA network object returned by blockwiseModules.
#   mergedColors : vector of module colors, one per gene.
#
# Returns:
#   None (plot is generated).
############################################
plot_dendogram <- function(netwk, mergedColors){
  plotDendroAndColors(
    netwk$dendrograms[[1]],
    mergedColors[netwk$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05
  )
}

############################################
# Function: calculateResiduals
#
# Description:
#   Fit a linear model for each gene and return residuals after adjusting
#   for a given set of covariates.
#
# Arguments:
#   expr_matrix : numeric matrix of expression data (genes x samples).
#   pheno       : data frame with phenotype / covariate data (samples x covs).
#   covs        : character vector of column names in 'pheno' to include
#                 as covariates in the linear model.
#
# Returns:
#   residuals   : numeric matrix of residual expression values (genes x samples).
############################################
calculateResiduals <- function(expr_matrix, pheno, covs){
  form_vars <- paste(covs, collapse = " + ")
  
  residuals <- as.matrix(
    do.call(
      rbind,
      lapply(rownames(expr_matrix), function(row){
        formula_lm <- paste("expr_matrix[row,] ~ ", form_vars)
        fit <- lm(
          formula(formula_lm),
          data = pheno
        )
        return(fit$residuals)
      })
    )
  )
  
  return(residuals)
}

############################################
# Function: heatmap.WGCNA
#
# Description:
#   Create module–trait correlation heatmap and dendrogram plot, and
#   return module assignments, eigengenes, and trait matrix.
#
# Arguments:
#   netwk       : WGCNA network object returned by blockwiseModules.
#   pheno_data  : data frame of phenotype data (rows = samples).
#   residuals   : matrix of residual expression values (genes x samples).
#   trait_list  : character vector of trait column names in 'pheno_data'
#                 to correlate with modules.
#   path        : character, path prefix to save plots.
#   title       : character, file name for the heatmap (unused if dev.off commented).
#
# Returns:
#   list(
#     module_df : data frame with gene IDs and module colors,
#     MEs0      : matrix of module eigengenes (samples x modules),
#     trait_df  : numeric trait matrix aligned to MEs0
#   )
############################################
heatmap.WGCNA <- function(netwk, pheno_data, residuals, trait_list, path, title){
  # Convert labels to colors for plotting
  mergedColors = labels2colors(netwk$colors)
  
  # Plot dendrogram with module colors
  png(paste0(path, "dendpogram.png"), width = 3000, height = 2000, res = 300)
  plot_dendogram(netwk, mergedColors)
  dev.off()
  
  # Data frame with module assignment per gene
  module_df <- data.frame(
    gene_id = names(netwk$colors),
    colors  = labels2colors(netwk$colors)
  )
  
  # Module eigengenes for residual expression
  MEs0 <- moduleEigengenes(t(residuals), mergedColors)$eigengenes
  
  # Reorder modules so similar modules are next to each other
  MEs0 <- orderMEs(MEs0)
  module_order <- names(MEs0) %>% gsub("ME", "", .)
  colnames(MEs0) <- module_order
  
  # Trait matrix
  trait_df <- pheno_data[, trait_list]
  trait_df <- apply(trait_df, 2, as.numeric)
  row.names(trait_df) <- pheno_data$BBNId
  
  # Align trait matrix and eigengenes
  if (!identical(rownames(trait_df), rownames(MEs0))){
    print("IDs have different order")
    trait_df <- trait_df[match(rownames(MEs0), rownames(trait_df)), ]
  }
  
  # Remove grey module (unassigned genes)
  MEs0 <- MEs0[, -which(colnames(MEs0) == "grey")]
  
  # Correlation between module eigengenes and traits
  moduleTraitCor    <- stats::cor(MEs0, trait_df, use = "p", method = "pearson")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(MEs0))
  
  # Prepare matrix with correlation and p-value
  textMatrix <- paste(
    signif(moduleTraitCor, 2),
    "\n(",
    signif(moduleTraitPvalue, 1),
    ")",
    sep = ""
  )
  dim(textMatrix) <- dim(moduleTraitCor)
  
  # Plot heatmap
  par(mar = c(10, 10, 5, 5), cex.main = 1.5, cex.axis = 1, cex.lab = 0.7)
  labeledHeatmap(
    Matrix       = moduleTraitCor,
    xLabels      = colnames(trait_df),
    yLabels      = names(MEs0),
    ySymbols     = names(MEs0),
    colorLabels  = FALSE,
    colors       = colorRampPalette(c("blue", "white", "red"))(50),
    textMatrix   = textMatrix,
    setStdMargins = TRUE,
    plotLegend   = TRUE,
    cex.text     = 0.7,
    zlim         = c(-1, 1),
    main         = "Module-trait relationships"
  )
  
  return(list(module_df, MEs0, trait_df))
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
# Function: moduleMembership
#
# Description:
#   Compute module membership (MM) and gene-trait significance (GS) for a
#   given trait, generate MM–GS scatterplots per module, and identify hub
#   genes based on quantile thresholds of MM and GS.
#
# Arguments:
#   trait                : character, name of the trait (column in trait_df).
#   network              : WGCNA network object.
#   trait_df             : data frame with trait values (rows = samples).
#   MEs0                 : matrix of module eigengenes (samples x modules).
#   residuals            : matrix of residual expression (genes x samples).
#   path                 : character, directory to save MM–GS plots and hub genes.
#   quantile_threshold_GS: numeric, GS quantile used as cutoff (default 0.75).
#   quantile_threshold_MM: numeric, MM quantile used as cutoff (default 0.75).
#
# Returns:
#   moduleGeneSignificance : data frame containing MM and GS values for each
#                            gene; has an attribute "hub_genes" which is a list
#                            of hub gene IDs per module.
############################################
moduleMembership <- function(
    trait,
    network,
    trait_df,
    MEs0,
    residuals,
    path,
    quantile_threshold_GS = 0.75,
    quantile_threshold_MM = 0.75
) {
  
  weight <- as.data.frame(trait_df[, trait])
  names(weight) <- "weight"
  
  # Module names (colors)
  modNames <- names(MEs0)
  
  # Gene-module membership
  geneModuleMembership <- as.data.frame(cor(t(residuals), MEs0, use = "p"))
  MMPvalue             <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), ncol(residuals)))
  
  # Gene-trait significance
  geneTraitSignificance <- as.data.frame(cor(t(residuals), weight, use = "p"))
  GSPvalue              <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), ncol(residuals)))
  
  # Convert Ensembl ID (with version) to plain Ensembl ID
  gene_list  <- sapply(strsplit(row.names(GSPvalue), "\\."), `[`, 1)
  GSPvalue_id <- data.frame(gene_list, GSPvalue)
  
  mergedColors <- as.matrix(labels2colors(network$colors))
  colorlevels  <- unique(mergedColors[, 1])
  
  hub_genes_list <- list()  # store hub genes for each module
  
  for (module in colorlevels) {
    png(paste0(path, trait, "/mm_", module, ".png"))
    
    # Match module name to modNames
    column <- match(module, modNames)
    if (is.na(column)) {
      warning(paste("Skipping module:", module))
      next
    }
    
    moduleGenes <- mergedColors[, 1] == module
    
    x <- geneModuleMembership[moduleGenes, column]
    if (!is.numeric(x)) x <- as.numeric(as.character(x))
    
    y <- geneTraitSignificance[moduleGenes, 1]
    if (!is.numeric(y)) y <- as.numeric(as.character(y))
    
    verboseScatterplot(
      abs(x),
      abs(y),
      xlab = paste("Module Membership in", module, "module"),
      ylab = paste0("Gene significance for ", trait),
      main = "Module membership vs. gene significance\n",
      cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module
    )
    
    # Quantile thresholds for MM and GS
    mm_cutoff <- quantile(abs(x), probs = quantile_threshold_MM, na.rm = TRUE)
    print(mm_cutoff)
    gs_cutoff <- quantile(abs(y), probs = quantile_threshold_GS, na.rm = TRUE)
    print(gs_cutoff)
    
    # Threshold lines
    abline(h = gs_cutoff, col = "red",  lty = 2)
    abline(v = mm_cutoff, col = "blue", lty = 2)
    
    dev.off()
    
    # Select hub genes
    hub_genes <- rownames(geneModuleMembership)[moduleGenes][
      abs(x) >= mm_cutoff & abs(y) >= gs_cutoff
    ]
    hub_genes_list[[module]] <- hub_genes
    
    # Save hub genes to file
    write.table(
      hub_genes,
      file = paste0(path, "hub_genes_", module, ".txt"),
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
  }
  
  # Combine MM and GS in a single object
  moduleGeneSignificance <- merge(geneModuleMembership, geneTraitSignificance, by = "row.names")
  
  # Attach hub genes
  attr(moduleGeneSignificance, "hub_genes") <- hub_genes_list
  
  return(moduleGeneSignificance)
}

################################################################################
##################################### MAIN #####################################
################################################################################

############################ Respiratory infection cohort ############################
# Import data and convert data types

pheno_resp_infection <- read.csv("input\\data\\Metadata\\pheno_resp_infection_no_out.csv", row.names = 1)
gene_expr_resp_inf   <- as.matrix(read.csv("input\\data\\Expression_data\\gene_expr_resp_inf_no_out.csv", row.names = 1))
limma_AD_CTR_with_infection <- read.table(
  "input\\data\\sign_AD_vs_Control_inf.txt",
  header = TRUE
)

factor_list <- c(
  "AD", "Infection", "Gender",
  "Brain.Bank_Bristol", "Brain.Bank_Edinburgh",
  "Brain.Bank_Manchester", "Brain.Bank_Queen_Square",
  "Brain.Bank_Newcastle", "Batch"
)

numeric_list <- c(
  "Age", "PMD", "TIN.median.",
  "fpkm_astrocytes", "fpkm_endothelial", "fpkm_microglia",
  "fpkm_neurons", "fpkm_OPC", "fpkm_oligodendrocytes",
  "fpkm_fetal_quiescent", "fpkm_fetal_replicating"
)

# Expand Group into dummy variables
pheno_resp_infection <- dummy_cols(
  pheno_resp_infection,
  select_columns = c("Group")
)

colnames(pheno_resp_infection)[colnames(pheno_resp_infection) == "Group_0"]  <- "CTR No Infection"
colnames(pheno_resp_infection)[colnames(pheno_resp_infection) == "Group_1"]  <- "CTR with Infection"
colnames(pheno_resp_infection)[colnames(pheno_resp_infection) == "Group_10"] <- "AD without Infection"
colnames(pheno_resp_infection)[colnames(pheno_resp_infection) == "Group_11"] <- "AD with Infection"

# Convert factor and numeric variables to appropriate types
for (i in factor_list){
  print(i)
  pheno_resp_infection[, i] <- as.factor(pheno_resp_infection[, i])
}

for (i in numeric_list){
  print(i)
  pheno_resp_infection[, i] <- as.numeric(pheno_resp_infection[, i])
}

############################ Cell proportion regression ############################
# Regress out cell-type proportions and other covariates

covs_cell_no_fetal <- c(
  "Age", "Gender",
  "Brain.Bank_Bristol", "Brain.Bank_Edinburgh",
  "Brain.Bank_Manchester", "Brain.Bank_Queen_Square", "Brain.Bank_Newcastle",
  "fpkm_astrocytes", "fpkm_endothelial", "fpkm_microglia",
  "fpkm_neurons", "fpkm_OPC", "fpkm_oligodendrocytes"
)

residuals_cell_no_fetal <- calculateResiduals(
  gene_expr_resp_inf,
  pheno_resp_infection,
  covs_cell_no_fetal
)

row.names(residuals_cell_no_fetal) <- row.names(gene_expr_resp_inf)
colnames(residuals_cell_no_fetal) <- colnames(gene_expr_resp_inf)

# Soft-threshold selection and network construction
power_cell_no_fetal <- calculate_power(residuals_cell_no_fetal)
netwk_cell_no_fetal <- perform_WCGNA(residuals_cell_no_fetal, 5, 20000)

trait_cell <- c("AD", "Infection")

module_df_cell_no_fetal <- heatmap.WGCNA(
  netwk_cell_no_fetal,
  pheno_resp_infection,
  residuals_cell_no_fetal,
  trait_cell,
  "output\\figures\\module_trait_matrix",
  "mtr_cell_nofetal.png"
)

save("module_df_cell_no_fetal.RData")

# Export module gene lists for selected colors
write.csv(
  module_df_cell_no_fetal[[1]][module_df_cell_no_fetal[[1]]$colors == "tan", "gene_id"],
  "sncEnrichment\\input\\data\\tan_module.csv"
)
write.csv(
  module_df_cell_no_fetal[[1]][module_df_cell_no_fetal[[1]]$colors == "yellow", "gene_id"],
  "sncEnrichment\\input\\data\\yellow_module.csv"
)
write.csv(
  module_df_cell_no_fetal[[1]][module_df_cell_no_fetal[[1]]$colors == "red", "gene_id"],
  "sncEnrichment\\input\\data\\red_module.csv"
)

# Overlap between modules and DE genes
color_intersection <- lapply(
  unique(module_df_cell_no_fetal[[1]]$colors),
  function(color){
    intersect(
      c(limma_AD_CTR_with_infection$x),
      c(module_df_cell_no_fetal[[1]][module_df_cell_no_fetal[[1]]$colors == color, "gene_id"])
    )
  }
)

############################ Convert Ensembl IDs to HGNC symbols ############################

red_module_symbol    <- convertToSymbol(module_df_cell_no_fetal[[1]][module_df_cell_no_fetal[[1]]$colors == "red",    "gene_id"])
yellow_module_symbol <- convertToSymbol(module_df_cell_no_fetal[[1]][module_df_cell_no_fetal[[1]]$colors == "yellow", "gene_id"])
tan_module_symbol    <- convertToSymbol(module_df_cell_no_fetal[[1]][module_df_cell_no_fetal[[1]]$colors == "tan",    "gene_id"])

red_module_symbol$ensembl_gene_id    <- unlist(removeVersion(red_module_symbol$ensembl_gene_id_version))
yellow_module_symbol$ensembl_gene_id <- unlist(removeVersion(yellow_module_symbol$ensembl_gene_id_version))
tan_module_symbol$ensembl_gene_id    <- unlist(removeVersion(tan_module_symbol$ensembl_gene_id_version))

write.csv(yellow_module_symbol, "output//data//yellow_module_symbol.csv", row.names = FALSE)
write.csv(tan_module_symbol,    "output//data//tan_module_symbol.csv",    row.names = FALSE)
write.csv(red_module_symbol,    "output//data//red_module_symbol.csv",    row.names = FALSE)

write.csv(red_module_symbol$ensembl_gene_id,                           "output//data//yellow_module_ensembl.csv", row.names = FALSE)
write.csv(removeVersion(tan_module_symbol$ensembl_gene_id_version),    "output//data//tan_module_ensembl.csv",    row.names = FALSE)
write.csv(tan_module_symbol$ensembl_gene_id,                           "output//data//red_module_ensembl.csv",    row.names = FALSE)

############################ ANOVA on module eigengenes ############################

MEs0_resp_inf_tot <- module_df_cell_no_fetal[[2]]

if (!identical(rownames(MEs0_resp_inf_tot), pheno_resp_infection$BBNId)){
  errorCondition("Not matching IDs")
  pheno_resp_infection <- pheno_resp_infection[
    match(rownames(MEs0_resp_inf_tot), pheno_resp_infection$BBNId),
  ]
}

anova_eigengenes_p_value <- lapply(colnames(MEs0_resp_inf_tot), function(col) {
  anova(aov(MEs0_resp_inf_tot[, col] ~ AD * Infection, data = pheno_resp_infection))
})
names(anova_eigengenes_p_value) <- colnames(MEs0_resp_inf_tot)

############################ Venn diagrams with WGCNA and DE results ############################

results_WGCNA_limma <- list(
  yellow_module = module_df_cell_no_fetal[[1]][module_df_cell_no_fetal[[1]]$colors == "yellow", "gene_id"],
  DE_genes     = limma_AD_CTR_with_infection$x,
  tan_module   = module_df_cell_no_fetal[[1]][module_df_cell_no_fetal[[1]]$colors == "tan", "gene_id"]
)

A <- unique(results_WGCNA_limma$yellow_module)
B <- unique(results_WGCNA_limma$DE_genes)
C <- unique(results_WGCNA_limma$tan_module)

png("output/figures/Venn_diagram/WGCNA_ANOVA_venn.png", width = 4500, height = 4500, res = 300)
grid.newpage()
draw.triple.venn(
  area1   = length(A),
  area2   = length(B),
  area3   = length(C),
  n12     = length(intersect(A, B)),
  n23     = length(intersect(B, C)),
  n13     = length(intersect(A, C)),
  n123    = length(Reduce(intersect, list(A, B, C))),
  category = c("Yellow module", "DEG AD vs. CTR with infection", "Tan module"),
  fill     = c("#0073C2FF", "green", "pink"),
  lty      = "blank",
  cex      = 3,
  cat.cex  = 1.2,
  euler.d  = FALSE,
  scaled   = FALSE
)
dev.off()

# Venn diagram for two sets (green module vs DE genes)
tiff("output/figures/Venn_diagram/venn_green_limma.tif", height = 1500, width = 2000)
venn.plot <- draw.pairwise.venn(
  area1      = length(module_df_cell_no_fetal[[1]][module_df_cell_no_fetal[[1]]$colors == "green", "gene_id"]),
  area2      = length(limma_AD_CTR_with_infection$ensembl_gene_id_version),
  cross.area = length(
    intersect(
      module_df_cell_no_fetal[[1]][module_df_cell_no_fetal[[1]]$colors == "green", "gene_id"],
      limma_AD_CTR_with_infection$ensembl_gene_id_version
    )
  ),
  category = c("Green module", "DE genes AD vs CTR with infection"),
  fill     = c("skyblue", "yellow"),
  alpha    = 0.5,
  cat.pos  = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  scaled   = FALSE,
  print.mode = c("raw", "none", "raw"),
  cex      = 5,
  cat.cex  = 4
)
dev.off()

# Intersection of gene symbols in common
intersect_genes <- intersect(green_module_symbol$hgnc_symbol, limma_AD_CTR_with_infection$hgnc_symbol)

# Annotate intersection with gene symbols (coordinates may need adjustment)
grid.text(
  paste(intersect_genes, collapse = "\n"), 
  x = 0.55, y = 0.5, 
  gp = gpar(fontsize = 15)
)
dev.off()

############################ Distribution of eigengenes across groups ############################

pheno_MEs0_resp_inf_tot <- merge(
  pheno_resp_infection,
  MEs0_resp_inf_tot,
  by.x = "BBNId",
  by.y = "row.names"
)

pheno_MEs0_resp_inf_tot$Group_condition <- ""
pheno_MEs0_resp_inf_tot[pheno_MEs0_resp_inf_tot$Group %in% "0",  "Group_condition"] <- "CTR"
pheno_MEs0_resp_inf_tot[pheno_MEs0_resp_inf_tot$Group %in% "1",  "Group_condition"] <- "CTR with infection"
pheno_MEs0_resp_inf_tot[pheno_MEs0_resp_inf_tot$Group %in% "10", "Group_condition"] <- "AD"
pheno_MEs0_resp_inf_tot[pheno_MEs0_resp_inf_tot$Group %in% "11", "Group_condition"] <- "AD with infection"

ggplot(pheno_MEs0_resp_inf_tot, aes(x = as.factor(Group_condition), y = green)) +
  geom_boxplot(aes(fill = as.factor(Group_condition)))

ggplot(pheno_MEs0_resp_inf_tot, aes(x = as.factor(Group_condition), y = yellow)) +
  geom_boxplot(aes(fill = as.factor(Group_condition)))

############################ Hub genes (legacy section) ############################

labels    <- module_df_resp_inf_tot[match(rownames(residuals_lm_resp_inf_tot), module_df_resp_inf_tot$gene_id), "colors"]
hub_genes <- chooseTopHubInEachModule(t(residuals_lm_resp_inf_tot), labels, 4)

############################ Identify hub genes with MM/GS thresholds ############################

path <- "output/figures/mm/"
moduleGeneSignificanceInfection <- moduleMembership(
  "Infection",
  netwk_cell_no_fetal,
  module_df_cell_no_fetal[[3]],
  module_df_cell_no_fetal[[2]],
  residuals_cell_no_fetal,
  path,
  0.75,
  0.5
)
moduleGeneSignificanceAD <- moduleMembership(
  "AD",
  netwk_cell_no_fetal,
  module_df_cell_no_fetal[[3]],
  module_df_cell_no_fetal[[2]],
  residuals_cell_no_fetal,
  path,
  0.75,
  0.5
)

hub_genes_yellow_infection <- attr(moduleGeneSignificanceInfection, "hub_genes")$yellow
hub_genes_yellow_AD        <- attr(moduleGeneSignificanceAD,        "hub_genes")$yellow
hub_genes_yellow           <- intersect(hub_genes_yellow_AD, hub_genes_yellow_infection)

write.table(
  hub_genes_yellow,
  "output/data/hub_genes_yellow.txt",
  quote = FALSE,
  row.names = FALSE
)
