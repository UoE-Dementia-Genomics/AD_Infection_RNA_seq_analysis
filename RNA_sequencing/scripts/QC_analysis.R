########################
# @author: Giulia Pegoraro
#
# RNA-seq QC and preprocessing pipeline – STEP 1
#
# This script:
#  - Loads gene-level quantification files and harmonises sample IDs
#  - Imports and merges phenotype, RIN, TIN, cell-type proportion and batch metadata
#  - Filters lowly expressed genes and computes TMM-normalised log2-CPM counts
#  - Performs sample-level QC (hierarchical clustering + Mahalanobis distance)
#    to identify and remove outlier samples
#  - Runs PCA on normalised counts and correlates PCs with technical/biological covariates
#  - Computes correlation structure among phenotype covariates
#  - Iteratively removes batch effects (library size, sequencing plate, TIN) with ComBat
#  - Exports cleaned expression and phenotype matrices for downstream analyses
########################

##########################################################
######################## LIBRARIES #######################
##########################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
if (!requireNamespace("readr", quietly = TRUE))
  install.packages("readr")
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")
if (!requireNamespace("edgeR", quietly = TRUE))
  BiocManager::install("edgeR")
if (!requireNamespace("Hmisc", quietly = TRUE))
  install.packages("Hmisc")
if (!requireNamespace("corrplot", quietly = TRUE))
  install.packages("corrplot")
if (!requireNamespace("fastDummies", quietly = TRUE))
  install.packages("fastDummies")
if (!requireNamespace("polycor", quietly = TRUE))
  install.packages("polycor")
if (!requireNamespace("parallel", quietly = TRUE))
  install.packages("parallel")
if (!requireNamespace("data.table", quietly = TRUE))
  install.packages("data.table")
if (!requireNamespace("sva", quietly = TRUE))
  BiocManager::install("sva")
if (!requireNamespace("Haplin", quietly = TRUE))
  install.packages("Haplin")
if (!requireNamespace("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")
if (!requireNamespace("WGCNA", quietly = TRUE))
  BiocManager::install("WGCNA")
if (!requireNamespace("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("tidyr", quietly = TRUE))
  install.packages("tidyr")
if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")
if (!requireNamespace("limma", quietly = TRUE))
  BiocManager::install("limma")
if (!requireNamespace("dendextend", quietly = TRUE))
  install.packages("dendextend")

library(readr)
library(dplyr)
library(limma)
library(edgeR)
library(Hmisc)
library(corrplot)
library(fastDummies)
library(sva)
library(polycor)
library(parallel)
library(data.table)
library(biomaRt)
library(WGCNA)
library(tidyr)
library(purrr)
library(readxl)
library(DESeq2)
library(ggplot2)
library(dendextend)
library(ggfortify)

##########################################################
######################## FUNCTIONS #######################
##########################################################

############################################
# Function: pca.correlation
#
# Description:
#   Compute correlations between selected principal components and
#   phenotype variables, and save a correlation heatmap to file.
#
# Arguments:
#   pca_score : numeric matrix/data frame of PCA scores (samples x PCs).
#   pheno_df  : data frame of phenotype variables (samples in rows).
#   n_PC      : integer, number of principal components to include.
#   file_name : character, suffix for the output file name.
#
# Returns:
#   None (PNG plot is written to disk).
############################################
pca.correlation <- function(pca_score, pheno_df, n_PC, file_name){
  
  sub_pheno <- pheno_df
  princcomps <- as.data.frame(pca_score[,1:n_PC])
  princcomps_PCA <- merge(sub_pheno, princcomps,  by ="row.names")
  princcomps_PCA <- princcomps_PCA[,-1]
  
  # Correlation matrix for p-values
  cor_p <- cor.mtest(princcomps_PCA, conf.level=0.95)
  p.val_df <- cor_p$p
  p.val_df <- p.val_df[1:ncol(sub_pheno),(ncol(sub_pheno)+1):ncol(p.val_df)]
  
  # Raw correlation matrix (PCs vs phenotype variables)
  cor_PCA <- cor(princcomps_PCA[,(ncol(sub_pheno)+1):ncol(princcomps_PCA)], 
                 princcomps_PCA[,1:ncol(sub_pheno)],
                 use = "complete.obs")
  colnames(cor_PCA) <- colnames(sub_pheno)
  
  png(paste0("C:\\Users\\gp487\\OneDrive - University of Exeter\\AD_infection\\Analyses\\RNA_sequencing\\output\\figures\\correlation",file_name),
      width=2400, height=1500, res=300)
  FullPCA_Corr <- corrplot(
    cor_PCA,
    p.mat = t(p.val_df),
    addrect = 2,
    method = "color",
    insig = "label_sig",
    sig.level = c(0.001, 0.01, 0.05), 
    tl.cex = 0.6,
    cl.cex = 0.4,
    cl.offset = 0.2,
    cl.align.text	= "l",
    pch.cex = 0.6,
    pch.col = "grey20"
  )
  dev.off()
}

############################################
# Function: variable.correlation
#
# Description:
#   Compute correlations between all principal components and
#   phenotype variables, and save a correlation heatmap to file.
#
# Arguments:
#   pca_score : numeric matrix/data frame of PCA scores (samples x PCs).
#   pheno_df  : data frame of phenotype variables (samples in rows).
#   file_name : character, output file name.
#
# Returns:
#   None (PNG plot is written to disk).
############################################
variable.correlation <- function(pca_score, pheno_df, file_name){
  
  sub_pheno <- pheno_df
  princcomps_PCA <- merge(sub_pheno, pca_score,  by ="row.names")
  princcomps_PCA <- princcomps_PCA[,-1]
  
  # Correlation matrix for p-values
  cor_p <- cor.mtest(princcomps_PCA, conf.level=0.95)
  p.val_df <- cor_p$p
  p.val_df <- p.val_df[1:ncol(sub_pheno),(ncol(sub_pheno)+1):ncol(p.val_df)]
  
  # Raw correlation matrix (PCs vs phenotype variables)
  cor_PCA <- cor(princcomps_PCA[,(ncol(sub_pheno)+1):ncol(princcomps_PCA)], 
                 princcomps_PCA[,1:ncol(sub_pheno)],
                 use = "complete.obs")
  colnames(cor_PCA) <- colnames(sub_pheno)
  
  png(paste0("C:\\Users\\gp487\\OneDrive - University of Exeter\\Documents\\PhD_project\\plots\\correlation\\",file_name),
      width=2400, height=1500, res=300)
  FullPCA_Corr <- corrplot(
    cor_PCA,
    p.mat = t(p.val_df),
    addrect = 2,
    method = "color",
    insig = "label_sig",
    sig.level = c(0.001, 0.01, 0.05), 
    tl.cex = 0.6,
    cl.cex = 0.4,
    cl.offset = 0.2,
    cl.align.text	= "l",
    pch.cex = 0.6,
    pch.col = "grey20"
  )
  dev.off()
}

############################################
# Function: pca.plot
#
# Description:
#   Plot the first two principal components, coloured by a specified
#   sample variable, and save the plot to file.
#
# Arguments:
#   pca_res     : prcomp object containing PCA results.
#   sample_info : data frame with sample metadata, including BBNId.
#   title       : character, title of the plot and prefix for the output file.
#   variable    : character, aesthetic mapped to point colour (column in sample_info).
#
# Returns:
#   None (PNG plot is written to disk).
############################################
pca.plot <- function(pca_res, sample_info, title, variable){
  setwd("C:\\Users\\gp487\\OneDrive - University of Exeter\\AD_infection\\Analyses\\RNA_sequencing\\output\\figures\\pca_plot")
  pc_sum <- summary(pca_res)$importance
  PC1_var <- signif(as.numeric(pc_sum["Proportion of Variance","PC1"])*100,2)
  PC2_var <- signif(as.numeric(pc_sum["Proportion of Variance","PC2"])*100,2)
  pca_res$x %>% 
    as_tibble(rownames = "BBNId") %>% 
    full_join(sample_info, by = "BBNId") %>% 
    ggplot(aes(x = PC1, y = PC2, colour = variable)) +
    geom_point() + 
    geom_text(
      label=row.names(pca_res$x), 
      nudge_x = 0.15, nudge_y = 0.15, 
      check_overlap = F
    )
  xlab(paste0("PC1(",PC1_var,"%)")) +
    ylab(paste0("PC2(",PC2_var,"%)")) +
    labs(colour = title)
  ggsave(file = paste0(title,"_pca_plot",".png"))
}

############################################
# Function: mahalanobis.outlier
#
# Description:
#   Identify multivariate outliers based on Mahalanobis distance computed
#   on either PCA coordinates or t-SNE embedding, and generate diagnostic plots.
#
# Arguments:
#   Data        : numeric matrix (features x samples).
#   method      : character, "pca" or "tsne" (dimension reduction method).
#   plot.title  : character, main title for all diagnostic plots (default NA).
#   tsne.seed   : numeric, random seed for t-SNE reproducibility (optional).
#   pca.scale   : logical, whether to scale variables before PCA (default TRUE).
#   pca.center  : logical, whether to center variables before PCA (default TRUE).
#
# Returns:
#   list with:
#     $Data.2D   : data frame with 2D coordinates, distances, p-values, and outlier flag.
#     $Plot.2D   : ggplot object of 2D embedding with outliers highlighted.
#     $Plot.Dist : ggplot object of MD vs chi-square probability.
#     $Plot.QQ   : ggplot object of MD vs chi-square quantiles (QQ-plot).
############################################
mahalanobis.outlier <- function(Data , method = "pca", plot.title=NA , tsne.seed = NA, pca.scale=T , pca.center=T){
  
  suppressMessages(library(car))
  suppressMessages(library(dplyr))
  suppressMessages(library(ggplot2))
  
  if(!is.na(tsne.seed)){
    set.seed(seed = tsne.seed) 
  }
  method = match.arg(arg = method , choices = c("pca" , "tsne") , several.ok = F)
  
  if(method == "tsne"){
    suppressMessages(library(Rtsne))
    tsne_out <- Rtsne(t(Data),dims = 2,)
    Data.2D <- data.frame(D1 = tsne_out$Y[,1], 
                          D2 = tsne_out$Y[,2])
    rownames(Data.2D) <- colnames(Data)
  }
  if(method == "pca"){
    pc <- prcomp(t(Data),scale. = pca.scale,center = pca.center , rank. =2)
    pc.importance <- summary(pc)$importance[2,]
    Data.2D <- as.data.frame(pc$x)
  }
  
  center_ <- colMeans(Data.2D)
  cov_ <- cov(Data.2D)
  
  # Calculate squared Mahalanobis distance
  Data.2D$mdist <- mahalanobis(
    x = Data.2D,
    center = center_,
    cov = cov_
  )
  
  cutoff <- qchisq(p = 0.99, df = 2)
  R <- sqrt(cutoff)
  
  ellipse_ <- car::ellipse(
    center = center_[1:2],
    shape = cov_[1:2,1:2],
    radius = R,
    segments = 150,
    draw = FALSE
  )
  ellipse_ <- as.data.frame(ellipse_)
  colnames(ellipse_) <- colnames(Data.2D)[1:2]
  
  Data.2D$pchisq <- pchisq(Data.2D$mdist, df = 2, lower.tail = FALSE)
  
  Data.2D <- Data.2D %>%
    mutate(Outlier = ifelse(mdist > cutoff, "Yes", "No"))
  
  if(method == "pca"){
    p1 <- ggplot(Data.2D, aes(x = PC1 , y = PC2, color = Outlier))+
      xlab(paste0("PC1 ",round(pc.importance[1],digits = 2)*100,"%"))+
      ylab(paste0("PC2 ",round(pc.importance[2],digits = 2)*100,"%"))
    plot.subtitle = "Dimensionality reduction method: PCA"
  }
  if(method == "tsne"){
    p1 <- ggplot(Data.2D, aes(x = D1 , y = D2, color = Outlier))
    plot.subtitle = "Dimensionality reduction method: tSNE"
  }
  if(is.na(plot.title)){
    plot.title = ""
  }
  p1 <- p1 +
    geom_point(size = 3) +
    geom_point(aes(center_[1], center_[2]) , size = 5 , color = "blue") +
    geom_polygon(data = ellipse_, fill = "white", color = "black", alpha = 0.3) +
    scale_color_manual(values = c("gray44", "red")) +
    labs(title = plot.title, subtitle = paste0("Outliers in 2D Plot, ",plot.subtitle)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0))
  
  p2 <- ggplot() +
    geom_point(data = Data.2D , aes(x = mdist , y = pchisq , color = Outlier)) +
    theme_bw() +
    scale_color_manual(values = c("black", "red"))+
    geom_vline(xintercept = cutoff , color="red")+
    ylab("Chi-Square probability")+
    xlab("Square Mahalanobis distance")+
    labs(title = plot.title, subtitle = paste0("Outliers in Chi-Square Plot, ",plot.subtitle))
  
  Data.2D$qchisq = qchisq(ppoints(length(Data.2D$mdist)), df = 2)
  p3 <- ggplot() +
    geom_point(data = Data.2D, aes(x = sort(qchisq), y = sort(mdist))) +
    theme_bw() +
    geom_abline(aes(slope = 1, intercept = 0),color="red")+
    xlab("Chi-Square quantiles")+
    ylab("Square Mahalanobis distance quantiles")+
    labs(title = plot.title, subtitle = paste0("QQ Plot, ",plot.subtitle))
  
  return(list(Data.2D = Data.2D , Plot.2D = p1 , Plot.Dist = p2 , Plot.QQ = p3))
}

##########################################################
########################## MAIN ##########################
##########################################################

############################ Input expression data ############################
# Import pre-processed RNA-seq quantification files

dir_names <- list.files("input\\data\\Expression_data\\rimmed")
files <- file.path("input","data","Expression_data", "rimmed",dir_names,"quant.genes.sf")
if(any(file.exists(files)==F)) errorCondition("One of the files doesn't exist")

############################ Input metadata ############################
# Import phenotype and technical covariate data

pheno <- read.csv("input\\data\\Metadata\\pheno_data_complete.csv", row.names = 1)
RIN_score_df <- read_excel("input\\data\\Metadata\\infection RNA RIN scores.xlsx")
cell_prop_df <- read.csv("input\\data\\Metadata\\Pheno_CellProp.csv")
Nature_cell_prop <- read.csv("input\\data\\Metadata\\Nature_cibersort_results.csv")
batch_info <- as.data.frame(read_excel("input\\data\\Metadata\\infection extraction batch info.xlsx"))

############################ Import quant.sf counts ############################
separate_counts <- lapply(files, read.table)
names(separate_counts) <- dir_names

count_df=as.data.frame(separate_counts[[1]][,1])

for (i in 1:length(separate_counts)) {
  count_df=cbind(count_df,separate_counts[[i]][,5])
}

count_df <- count_df[-1,]

############################ Harmonise sample IDs ############################
ids <- gsub("10570_", "", dir_names)
ids <- gsub("_", "", ids)
ids <- gsub("S.*", "",ids)

RIN_score_df$ID <- gsub("_","",RIN_score_df$ID)
cell_prop_df$BBNId <- gsub("\\.|_|-","",cell_prop_df$BBNId)
batch_info$BBN_ID <- gsub("\\.|_|-","",batch_info$BBN_ID)

############################ Build count matrix ############################
colnames(count_df) <- c("genes_id",ids)
rownames(count_df) <- count_df[,1]

# Convert to numeric matrix
count <- apply(count_df[,-1],2,as.numeric)
rownames(count) <- count_df$genes_id

############################ TIN scores ############################
# Add transcript integrity (TIN) scores to phenotype data

tin_files <- list.files("input\\data\\Expression_data\\tin_txt", pattern = "*.txt", full.names = F)
tin_files_path <- file.path("input","data","Expression_data", "tin_txt",tin_files)
tin_df_list <- lapply(tin_files_path, read_table)

# Extract TIN from each file
tin_df<- do.call(rbind,lapply(1:length(tin_df_list), function(x){
  return(tin_df_list[[x]][,c(1,3)])
}))

# Harmonise TIN IDs
ids_tin <- gsub("Aligned.sortedByCoord.out.summary.txt","",tin_df$Bam_file)
ids_tin <- gsub("0570_", "", ids_tin)
ids_tin <- gsub("_", "", ids_tin)
ids_tin <- gsub("S.*", "",ids_tin)

tin_df$Bam_file <- ids_tin
tin_df <- as.data.frame(tin_df)
write.csv(tin_df, "input\\data\\Expression_data\\tin_df.csv")

############################ Subset to samples in metadata ############################
count_df_meta <- count[,which(colnames(count) %in% pheno$BBNId)]

############################ Normalisation (CPM + log2) ############################

keep = edgeR::filterByExpr(count_df_meta,min.count = 10, min.prop = 0.8)
filt_count_df_meta <- count_df_meta[keep,]
write.csv(filt_count_df_meta, "input\\data\\Expression_data\\filt_count_df_meta.csv")

dge <- DGEList(filt_count_df_meta)
dge <- calcNormFactors(dge,method =c("TMM"))
log2.cpm <- edgeR::cpm(dge, log=TRUE)

if (!all(is.finite(log2.cpm))){
  geterrmessage("Infinite values in the normalised data")
}

# Plot density distribution of normalised counts
png("density_plot.png")
plot(density(log2.cpm))
dev.off()

log2.cpm <- log2.cpm[,pheno$BBNId]
write.csv(log2.cpm, "input\\data\\Expression_data\\RNA_normalised_counts.csv")

if(!identical(colnames(log2.cpm), pheno$BBNId)){
  print("Ids label not identical")
  pheno <- pheno[match(colnames(log2.cpm),pheno$BBNId),]
}

rownames(pheno) <- pheno$BBNId
write.csv(pheno, "input\\data\\Expression_data\\pheanotype_infection.csv")

############################ Merge phenotype and technical covariates ############################

# library size (sequencing depth) from DGEList
lib_df <- dge[["samples"]] # library size is equal to sequencing depth (sum of counts per sample)
lib_df$ID <- row.names(lib_df)

pheno_compl <- merge(pheno, RIN_score_df, by.x = "BBNId", by.y ="ID")
pheno_compl <- merge(pheno_compl, lib_df[,which(colnames(lib_df) %in% c("lib.size","ID"))], by.x = "BBNId", by.y = "ID" )
#pheno_compl <- merge(pheno_compl, cell_prop_df[,which(colnames(cell_prop_df) %in% c("BBNId","DoubleN","NeuNP","Sox10P"))], by = "BBNId", all.x = T)
pheno_compl <- merge(pheno_compl, tin_df, by.x = "BBNId", by.y = "Bam_file", all.x = T)
pheno_compl <- merge(pheno_compl, Nature_cell_prop, by.x = "BBNId", by.y = "X")
pheno_compl <- merge(pheno_compl, batch_info, by.x = "BBNId", by.y = "BBN_ID")

pheno_compl <- pheno_compl[, -which(colnames(pheno_compl) %in% 
                                      c("Immunised","Seq Well", "Basename", "X", "Sentrix_ID","Position"))]
pheno_compl <- pheno_compl[match(colnames(log2.cpm),pheno_compl$BBNId),]
if(!identical(colnames(log2.cpm), pheno_compl$BBNId)){errorCondition("Ids label not identical")}
write.csv(pheno_compl, "input\\data\\Metadata\\pheno_with_batch_data.csv")

############################ Hierarchical clustering (sample QC) ############################
sampleTree = hclust(dist(t(log2.cpm)), method = "average")

sizeGrWindow(12,9)
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 200, col = "red")

sampleTree_cut = hclust(dist(t(log2.cpm)), method = "average")
pheno_compl <- pheno_compl[which(pheno_compl$BBNId %in% colnames(log2.cpm)),]

dend <- as.dendrogram(sampleTree)

par(cex = 0.5);
par(mar=c(15,3,1,1))
colors_AD <- ifelse(pheno_compl$AD==0, "white", "darkblue")
colors_Infection <- ifelse(pheno_compl$Infection==0, "white", "darkblue")
colors_Gender <- ifelse(pheno_compl$Gender=="M", "lightblue", "pink")
colors_batch <- colorspace::rainbow_hcl(length(unique(pheno_compl$Batch)), l = 70)
colors_group <- colorspace::heat_hcl(length(unique(pheno_compl$Group)))
colors_braak <- colorspace::diverge_hcl(length(unique(pheno_compl$Braak)))

plot(dend, main = "Sample clustering", sub="", xlab="", cex.lab = 2,
     cex.axis = 1.5, cex.main = 3)
colored_bars(cbind(colors_braak,colors_group,colors_batch,colors_Gender,colors_Infection,colors_AD), dend = dend, 
             rowLabels = c("Braak","Group","Batch","Gender","Infection","AD"), text_shift = -0.3,
             cex.rowLabels = 2, y_shift = -3)

if(!identical(colnames(log2.cpm), pheno_compl$BBNId)){
  print("Ids label not identical")
  pheno_compl <- pheno_compl[match(colnames(log2.cpm),pheno_compl$BBNId),]
}

############################ Outlier removal (Mahalanobis distance) ############################

# Mahalanobis distance-based outlier detection
rownames(pheno_compl) <- pheno_compl$BBNId
outliers_results <- mahalanobis.outlier(log2.cpm)
outliers <- rownames(outliers_results$Data.2D[outliers_results$Data.2D$Outlier == "Yes", ])

pheno_compl <- pheno_compl[-(which(pheno_compl$BBNId %in% outliers)),]
log2.cpm <- log2.cpm[,which(colnames(log2.cpm) %in% pheno_compl$BBNId)]

# PCA and first 2 PCs coloured by infection status
pca_res <- prcomp(t(log2.cpm))
pheno_compl$Infection <- as.factor(pheno_compl$Infection)
autoplot(pca_res, data = pheno_compl, colour = ("Infection"), label = T, label.size = 3)

write.csv(pheno_compl,"input\\data\\Metadata\\Filtered_pheno_RNA.csv")

############################ Prepare phenotype matrix for correlations ############################

# Remove unwanted columns
pheno_cor <- pheno_compl[, -which(colnames(pheno_compl) %in% c("BBNId", "Braak_Tangle",  "Group"))]
# Convert gender to 0/1
pheno_cor$Gender <- ifelse(pheno_cor$Gender == "M", 1,0)
# Create dummy variables
pheno_cor <- dummy_cols(pheno_cor,select_columns=c("Brain.Bank", "Batch"))
# Remove redundant columns
pheno_cor <- pheno_cor[,-which(colnames(pheno_cor) %in% c("Brain.Bank", "Brain.Bank_Southampton", "Pos") )]
# Convert to numeric
pheno_cor <- as.data.frame(apply(pheno_cor, 2, as.numeric))
#pheno_cor <- pheno_cor[,-grep("Fetal", colnames(pheno_cor), ignore.case = T)]
row.names(pheno_cor) <- pheno_compl$BBNId
# Check labels between expression matrix and pheno data
if(!identical(colnames(log2.cpm), row.names(pheno_cor))){print("Ids label not identical")}

############################ Correlation among phenotype variables ############################

var_cor <- cor(pheno_cor, use = "complete.obs")
pheno_p <- cor.mtest(pheno_cor, conf.level=0.95)
png("plots\\correlation\\Pheno_correlation_matrix.png", height = 2000, width = 2000)
corrplot(var_cor, p.mat = t(pheno_p$p), type = "upper",  method = "color",
         insig = "label_sig", sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 2.0,tl.cex = 2.0 )
dev.off()

############################ PCA on normalised counts ############################

pca_res <- prcomp(t(log2.cpm))

# Proportion of variance explained
prop_var <- ((pca_res$sdev)^2)/sum((pca_res$sdev)^2)

# Scree plot
plot((prop_var), xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")

png("output\\figures\\PCA\\Cumulative_prop_var.png")
plot(cumsum(prop_var), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b")
dev.off()

############################ PCA–phenotype correlations ############################

pca.correlation(pca_res$x, pheno_cor, 10, "PCA_correlations_RNA_seq_AD_inf.png")

############################ Batch-effect removal (ComBat) ############################
# Convert numeric covariates to categorical quantiles for ComBat

TIN_int <- quantile(pheno_compl$`TIN(median)`)

pheno_compl[which(pheno_compl$`TIN(median)` >= TIN_int[[1]] & pheno_compl$`TIN(median)` < TIN_int[[2]]), "TIN_quant"] <- 1
pheno_compl[which(pheno_compl$`TIN(median)` >= TIN_int[[2]] & pheno_compl$`TIN(median)` < TIN_int[[3]]), "TIN_quant"] <- 2
pheno_compl[which(pheno_compl$`TIN(median)` >= TIN_int[[3]] & pheno_compl$`TIN(median)` < TIN_int[[4]]), "TIN_quant"] <- 3
pheno_compl[which(pheno_compl$`TIN(median)` >= TIN_int[[4]] & pheno_compl$`TIN(median)` <= TIN_int[[5]]), "TIN_quant"] <- 4

pheno_compl$TIN_quant <- as.factor(pheno_compl$TIN_quant)

lib.size_int <- quantile(pheno_compl$lib.size)

pheno_compl[which(pheno_compl$lib.size >= lib.size_int[[1]] & pheno_compl$lib.size < lib.size_int[[2]]), "lib.size_quant"] <- 1
pheno_compl[which(pheno_compl$lib.size >= lib.size_int[[2]] & pheno_compl$lib.size < lib.size_int[[3]]), "lib.size_quant"] <- 2
pheno_compl[which(pheno_compl$lib.size >= lib.size_int[[3]] & pheno_compl$lib.size < lib.size_int[[4]]), "lib.size_quant"] <- 3
pheno_compl[which(pheno_compl$lib.size >= lib.size_int[[4]] & pheno_compl$lib.size <= lib.size_int[[5]]), "lib.size_quant"] <- 4

pheno_compl$lib.size_quant <- as.factor(pheno_compl$lib.size_quant)

# Remove library size effect
log2.cpm_whitout_lib.size <- ComBat(dat = log2.cpm, batch = pheno_compl$lib.size_quant)
if(!identical(colnames(log2.cpm_whitout_lib.size), row.names(pheno_cor))){print("Ids label not identical")}

pca_without_lib.size <- prcomp(t(log2.cpm_whitout_lib.size))
pca.correlation(pca_without_lib.size$x, pheno_cor, 10, "PCA_correlations_no_lib.size.png")

# Remove sequencing plate effect
pheno_compl$`Seq Plate` <- as.factor(pheno_compl$`Seq Plate`)
log2.cpm_whitout_seq.plate <- ComBat(dat = log2.cpm_whitout_lib.size, batch = pheno_compl$`Seq Plate`)
if(!identical(colnames(log2.cpm_whitout_seq.plate), row.names(pheno_cor))){print("Ids label not identical")}

pca_without_seq.plate <- prcomp(t(log2.cpm_whitout_seq.plate))
pca.correlation(pca_without_seq.plate$x, pheno_cor, 10, "PCA_correlations_no_log2.cpm_seq.plate.png")

# Remove TIN effect
log2.cpm_whitout_be <- ComBat(dat = log2.cpm_whitout_seq.plate, batch = pheno_compl$TIN_quant)
if(!identical(colnames(log2.cpm_whitout_be), row.names(pheno_cor))){print("Ids label not identical")}

pca_without_TIN <- prcomp(t(log2.cpm_whitout_be))
pca.correlation(pca_without_TIN$x, pheno_cor, 10, "PCA_correlations_no_be_inf.png")

# Export final batch-corrected expression matrix
write.csv(log2.cpm_whitout_be, "input\\data\\Expression_data\\log2.cpm_whitout_be.csv")
