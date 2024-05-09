getwd()

setwd('/gstore/scratch/u/gaos38/sc_T')
list.files()

library(scuttle)
library(scran)
library(MAST)
library(data.table)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# adata v1.2 has been exported as RData for MAST analysis
load('sce_v3_hash.RData')

counts(sce_v3) <- assay(sce_v3, "X")
assay(sce_v3)
head(colData(sce_v3))

id_cells_C0 <- which(colData(sce_v3)$hash == 'LV_2KO')

# identify cells belonging to rest of the clusters (i.e. not in Cluster 0):
id_cells_rest <- which(colData(sce_v3)$hash == 'LV_3KO')

# Check the length of the two sets
print(length(id_cells_C0))
print(length(id_cells_rest))


#df1 lv_2ko
#df2 lv_3ko
df1 <- t(data.frame(counts(sce_v3)[, id_cells_C0])) # transpose because in sce genes are rows
df2 <- t(data.frame(counts(sce_v3)[, id_cells_rest])) # transpose because in sce genes are rows

id_cells_lvwt <- which(colData(sce_v3)$hash == 'LV_WT')
df_lvwt <- t(data.frame(counts(sce_v3)[, id_cells_rest]))


#source("run_MAST.r")
#pairwise_de(df1, df2, 'Cluster_0_vs_rest', '/gstore/scratch/u/gaos38')

##############
# ORdered data frame
o_df <- rbind(df1, df2)
condition <- rep(c(1, 2), times=c(nrow(df1), nrow(df2)))

# Prepare for MAST
wellKey <- rownames(o_df)
cdata <- data.frame(cbind(wellKey=wellKey, condition=condition))
fdata <- data.frame(primerid=colnames(o_df))

# SCa data
sca <- FromMatrix( t(o_df), cdata, fdata)
cdr2 <-colSums(assay(sca)>0)
colData(sca)$cngeneson <- scale(cdr2)
colData(sca)$cond <- factor(unlist(as.list(condition)))
colData(sca)$cond <- relevel(colData(sca)$cond, 2)

# Fits (most time consuming 2 steps)
zlmCond <- zlm(~cond + cngeneson, sca)
summaryDt <- summary(zlmCond, doLRT='cond1')$datatable

# Significance table
fcHurdle <- merge(summaryDt[contrast=='cond1' & component=='H',.(primerid, `Pr(>Chisq)`)], 
                  summaryDt[contrast=='cond1' & component=='logFC', .(primerid,coef, ci.hi, ci.lo)], by='primerid') 

# FDR
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
setorder(fcHurdle, fdr)

# Export data
write.csv(as.data.frame(fcHurdle), sprintf("%s%s.csv", "/gstore/scratch/u/gaos38/sc_T/", 'MAST_lv2ko3ko'), quote=FALSE)

#####################

for (j in unique(colData(sce_v2)$phenograph_label)) {
  print(paste0('Analysing Cluster', toString(j)))
  # identify cells belonging a specific cluster
  id_cells <- which(colData(sce_v2)$phenograph_label == j)
  # identify cells belonging to rest of the clusters
  id_cells_rest <- which(colData(sce_v2)$phenograph_label != j)
  # Create two dataframes: 
  df1 <- t(data.frame(counts(sce_v2)[, id_cells])) # transpose because in sce genes are rows
  df2 <- t(data.frame(counts(sce_v2)[, id_cells_rest])) # transpose because in sce genes are rows
  
  # use this function
  source("run_MAST.r")
  
  file_name_temp = paste0('Cluster_', toString(j), '_vs_rest')
  # Syntax: pairwise_de(dataframe1, dataframe2, 'output_filename', 'output_folder')
  pairwise_de(df1, df2, file_name_temp, outbase)
}


########## clean code for volcano
df_MAST_lv2ko3ko <- read.csv("MAST_lv2ko3ko.csv")
head(df_MAST_lv2ko3ko)

data <- df_MAST_lv2ko3ko %>% select(primerid, coef, fdr)
head(data)

#rename
colnames(data)[colnames(data) == "primerid"] <- "gene"
colnames(data)[colnames(data) == "coef"] <- "log2FC"
colnames(data)[colnames(data) == "fdr"] <- "adj_p"

data$neg_log_p_value <- -log10(data$adj_p)
# change the inf to the max -logp value
max_finite <- max(data$neg_log_p_value[is.finite(data$neg_log_p_value)], na.rm = TRUE)
data$neg_log_p_value[!is.finite(data$neg_log_p_value)] <- max_finite

# characterize sig DEG, up, down
data$diffexpressed <- "NO"
data$diffexpressed[data$log2FC > log2(1.4) & data$adj_p < 0.05] <- "UP"
data$diffexpressed[data$log2FC < -log2(1.4) & data$adj_p < 0.05] <- "DOWN"

# prepare the lable for DEG
data$delabel <- NA
data$delabel[data$diffexpressed != "NO"] <- data$gene[data$diffexpressed != "NO"]

#plot
ggplot(data, aes(x=log2FC, y=neg_log_p_value, col=diffexpressed, label=delabel)) + 
  geom_point(alpha = 0.5) + 
  theme_minimal() +
  geom_text_repel(size = 3,
                  max.overlaps = Inf) +
  scale_colour_manual(values = mycolors) +
  geom_vline(xintercept=c(-log2(1.4), log2(1.4)), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

######################## batch plot for volcano
getwd()
setwd('/Users/gaos1/Desktop/scRNA_temp/MAST_result_01')

library(scuttle)
library(scran)
library(data.table)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggrepel)

########
prefix = 'MAST_LN_2KO_vs_'
files <- c('LN_WT', 'LN_RII', 'LN_PD1', 'LN_3KO', 'LV_2KO')
#
sig_log2FC = 1.4
sig_p = 0.05
mycolors <- c("blue", "black", "red")
#
for (file in files){
  filename = paste0(prefix, file, '.csv')
  print(paste0('working on ', filename))
  df = read.csv(filename)
  data <- df %>% select(primerid, coef, fdr)
  
  #rename
  colnames(data)[colnames(data) == "primerid"] <- "gene"
  colnames(data)[colnames(data) == "coef"] <- "log2FC"
  colnames(data)[colnames(data) == "fdr"] <- "adj_p"
  
  data$neg_log_p_value <- -log10(data$adj_p)
  # change the inf to the max -logp value
  max_finite <- max(data$neg_log_p_value[is.finite(data$neg_log_p_value)], na.rm = TRUE)
  data$neg_log_p_value[!is.finite(data$neg_log_p_value)] <- max_finite
  
  # characterize sig DEG, up, down
  data$diffexpressed <- "NO"
  data$diffexpressed[data$log2FC > log2(sig_log2FC) & data$adj_p < sig_p] <- "UP"
  data$diffexpressed[data$log2FC < -log2(sig_log2FC) & data$adj_p < sig_p] <- "DOWN"
  
  # prepare the lable for DEG
  data$delabel <- NA
  data$delabel[data$diffexpressed != "NO"] <- data$gene[data$diffexpressed != "NO"]
  
  #plot
  p <- ggplot(data, aes(x=log2FC, y=neg_log_p_value, col=diffexpressed, label=delabel)) + 
    geom_point(alpha = 0.5) + 
    theme_minimal() +
    geom_text_repel(size = 3,
                    max.overlaps = Inf) +
    scale_colour_manual(values = mycolors) +
    geom_vline(xintercept=c(-log2(sig_log2FC), log2(sig_log2FC)), col="red") +
    geom_hline(yintercept=-log10(sig_p), col="red") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  plot_name = paste0(prefix, file, '_DEG_volcano.pdf')
  ggsave(plot_name, plot = p, width = 10, height = 8)
}  

########
prefix = 'MAST_LV_2KO_vs_'
files <- c('LV_WT', 'LV_RII', 'LV_PD1', 'LV_3KO')
#
sig_log2FC = 1.4
sig_p = 0.05
mycolors <- c("blue", "black", "red")
#
for (file in files){
  filename = paste0(prefix, file, '.csv')
  print(paste0('working on ', filename))
  df = read.csv(filename)
  data <- df %>% select(primerid, coef, fdr)
  
  #rename
  colnames(data)[colnames(data) == "primerid"] <- "gene"
  colnames(data)[colnames(data) == "coef"] <- "log2FC"
  colnames(data)[colnames(data) == "fdr"] <- "adj_p"
  
  data$neg_log_p_value <- -log10(data$adj_p)
  # change the inf to the max -logp value
  max_finite <- max(data$neg_log_p_value[is.finite(data$neg_log_p_value)], na.rm = TRUE)
  data$neg_log_p_value[!is.finite(data$neg_log_p_value)] <- max_finite
  
  # characterize sig DEG, up, down
  data$diffexpressed <- "NO"
  data$diffexpressed[data$log2FC > log2(sig_log2FC) & data$adj_p < sig_p] <- "UP"
  data$diffexpressed[data$log2FC < -log2(sig_log2FC) & data$adj_p < sig_p] <- "DOWN"
  
  # prepare the lable for DEG
  data$delabel <- NA
  data$delabel[data$diffexpressed != "NO"] <- data$gene[data$diffexpressed != "NO"]
  
  #plot
  p <- ggplot(data, aes(x=log2FC, y=neg_log_p_value, col=diffexpressed, label=delabel)) + 
    geom_point(alpha = 0.5) + 
    theme_minimal() +
    geom_text_repel(size = 3,
                    max.overlaps = Inf) +
    scale_colour_manual(values = mycolors) +
    geom_vline(xintercept=c(-log2(sig_log2FC), log2(sig_log2FC)), col="red") +
    geom_hline(yintercept=-log10(sig_p), col="red") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  plot_name = paste0(prefix, file, '_DEG_volcano.pdf')
  ggsave(plot_name, plot = p, width = 10, height = 8)
}  


################################## ploting for clusters
prefix = 'wilcoxon_cluster_'
#
sig_log2FC = 8
sig_p = 0.05
mycolors <- c("blue", "black", "red")
#
for (i in 0:20){
  file = as.character(i)
  filename = paste0(prefix, file, '_vs_rest.csv')
  print(paste0('working on ', filename))
  df = read.csv(filename)
  
  data <- df %>% select(names, logfoldchanges, pvals_adj)
  
  #rename
  colnames(data)[colnames(data) == "names"] <- "gene"
  colnames(data)[colnames(data) == "logfoldchanges"] <- "log2FC"
  colnames(data)[colnames(data) == "pvals_adj"] <- "adj_p"
  
  data$neg_log_p_value <- -log10(data$adj_p)
  # change the inf to the max -logp value
  max_finite <- max(data$neg_log_p_value[is.finite(data$neg_log_p_value)], na.rm = TRUE)
  data$neg_log_p_value[!is.finite(data$neg_log_p_value)] <- max_finite
  
  # characterize sig DEG, up, down
  data$diffexpressed <- "NO"
  data$diffexpressed[data$log2FC > log2(sig_log2FC) & data$adj_p < sig_p] <- "UP"
  data$diffexpressed[data$log2FC < -log2(sig_log2FC) & data$adj_p < sig_p] <- "DOWN"
  
  # prepare the lable for DEG
  data$delabel <- NA
  data$delabel[data$diffexpressed != "NO"] <- data$gene[data$diffexpressed != "NO"]
  
  #plot
  p <- ggplot(data, aes(x=log2FC, y=neg_log_p_value, col=diffexpressed, label=delabel)) + 
    geom_point(alpha = 0.5) + 
    theme_minimal() +
    geom_text_repel(size = 3,
                    max.overlaps = Inf) +
    scale_colour_manual(values = mycolors) +
    geom_vline(xintercept=c(-log2(sig_log2FC), log2(sig_log2FC)), col="red") +
    geom_hline(yintercept=-log10(sig_p), col="red") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  plot_name = paste0(prefix, file, '_DEG_volcano.pdf')
  ggsave(plot_name, plot = p, width = 10, height = 8)
} 