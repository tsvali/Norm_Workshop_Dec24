# a simple example to demonstrate how batch effects can be created 
# during data processing with slight variations in the processing workflow.

# load libraries
library(RColorBrewer)
library(pheatmap)

# add installed custom packages in Puhti
.libPaths(c("/scratch/project_2009164/2_OULANKA/Tommi/normalization_workshop/RPackages", .libPaths()))
library(ROTS)

# load the needed data
load("example5_data.RData")

# both the object "tpm_batch" and "tpm_data" contain gene count data from the contigs 
# aggregated to the KEGG KO group level.

# the contigs have been generated from two coassemblies, one from Condition 1 and one from
# Condition 2, which were observed to be distinct / different based on exploration of the 
# 16S data from the samples as well as the contig data from single and multiple different coassemblies.

# in the "tpm_batch" data, only the samples used in the coassemblies have been mapped to the
# contigs with gene counts for the contigs genes generated only for those samples, aggregated
# to KO level and then the KO level data from the samples from the different coassemblies (Condition 1 / 2)
# have been combined. 

# whereas in the "tpm_data", the contigs from the different coassemblies have been concantenated
# into one contig database to which all the samples have been mapped to (so all samples are mapped to genes
# from both coassemblies) and then aggregated to KO level.

# make sure everything is in the right order, sample names are already parsed.
tpm_batch <- tpm_batch[,rownames(sample_table)]
tpm_data <- tpm_data[,rownames(sample_table)]

# let's plot the datas in a few different days.
# first, let's look at the total TPM in the samples in both datas

# color by condition / coassembly
cols_bar <- rep("blue", nrow(sample_table))
cols_bar[which(sample_table$Condition=="Condition2")] <- "red"

barplot(colSums(tpm_batch), col=cols_bar, las=2, ylab = "TPM")
barplot(colSums(tpm_data), col=cols_bar, las=2, ylab = "TPM")
# is there a difference?

# look at sample clustering 

# PCA
pca_batch <- prcomp(t(log2(tpm_batch+1)))
pca_data <- prcomp(t(log2(tpm_data+1)))

# color by condition
cols_pca <- rep("blue", nrow(sample_table))
cols_pca[which(sample_table$Condition=="Condition2")] <- "red"

plot(pca_batch$x[,1:2], col = cols_pca, pch=19, cex=2, xlab = paste("PC1 ", "(",(100*round((pca_batch$sdev[1])^2/(sum(pca_batch$sdev^2)),3)),"%)" , sep = ""),
     ylab = paste("PC2 ", "(",(100*round((pca_batch$sdev[2])^2/(sum(pca_batch$sdev^2)),3)),"%)" , sep = ""), main = "Condition")

plot(pca_data$x[,1:2], col = cols_pca, pch=19, cex=2, xlab = paste("PC1 ", "(",(100*round((pca_data$sdev[1])^2/(sum(pca_data$sdev^2)),3)),"%)" , sep = ""),
     ylab = paste("PC2 ", "(",(100*round((pca_data$sdev[2])^2/(sum(pca_data$sdev^2)),3)),"%)" , sep = ""), main = "Condition")


# sample clustering based on correlation distance
sample_cors_batch <- cor(log2(tpm_batch+1), use = "pairwise.complete.obs", method = "pearson")
sample_cors_data <- cor(log2(tpm_data+1), use = "pairwise.complete.obs", method = "pearson")

# hierarchical clustering with complete linkage
clust_complete_batch <- hclust(d = as.dist(1-sample_cors_batch), method = "complete")
clust_complete_data <- hclust(d = as.dist(1-sample_cors_data), method = "complete")

# plot with heatmap
# annotations for plotting
annotation_col <- sample_table[,c(1,2)]
for(i in 1:ncol(annotation_col)){annotation_col[,i] <- as.character(annotation_col[,i])}

# colors for plotting
annotation_colors <- list(
  Block = c("1"= "purple", "2"="cornflowerblue", "3"="mistyrose2", "4"="darkgrey", "5"="deeppink", "6"="brown",
            "7"="violet", "8"="firebrick1", "9"="chocolate3", "10"="forestgreen", "11"="darkslateblue", "12"="darkolivegreen2"),
  Condition = c("Condition1"="blue", "Condition2"="red")
)

r <- c(min(sample_cors_batch),quantile(sample_cors_batch, 0.95))
breakslist <- seq(from=r[1], to=r[2], length.out=11)
cols_heat <- rev(brewer.pal(11,"RdBu"))
pheatmap(mat = sample_cors_batch, cluster_rows = clust_complete_batch, cluster_cols = clust_complete_batch, annotation_col = annotation_col, annotation_colors = annotation_colors, scale="none", cellwidth = 10, cellheight = 10, breaks=breakslist, color = cols_heat, fontsize_row = 8, fontsize_col = 8, filename = NA, main = "Complete linkage")

r <- c(min(sample_cors_data),quantile(sample_cors_data, 0.95))
breakslist <- seq(from=r[1], to=r[2], length.out=11)
cols_heat <- rev(brewer.pal(11,"RdBu"))
pheatmap(mat = sample_cors_data, cluster_rows = clust_complete_data, cluster_cols = clust_complete_data, annotation_col = annotation_col, annotation_colors = annotation_colors, scale="none", cellwidth = 10, cellheight = 10, breaks=breakslist, color = cols_heat, fontsize_row = 8, fontsize_col = 8, filename = NA, main = "Complete linkage")

# how does it look?

# let's also do some differential expression and look at thre results.
# batch effect data
ROTS_batch <- ROTS(data = log2(tpm_batch+1), groups = c(rep(1,18), rep(0,18)), B = 1000, K = nrow(tpm_batch)/4, paired = F, progress = T, log = T, seed = 123)
ROTS_batch_df <- data.frame(d = ROTS_batch$d, logfc = ROTS_batch$logfc, pvalue = ROTS_batch$pvalue, fdr = ROTS_batch$FDR)

# order according to p-value
ROTS_batch_df <- ROTS_batch_df[order(as.numeric(ROTS_batch_df$pvalue)),]

# concatenated data
ROTS_data <- ROTS(data = log2(tpm_data+1), groups = c(rep(1,18), rep(0,18)), B = 1000, K = nrow(tpm_data)/4, paired = F, progress = T, log = T, seed = 123)
ROTS_data_df <- data.frame(d = ROTS_data$d, logfc = ROTS_data$logfc, pvalue = ROTS_data$pvalue, fdr = ROTS_data$FDR)

# order according to p-value
ROTS_data_df <- ROTS_data_df[order(as.numeric(ROTS_data_df$pvalue)),]

# how much do we have differential expression in each data set? How many KO groups
# show differential expression between Condition 1 and Condition according to ROTS?
length(which(ROTS_batch_df$fdr<=0.05))
length(which(ROTS_data_df$fdr<=0.05))

# how much is that in proportions, what is the percentage of KO groups showing
# differential expression between Condition 1 and Condition 2?
round((length(which(ROTS_batch_df$fdr<=0.05)) / nrow(ROTS_batch_df))*100)
round((length(which(ROTS_data_df$fdr<=0.05)) / nrow(ROTS_data_df))*100)

# which one seems more likely? Which data seems better? 
# basicly the datasets are resulting from the same coassemblies but are just 
# postprocessed slightly differently. 
# can you see how a huge rule also the data processing steps can play to data quality?

