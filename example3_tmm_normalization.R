# a small example on how to perform GeTMM normalization (TMM normalization for RPK data),
# related data processing and some preliminary exploration of the normalized and non-normalized data

# load libraries
library(edgeR)
library(umap)
library(pheatmap)
library(RColorBrewer)

# load the needed data
load("example3_data.RData")

# now we already have TPM data. This time the data is summarized to KEGG orthologous group (KO) level
# we have KO groups as rows and samples as columns. 

# first, filter the existing TPM data to remove lowly expressed genes / features with little signal
# as these can likely to be just technical noise and such spurious counts and will bias the
# the following downstream analysis

# in this data, for a gene / feature not to be filtered away, it is required to have at least some expression in enough samples
# in this case, the feature is required to have a tpm value corresponding to 5 counts in at least 12 samples 
# (12 being the size of the smallest sample group in this data)

# one popular threshold you see in RNASeq literature a lot is CPM = 1 (CPM data is quite commonly used)
# however, it is hard to say how strict this filtering is in practice when thinking about the raw data,
# as it can mean anything from 0.5 to 50 counts (or more or less), depending on your library sizes
# of course it has an intuitive interpration as 1 counts per million is always 1 counts per million
# this could be adopted to TPM data as well (e.g. by thinking that a feature needs to have at least 1 TPM of 
# expression in at least n samples)

# however, here we will apply a slightly different strategy for this data
# we want to use the sample with the smallest library size here for filtering, as it will produce larger TPM values
# for similar RPK and count values to other samples, and we want to be rather strict in our filtering (typically I think
# being strict with filtering is better than lenient, it makes the remaining signal in your data to be more reliable, but
# of course this depends on the data and the project and may require some exploration of different alternatives)

# so we are filtering away such features, which based on the sample with the smallest library size, doesn't have
# at least 5 counts (or the corresponding TPM) in at least 12 samples, which is the size of the 
# smallest sample group / condition in this data

# this method is similar to how filtering low expression genes is performed for RNASeq data in the edgeR manual

# first, we will determine the sample with the smallest library size
loc_min_sample <- which.min(rpk_scaling_factors)

# then we will get the features with 5 counts in the sample with the smallest library size
locs_limit_vals <- which(count_data[,loc_min_sample]==5)

# then we will determine a TPM value as the median TPM value for these features in that sample
tpm_thres <- median(tpm_data[locs_limit_vals, loc_min_sample])

# keep only those genes having at least that much expression in 12 samples
# don't predetermine the samples where you need to have that much expression,
# as it too can bias your results. Just determine the number minimun number of samples needed

keep_genes <- rowSums(tpm_data > tpm_thres ) >= 12 
tpm_filt <- tpm_data[keep_genes,]
rpk_filt <- rpk_data[keep_genes,]

# how does this data look like, explore:
boxplot(log2(tpm_filt+1))

# compare it to the unfiltered data
boxplot(log2(tpm_data+1)) 
# mostly zeroes and very lowly expressed features

# the filtering seems to have been effective, we got rid most of the lowly expressed features
# however, by looking at the data, some extra normalization could be appropriate? 
# at least we could try out a normalization and see what it results in

# normalize using GeTMM
# this is basicly performing TMM normalization on RPK transformed data and then dividing by the rpk per million scaling factor
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6013957/

# use the edgeR package for this.
# make a DGEList object with the rpk data and the library sizes (sample total rpks) for normalization.
tpm_norm <- DGEList(rpk_filt, lib.size = c(rpk_scaling_factors)*10^6) # as rpk_scaling_factors were just total rpks / 10^6, we can get the total rpks from the scaling factors

# normalize with the TMM normalization and the default parameters
tpm_norm <- calcNormFactors(tpm_norm, method = "TMM")

# CPM (count per million) transform the normalized RPK data
# CPM transformation for the RPK data produces TPM data
tpm_norm <- cpm(tpm_norm, normalized.lib.sizes = T)

# how does this data look like, explore:
boxplot(log2(tpm_norm+1))
# does the data look different? 
# which data should be preferred?
# remember, the TMM normalization assumes that most genes / features are not differentially expressed between samples
# does this assumption hold for this data type? For KEGG orthologous gropus?

# parse the sample names of the datasets again to match the sample metadata
colnames(tpm_filt) <- paste("P",unlist(lapply(colnames(tpm_filt), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")
colnames(tpm_norm) <- paste("P",unlist(lapply(colnames(tpm_norm), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")

tpm_filt <- tpm_filt[,rownames(sample_table)]
tpm_norm <- tpm_norm[,rownames(sample_table)]

# let's do some PCA
pca_logged <- prcomp(t(log2(tpm_filt+1)))
pca_normed <- prcomp(t(log2(tpm_norm+1)))

# define some well differentiating colors
cols_diff <- c("purple","cornflowerblue","mistyrose2","darkgrey", "deeppink","brown","violet","firebrick1","chocolate3","forestgreen","darkslateblue","darkolivegreen2", "yellow")

# color according to experimental block
cols_pca <- cols_diff[as.numeric(sample_table$Block)]

plot(pca_logged$x[,1:2], col = cols_pca, pch=19, cex=2, xlab = paste("PC1 ", "(",(100*round((pca_logged$sdev[1])^2/(sum(pca_logged$sdev^2)),3)),"%)" , sep = ""),
     ylab = paste("PC2 ", "(",(100*round((pca_logged$sdev[2])^2/(sum(pca_logged$sdev^2)),3)),"%)" , sep = ""), main = "Block")

plot(pca_normed$x[,1:2], col = cols_pca, pch=19, cex=2, xlab = paste("PC1 ", "(",(100*round((pca_normed$sdev[1])^2/(sum(pca_normed$sdev^2)),3)),"%)" , sep = ""),
     ylab = paste("PC2 ", "(",(100*round((pca_normed$sdev[2])^2/(sum(pca_normed$sdev^2)),3)),"%)" , sep = ""), main = "Block")

# color according to the condition of interest
cols_pca <- rep("blue", nrow(sample_table))
cols_pca[which(sample_table$Condition=="Condition2")] <- "red"

plot(pca_logged$x[,1:2], col = cols_pca, pch=19, cex=2, xlab = paste("PC1 ", "(",(100*round((pca_logged$sdev[1])^2/(sum(pca_logged$sdev^2)),3)),"%)" , sep = ""),
     ylab = paste("PC2 ", "(",(100*round((pca_logged$sdev[2])^2/(sum(pca_logged$sdev^2)),3)),"%)" , sep = ""), main = "Condition")

plot(pca_normed$x[,1:2], col = cols_pca, pch=19, cex=2, xlab = paste("PC1 ", "(",(100*round((pca_normed$sdev[1])^2/(sum(pca_normed$sdev^2)),3)),"%)" , sep = ""),
     ylab = paste("PC2 ", "(",(100*round((pca_normed$sdev[2])^2/(sum(pca_normed$sdev^2)),3)),"%)" , sep = ""), main = "Condition")

# perhaps we can do some other ordination as well
# do UMAP
set.seed(123)
umap_logged <- umap(t(log2(tpm_filt+1)))

set.seed(123)
umap_normed <- umap(t(log2(tpm_norm+1)))

# plot
cols_pca <- cols_diff[as.numeric(sample_table$Block)]
plot(umap_logged$layout[,1], umap_logged$layout[,2], col=cols_pca, pch=19, xlab = "UMAP1", ylab =  "UMAP2", cex=2, main = "Block")
plot(umap_normed$layout[,1], umap_normed$layout[,2], col=cols_pca, pch=19, xlab = "UMAP1", ylab =  "UMAP2", cex=2, main = "Block")

cols_pca <- rep("blue", nrow(sample_table))
cols_pca[which(sample_table$Condition=="Condition2")] <- "red"
    
plot(umap_logged$layout[,1], umap_logged$layout[,2], col=cols_pca, pch=19, xlab = "UMAP1", ylab =  "UMAP2", cex=2, main = "Condition")
plot(umap_normed$layout[,1], umap_normed$layout[,2], col=cols_pca, pch=19, xlab = "UMAP1", ylab =  "UMAP2", cex=2, main = "Condition")

# look at clustering based on sample correlations
sample_cors_logged <- cor(log2(tpm_filt+1), use = "pairwise.complete.obs", method = "pearson")
sample_cors_normed <- cor(log2(tpm_norm+1), use = "pairwise.complete.obs", method = "pearson")

# hierarchical clustering with complete linkage
clust_complete_logged <- hclust(d = as.dist(1-sample_cors_logged), method = "complete")
clust_complete_normed <- hclust(d = as.dist(1-sample_cors_normed), method = "complete")

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

# define the scale and colors to be used
r <- c(min(sample_cors_logged),quantile(sample_cors_logged, 0.95))
breakslist <- seq(from=r[1], to=r[2], length.out=11)
cols_heat <- rev(brewer.pal(11,"RdBu"))

# plot
pheatmap(mat = sample_cors_logged, cluster_rows = clust_complete_logged, cluster_cols = clust_complete_logged, annotation_col = annotation_col, annotation_colors = annotation_colors, scale="none", cellwidth = 10, cellheight = 10, breaks=breakslist, color = cols_heat, fontsize_row = 8, fontsize_col = 8, filename = NA, main = "Complete linkage")

# GeTMM normalized data
r <- c(min(sample_cors_normed),quantile(sample_cors_normed, 0.95))
breakslist <- seq(from=r[1], to=r[2], length.out=11)
pheatmap(mat = sample_cors_normed, cluster_rows = clust_complete_normed, cluster_cols = clust_complete_normed, annotation_col = annotation_col, annotation_colors = annotation_colors, scale="none", cellwidth = 10, cellheight = 10, breaks=breakslist, color = cols_heat, fontsize_row = 8, fontsize_col = 8, filename = NA, main = "Complete linkage")

# cluster using euclidean distance with complete linkage and scale the data featurewise
breakslist <- seq(from=-2, to=2, length.out=11)
pheatmap(mat = log2(tpm_filt+1), cluster_rows = T, cluster_cols = T, annotation_col = annotation_col, annotation_colors = annotation_colors, scale="row", breaks=breakslist, show_rownames = F, color = cols_heat, fontsize_row = 8, fontsize_col = 8, filename = NA, main = "Complete linkage")
pheatmap(mat = log2(tpm_norm+1), cluster_rows = T, cluster_cols = T, annotation_col = annotation_col, annotation_colors = annotation_colors, scale="row", show_rownames = F, breaks=breakslist, color = cols_heat, fontsize_row = 8, fontsize_col = 8, filename = NA, main = "Complete linkage")

# is there anything that suggests which data you should use, the GeTMM normalized or non-normalized?
