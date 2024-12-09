# a small example on how to perform total sum and CLR normalizations for 16S rRNA count data

# load libraries
library(phyloseq)
library(microbiome)

# add installed custom packages in Puhti
# you have to modify the path in the below statement to the folder where you installed the normalization workshop
# files in Puhti
.libPaths(c("/scratch/project_2009164/2_OULANKA/Tommi/normalization_workshop/RPackages", .libPaths()))
library(ROTS)

# load the needed data
load("example4_data.RData")

# use phyloseq and microbiome pacakges to process the data for this example
# there are other options as well - e.g. the miaverse enviroment
# https://microbiome.github.io/

# make a phyloseq object
taxa <- tax_table(as.matrix(tax_table))
otus <- otu_table(as.matrix(otu_table), taxa_are_rows = TRUE)
samples <- sample_data(sample_table)

psq <- phyloseq(taxa, otus, samples)
print(psq)

# the analysis could be performed at various taxonomical levels
# in this example, we will process order level data

# aggregate taxa to the desired rank level
psq_agg <- aggregate_taxa(x = psq, level = "Order")
print(psq_agg)

# convert to relative abundances, total sum scale (TSS)
psq_agg_rel <- microbiome::transform(x = psq_agg, transform = "compositional")

# aggregate rare taxa
psq_agg_rel_filt <- aggregate_rare(psq_agg_rel, level="Order", detection = 0.1/100, prevalence = 2/36)
print(psq_agg_rel_filt)
# so we want the taxa kept in the data to have at least some expression/counts in at least two samples
# not a very strict requirement - could be stricter

# use this data
psq_rel <- psq_agg_rel_filt

# transform otu to clr - results in the same clr transformed data from the relative and count data (from otherwise similar objects) - so relative can be used here
psq_clr <- microbiome::transform(x = psq_rel, transform = "clr")

# get th TSS normalized OTU data
data_rel <- data.frame(otu_table(psq_rel), stringsAsFactors = F, check.names = F)
data_rel <- data_rel * 100 # just to scale up.

# get the CLR normalized OTU data
data_clr <- data.frame(otu_table(psq_clr), stringsAsFactors = F, check.names = F)

# get the unique OTU names for the order level
tax_order <- data.frame(tax_table(psq_rel), stringsAsFactors = F, check.names = F)
rownames(data_rel) <- rownames(data_clr) <- tax_order$unique

# explore both tables, "data_rel" and "data_clr"
# is there an intuitive real world connection to the values?
colSums(data_rel)

# look at the values in the tables
boxplot(log2(data_rel+0.001), ylab="Log2 Proportions", main="Proportions", xlab="", las=2)
boxplot(data_clr, ylab="CLR", main="CLR", xlab="", las=2)
# you can see the pseudocount added to both datasest as the lowest 
# values (micorbiome packages CLR transformation also adds a pseudocount before transformation)

# is the order of taxons within a sample preserved by the CLR transformation?
all(rank(data_rel[,1]) == rank(data_clr[,1]))

# explore these datasets a little bit
# let's do some PCA
pca_logged <- prcomp(t(log2(data_rel+0.001)))
pca_normed <- prcomp(t(data_clr))

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

# look at clustering based on sample correlations
sample_cors_logged <- cor(log2(data_rel+0.001), use = "pairwise.complete.obs", method = "pearson")
sample_cors_normed <- cor(data_clr, use = "pairwise.complete.obs", method = "pearson")

# hierarchical clustering with complete 
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

# get the scale and colors for the heatmap
r <- c(min(sample_cors_logged),quantile(sample_cors_logged, 0.95))
breakslist <- seq(from=r[1], to=r[2], length.out=11)
cols_heat <- rev(brewer.pal(11,"RdBu"))

# plot
pheatmap(mat = sample_cors_logged, cluster_rows = clust_complete_logged, cluster_cols = clust_complete_logged, annotation_col = annotation_col, annotation_colors = annotation_colors, scale="none", cellwidth = 10, cellheight = 10, breaks=breakslist, color = cols_heat, fontsize_row = 8, fontsize_col = 8, filename = NA, main = "Complete linkage")

# CLR - normalized data
r <- c(min(sample_cors_normed),quantile(sample_cors_normed, 0.95))
breakslist <- seq(from=r[1], to=r[2], length.out=11)
pheatmap(mat = sample_cors_normed, cluster_rows = clust_complete_normed, cluster_cols = clust_complete_normed, annotation_col = annotation_col, annotation_colors = annotation_colors, scale="none", cellwidth = 10, cellheight = 10, breaks=breakslist, color = cols_heat, fontsize_row = 8, fontsize_col = 8, filename = NA, main = "Complete linkage")

# cluster using euclidean distance with complete linkage and scale the data featurewise
breakslist <- seq(from=-2, to=2, length.out=11)
pheatmap(mat = log2(data_rel+0.001), cluster_rows = T, cluster_cols = T, annotation_col = annotation_col, annotation_colors = annotation_colors, scale="row", breaks=breakslist, show_rownames = F, color = cols_heat, fontsize_row = 8, fontsize_col = 8, filename = NA, main = "Complete linkage")
pheatmap(mat = data_clr, cluster_rows = T, cluster_cols = T, annotation_col = annotation_col, annotation_colors = annotation_colors, scale="row", show_rownames = F, breaks=breakslist, color = cols_heat, fontsize_row = 8, fontsize_col = 8, filename = NA, main = "Complete linkage")

# do some differential abundance testing and calculate some
# p-values for the differential abundance of OTUs between Condition 1 and Condition 2 
# using the non-parametric Wilcoxon signed-rank test
# there are very many ways to do statistical differential abundance testing between groups / conditions
# and this is possibly one of the simplest
wilcox_logged <- apply(log2(data_rel+0.001), 1, function(z){
  w_test <- wilcox.test(x = as.numeric(z[which(sample_table$Condition=="Condition1")]), y=as.numeric(z[which(sample_table$Condition=="Condition2")]), alternative = "two.sided")$p.value
})

# gather the results and also adjust the p-values for multiple hypothesis testing and calculate a false discovery rate
wilcox_logged <- data.frame(rownames(data_rel), wilcox_logged, p.adjust(wilcox_logged, method = "fdr"))
colnames(wilcox_logged) <- c("otu", "pvalue", "fdr")

# order the results according to the signifance values
wilcox_logged <- wilcox_logged[order(as.numeric(wilcox_logged$pvalue)),]

# for CLR data
wilcox_normed <- apply(data_clr, 1, function(z){
  w_test <- wilcox.test(x = as.numeric(z[which(sample_table$Condition=="Condition1")]), y=as.numeric(z[which(sample_table$Condition=="Condition2")]), alternative = "two.sided")$p.value
})

# gather the results and also adjust the p-values for multiple hypothesis testing and calculate a false discovery rate
wilcox_normed <- data.frame(rownames(data_clr), wilcox_normed, p.adjust(wilcox_normed, method = "fdr"))
colnames(wilcox_normed) <- c("otu", "pvalue", "fdr")

# order the results according to the signifance values
wilcox_normed <- wilcox_normed[order(as.numeric(wilcox_normed$pvalue)),]

# compare the results
# get the intersect of results

# top 10
length(intersect(rownames(wilcox_logged)[1:10], rownames(wilcox_normed)[1:10])) 

# top20
length(intersect(rownames(wilcox_logged)[1:20], rownames(wilcox_normed)[1:20])) 

wilcox_logged[1:10,]
wilcox_normed[1:10,]
# are the results similar?

# let's try another test, the Reproducibility optimized test statistic (ROTS)
# similar to the wilcoxon test, ROTS is t-test alternative suitable for different types of data
# with no strict assumptions or data normality requirements as compared for example to the standard student's t-test
# furthermore, ROTS is developed specifically for omics data and has performed well in many comparisons for omics data
# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005562

# unlike the wilcoxon test, where each feature is tested separately, ROTS processes the whole data at once
# and adjusts its parameters according to the data, so it's "data driven"

# let's try it out
ROTS_logged <- ROTS(data = log2(data_rel+0.001), groups = c(rep(1,18), rep(0,18)), B = 1000, K = nrow(data_rel)/4, paired = F, progress = T, log = T, seed = 1)
ROTS_logged_df <- data.frame(d = ROTS_logged$d, logfc = ROTS_logged$logfc, pvalue = ROTS_logged$pvalue, fdr = ROTS_logged$FDR)

# order according to p-value
ROTS_logged_df <- ROTS_logged_df[order(as.numeric(ROTS_logged_df$pvalue)),]

# CLR data
ROTS_normed <- ROTS(data = data_clr, groups = c(rep(1,18), rep(0,18)), B = 1000, K = nrow(data_clr)/4, paired = F, progress = T, log = T, seed = 1)
ROTS_normed_df <- data.frame(d = ROTS_normed$d, logfc = ROTS_normed$logfc, pvalue = ROTS_normed$pvalue, fdr = ROTS_normed$FDR)

# order according to p-value
ROTS_normed_df <- ROTS_normed_df[order(as.numeric(ROTS_normed_df$pvalue)),]

# compare the results
# get the intersect of results

# top 10
length(intersect(rownames(ROTS_logged_df)[1:10], rownames(ROTS_normed_df)[1:10])) 

# top20
length(intersect(rownames(ROTS_logged_df)[1:20], rownames(ROTS_normed_df)[1:20])) 

ROTS_logged_df[1:10,]
ROTS_normed_df[1:10,]

# do you notice some differences when compared to the wilcoxon test results?
# based on these explorations and tests, which normalization would be suitable for this data in your opinion?

