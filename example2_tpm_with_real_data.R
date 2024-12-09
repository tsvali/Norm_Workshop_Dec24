# a simple example to calculate TPM values for real data.

# load the needed data
load("example2_data.RData")

# In the "all_sample_data table", we have counts for the 51 metabolic marker genes 
# The 51 metabolic marker gene database has collected / compiled protein sequences for
# marker genes for different metabolic functions (e.g. methanogenesis, denitricication)

# for this project, the eukaryotic genes have first been removed from the 
# 51 metabolic marker gene database. Following the metagenomic read from all the samples
# have been aligned to the database with specific criteria

# Following, the results have been filtered with specific thresholds for increased reliability
# and the results for both read pairs for the same sample have been combined (in such a way, that
# the same gene is only counted once, if both reads align to the same gene. However if the paired
# reads align to entirely different functional genes, both hits are counted)

# reference for the database
# Leung, Pok Man; Greening, Chris (2021). Compiled Greening lab metabolic marker gene databases. Monash University. Online resource. https://doi.org/10.26180/14431208.v1
# https://doi.org/10.26180/14431208.v1

# load libraries
library(tidyr)

# explore the data
str(all_sample_data)

# out of this, we would like to have TPM - transformed marker gene data, aggregated to marker gene 
# level, with genes as rows and samples as columns

# first get RPK, read per kilobase, values
# simply divide counts with the gene lenghts in kilobases
all_sample_data$RPK <- all_sample_data$Count / (all_sample_data$Gene_length/1000) # divide length in bases with 1000 to get kilobases

# now we could aggregate the RPK data at metabolic marker gene and sample level
rpk_data <- aggregate(RPK ~ Sample+Gene, all_sample_data, sum)

# spread this data into gene X sample matrix
rpk_data <- spread(rpk_data, Sample, RPK)

# remove the gene column and use it as rownames
rownames(rpk_data) <- rpk_data$Gene
rpk_data <- rpk_data[,-which(colnames(rpk_data)%in%"Gene")]

# explore
str(rpk_data)
# looks good

# next, we want to adjust the RPK values in each sample with the total sum of RPK values
# of all genes in the sample 

# for this data, all the prokaryotic / bacterial genes in each sample have been estimated
# by aligning the metagenomic reads for samples to the KEGG prokaryotic database

# this gives us an estimate for the counts of all the bacterial genes in the samples
# following, the gene counts have been divided with the respective gene lengths to get the RPK for all the genes
# then for each sample, the RPK values for all genes have been summed together to create sample total RPKs

# these values are stored in the variable "kegg_prok_genes_sample_sum_rpks"
kegg_prok_genes_sample_sum_rpks

# we can use these values to normalize for differences in library sizes in 
# the gene length normalized RPK data.

# make sure these values are in same order as samples in the RPK data
kegg_prok_genes_sample_sum_rpks <- kegg_prok_genes_sample_sum_rpks[colnames(rpk_data)]

# formulate the per million scaling factors from these
rpk_scaling_factors <- kegg_prok_genes_sample_sum_rpks / 10^6

# get TPM data
tpm_data <- t(t(rpk_data) / rpk_scaling_factors)

# change NA:s to 0:s, could also be kept as NA:s, depends on how the data will be used / further processed
tpm_data[is.na(tpm_data)] <- 0

# how does this data look
colSums(tpm_data)
boxplot(log2(tpm_data+1)) # log2 transform for a nicer scale for the boxplot and add offset / pseudocount to deal with 0. Log2(0)=not defined, log2(1)=0
# looks pretty decent

# this looks now somewhat different from the tpm data in exercise 1. Not all the sample
# totals are equal. This is because in this data table, we only have a subset of all the genes
# (or estimated all genes), to which the expression of all the genes are normalized to

# if we would have all the genes from the KEGG database that were used in estimating the total RPKs for 
# each sample, all the sample sums would be 10^6.

# in the table "sample_table" we have some metadata for the study
# use it to do some simple exploration of the TPM data

# first parse sample names to match metadata
colnames(tpm_data) <- paste("P",unlist(lapply(colnames(tpm_data), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")

# organize similarly
tpm_data <- tpm_data[,rownames(sample_table)]

# do simple PCA
pca_logged <- prcomp(t(log2(tpm_data+1)))
# logged data more suitable for pca, more normally distributed, more equal variances for genes.

# define some well differentiating colors
cols_diff <- c("purple","cornflowerblue","mistyrose2","darkgrey", "deeppink","brown","violet","firebrick1","chocolate3","forestgreen","darkslateblue","darkolivegreen2", "yellow")

# color according to experimental block
cols_pca <- cols_diff[as.numeric(sample_table$Block)]

plot(pca_logged$x[,1:2], col = cols_pca, pch=19, cex=2, xlab = paste("PC1 ", "(",(100*round((pca_logged$sdev[1])^2/(sum(pca_logged$sdev^2)),3)),"%)" , sep = ""),
     ylab = paste("PC2 ", "(",(100*round((pca_logged$sdev[2])^2/(sum(pca_logged$sdev^2)),3)),"%)" , sep = ""), main = "Block")
# block seems to have some slight effect - which is natural as the samples from the same block are adjacent - but doesn't look too controlling


# color according to the condition of interest
cols_pca <- rep("blue", nrow(sample_table))
cols_pca[which(sample_table$Condition=="Condition2")] <- "red"

plot(pca_logged$x[,1:2], col = cols_pca, pch=19, cex=2, xlab = paste("PC1 ", "(",(100*round((pca_logged$sdev[1])^2/(sum(pca_logged$sdev^2)),3)),"%)" , sep = ""),
     ylab = paste("PC2 ", "(",(100*round((pca_logged$sdev[2])^2/(sum(pca_logged$sdev^2)),3)),"%)" , sep = ""), main = "Condition")
# condition seems to have a rather clear effect
# these were just examples of some preliminary explorations. In practice I would do a lot more data exploration

###
# what if we don't have the sample specific RPK totals (the length normalized counts for all the genes in samples)?
# or if we think that the estimate used here (or a similar estimate) is not appropriate?

# could we get TPM data just by using the sample totals in the used RPK data (for marker genes in this case)?
total_rpks <- colSums(rpk_data, na.rm = T)
rpk_scaling_factors <- total_rpks / 10^6

# get TPM data
tpm_data_rpk_sums <- t(t(rpk_data) / rpk_scaling_factors)

# change NA:s to 0:s, could also be kept as NA:s, depends on how the data will be used.
tpm_data_rpk_sums[is.na(tpm_data_rpk_sums)] <- 0

# how does this data look
colSums(tpm_data_rpk_sums)
# now all the sample totals are equal and correspond to one million.

# however, the interpration for this data is very different than for the first tpm data using sample total RPKs for scaling
# the first data was adjusted to an estimate for all the genes in the sample
# the latter data is adjusted to all the marker genes in the sample

# thus the abundance of a marker gene (e.g. AtpA) in a sample is proportional to an estimate of total (bacterial) gene
# expression in the sample in the first TPM data, while in the latter TPM data, the abundance of a marker gene in a sample
# is proportional to the expression of all the marker genes (instead of all genes) in the sample

# yet another way could be by using the total sequences in the samples
# one could think of using the total sequences in a sample as a surrogate for total gene counts in 
# a sample. Even though the number of sequences will typically be larger than the number of gene counts 
# (as each sequence / read in a sample does not necessarily map to a gene),
# there is typically a very good correlation between the two. For example, 0.99 in this data

# let's compose a tpm data using total sequences for each sample
total_sequences <- total_sequences[colnames(rpk_data)]

# to get total rpks for samples we could use median gene length as an estimate of gene length
# then total sequences / median gene length would be an estimate of total rpks in the samples
# or perhaps one could also use sample specific median lengths instead of a global common median
# for here, let's use the global common median for simplicity

# get median gene length
med_gene_length_kb <- median(all_sample_data$Gene_length/1000)

# get total rpks and rpk scaling factors
total_rpks <- total_sequences / med_gene_length_kb
rpk_scaling_factors <- total_rpks / 10^6

# get TPM data based on this
tpm_data_tot_seqs <- t(t(rpk_data) / rpk_scaling_factors)

# change NA:s to 0:s, could also be kept as NA:s, depends on how the data will be used.
tpm_data_tot_seqs[is.na(tpm_data_tot_seqs)] <- 0

# how does this data look
colSums(tpm_data_tot_seqs)
# this looks now more similar to the first TPM data as we are scaling relative to the total
# number of sequences in the sample

# compare the TPM datas a little bit

# parse sample names and organize all similarly
colnames(tpm_data_rpk_sums) <- paste("P",unlist(lapply(colnames(tpm_data_rpk_sums), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")
colnames(tpm_data_tot_seqs) <- paste("P",unlist(lapply(colnames(tpm_data_tot_seqs), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")

tpm_data_rpk_sums <- tpm_data_rpk_sums[,rownames(sample_table)]
tpm_data_tot_seqs <- tpm_data_tot_seqs[, rownames(sample_table)]

# correlate genes over samples in the different datasets
all_tpm_datas <- list(tpm_data, tpm_data_rpk_sums, tpm_data_tot_seqs)
names(all_tpm_datas) <- c("est_total_rpk_kegg", "marker_total_rpk", "total_seqs")

# the combination of datas we want to compare
all_combs <- combn(x = seq(1:3), m = 2)
gene_cors <- data.frame(matrix(nrow = 3, ncol = nrow(tpm_data)), stringsAsFactors = F)

namn <- character(0)
for(i in 1:ncol(all_combs)){
  for(j in 1:nrow(tpm_data)){
    gene_cors[i,j] <- cor(as.numeric(all_tpm_datas[[all_combs[1,i]]][j,]),as.numeric(all_tpm_datas[[all_combs[2,i]]][j,]), method = "pearson", use = "pairwise.complete.obs") 
  }
  namn <- c(namn, paste(names(all_tpm_datas)[all_combs[1,i]], " vs. ",names(all_tpm_datas)[all_combs[2,i]], sep=""))
}

# name
colnames(gene_cors) <- rownames(tpm_data)
rownames(gene_cors) <- namn

# what does this look like?
gene_cors
summary(as.numeric(gene_cors[1,]))
summary(as.numeric(gene_cors[2,]))
summary(as.numeric(gene_cors[3,]))

# for most genes the data looks really similar, but there are a few exceptions.

### yet another option could be to use the sample totals for 16S rRNA genes
# e.g. total counts for the 16S rRNA genes for the samples divided by the length of the 16S rRNA gene(s).
# however, as already discussed, the interpration might be somewhat different 