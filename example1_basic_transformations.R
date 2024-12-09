# in these examples / exercises you will need the following R packages:
# tidyr
# RColorBrewer
# pheatmap
# edgeR
# umap
# phyloseq
# microbiome
# ROTS

# if you are using RStudio-Server in Puhti, all other packages should be installed except ROTS,
# which you should find for R 4.3 in the normalization workshop directory RPackages

# if you are using local RStudio in your own computer, make sure you have these packages installed
# some are from cran, some are from BioConductor
# remember the helpful setRepositories() and select 1-5 before running install.packages() for
# the above packages

# e.g. 
# setRepositories(ind = c(1,2,3,4,5))
# install.packages(c("tidyr", "RColorBrewer", "pheatmap", "edgeR", "umap", "phyloseq", "microbiome", "ROTS"))

# a very simple transformation exercise to explore how some of the common transformations
# are performed for count-based data

# this example corresponds to the examples shown in the
# video from: https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/

# let's create a simple dummy count matrix
# this data will have 4 genes and 3 samples
gene_matrix <- data.frame(matrix(nrow = 4, ncol = 3), stringsAsFactors = F)
gene_matrix[1,] <- c(10,12,30)
gene_matrix[2,] <- c(20,25,60)
gene_matrix[3,] <- c(5,8,15)
gene_matrix[4,] <- c(0,0,1)

# we will call the genes A,B,C and D.
rownames(gene_matrix) <- c("A", "B", "C", "D")

# the samples are simply called Sample1, Sample2 and Sample3
colnames(gene_matrix) <- c("Sample1", "Sample2", "Sample3")

# we also know the lenghts for these genes
gene_lengths_kb <- c(2,4,1,10) # in kilobases (kb)

# how does the data look
gene_matrix

# is there a sample having more counts / recruited reads than others?
colSums(gene_matrix)

# what about genes and gene lengths, is there a pattern?
gene_lengths_kb
gene_matrix
rowSums(gene_matrix)

### let's do some transformations
### first, let's get RPM (reads per million) / CPM (counts per million) 

# we need all the read counts per sample - let's assume that in this dummy matrix all the genes 
# in the sample are accounted for
total_reads <- colSums(gene_matrix)

# normally, we would use a scaling factor of 10^6 to get the reads / counts "per million"
# however, for the sake of nice numbers in this practice, let's just use a scaling factor of 10
# instead of 10^6. No difference in the calculations, just the scale of the numbers will be different

# so our "per million" scaling factor for each sample will be:
sample_scaling_factors <- total_reads / 10

# to get the RPM / CPM data, we just divide the gene counts in each sample with 
# the corresponding "per million" scaling factor for that sample
rpm_data <- t(t(gene_matrix) / sample_scaling_factors) # transpose to divide the samples and then transpose back

# this transformation / "normalization" accounts for the differences in library sizes or sequencing
# depth between the samples.

# explore this data
rpm_data

# look at the gene abundances between samples (rows) now compared to the raw counts, did something happen?
gene_matrix

### Next let's get RPK (reads per kilobase) data

# to get RPK data, we just divide each gene by it's length in kilobases.
rpk_data <- gene_matrix / gene_lengths_kb

# explore this data
rpk_data

# look at the gene abundances within samples (columns) now compared to the raw counts, did something happen?
gene_matrix

### RPKM (reads per kilobase million)
# typically we want to adjust both for variations in gene length and library size
# RPKM was developed for this purpose for RNA sequencing data.

# to get RPKM, we use the RPM data and divide by gene lengths
rpkm_data <- rpm_data / gene_lengths_kb

# explore this data
rpkm_data

# look at the gene abundances both between and within samples (columns) now compared to the raw counts,
# did something happen?
gene_matrix

### TPM (Transcripts per million)
# while RPKM was developed to normalize for both gene length and differences in sequencing depth
# there are some issues due to which RPKM shouldn't perhaps used for between sample comparison 
# of gene abundances (we'll soon look into this a bit more)

# another very similar transformation to RPKM, is TPM, where the order of operations is slightly
# different resulting in data which should enable between sample comparison of gene abundances better 
# than RPKM

# to get TPM, we use the gene length normalized RPK data
# next, we want to normalize for differences in library sizes or sequencing depth

# we need the read counts for all genes per sample normalized to their respective gene lengths (in kilobases)
# we can get this information from the RPK data:
total_rpks <- colSums(rpk_data)

# we then again get the "per million" scaling factors for each sample 
# (except per million (10^6), we use "per ten"in this exercise to get nicer numbers, but
# nothing really changes except the scale of the resulting numbers)
rpk_sample_scaling_factors <- total_rpks / 10

# then, TPM data is simply:
tpm_data <- t(t(rpk_data) / rpk_sample_scaling_factors) # transpose to divide the samples and then transpose back

# explore this data
tpm_data

# again, look at the gene abundances both between and within samples (columns) now compared to the raw counts,
# did something happen?
gene_matrix

# what's the difference between RPKM and TPM then?
rpkm_data
tpm_data

# both RPKM and TPM adjust for biases caused by differences in 
# gene length and library sizes / sequencing depth

# but the sums of total normalized reads for each sample are very different:
colSums(rpkm_data)
colSums(tpm_data)
# for RPKM, we have a different value for each sample. For TPM we have the same value for each sample

# now consider gene A. In TPM data:
tpm_data[1,]
# gene A has a value of approximately 3.333 in sample 1 and 3.326 sample 3
# This means that the relative expression of gene A in sample 1, 3.333 / per million 
# (normally, now 10, because of the different scaling factor we used),
# is slightly larger than in sample 3, 3.326 / per million.
# Of all the reads that mapped to sample 1, a slightly larger proportion of them
# mapped to gene A as compared to the proportion of reads mapping to gene A of all the
# reads that mapped to sample 3

# if you look more closely, you see some reads mapping to gene D in sample 3 but none sample 1
tpm_data

# so as the samples are normalized to have the same overall total, different proportions mean 
# different relative abundance for a gene. 

# however, with RPKM, it is harder to compare the genes between samples as each sample has a
# different total 
colSums(rpkm_data)

# while in the RPKM data, gene A in sample 1 has a value of approximately 1.429 and a value
# of 1.415 in sample 3, sample 1 has a total of 4.285714 and sample 3 has a total of 4.254717. 
# Thus the proportions are not directly comparable as they are proportions of a different whole

# however... if one has RPKM data for all the genes in the data, it can be easily converted 
# to tpm data, similar to forming tpm data from RPK data:
total_rpkms <- colSums(rpkm_data)
rpkm_sample_scaling_factors <- total_rpkms / 10
tpm_data_rpkm_converted <- t(t(rpkm_data) / rpkm_sample_scaling_factors) # transpose to divide the samples and then transpose back

# compare 
tpm_data
tpm_data_rpkm_converted
