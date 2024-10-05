#mneon rnaseq analysis - processing and differential gene expression

#libraries
library(AnnotationHub)
library("DESeq2")
library(tidyverse)
library(GenomicFeatures)
library(tximport)

#directories
datadir <-  "C:/Users/Lizzy/OneDrive - King's College London/PhD/PhD_mNeonsortedRNAseq" 
resultsdir <- "C:/Users/Lizzy/OneDrive - King's College London/PhD/PhD_paperdrafts/rnaseqanalysis"

#####Import data using tximeta####

#sample info
coldata <- read.csv(paste(datadir, "sample_info.csv", sep = "/"), header = T)
colnames(coldata)[1] <- "names"
head(coldata)

#check all files exist
files <- file.path(datadir,coldata$names, "abundance.h5")
names(files) <- coldata$names
file.exists(files)

#setting up transcripts to genes tx2gene 
ah <- AnnotationHub()
query(ah, c("EnsDb", "v108", "Homo Sapiens"))  #latest version 108 used for alignment
txdb <- ah[["AH109336"]]

t <- transcripts(txdb, return.type = "DataFrame") #get transcript info
colnames(t) #preview what info is here
t <- as.data.frame(t)
# select format required for tx2gene (transcript ID then gene ID)
tx2gene <- dplyr::select(t, tx_id_version, gene_id)
saveRDS(tx2gene, paste(resultsdir, "tx2gene.rds", sep = "/")) #save for reference

#make gene key - to add gene symbol to results later 
#using gene info from earlier txdb select ens_id and symbol
g <- genes(txdb, return.type = "DataFrame")
colnames(g)
g <- as.data.frame(g)
gene_info <- dplyr::select(g, gene_id, symbol)
#save gene key for later
saveRDS(gene_info, paste(resultsdir, "gene_info.rds", sep = "/"))

txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = F) 
names(txi.kallisto)

ddsTxi <- DESeqDataSetFromTximport(txi.kallisto,
                                   colData = coldata,
                                   design = ~ Diff + mNEON)

#pre-filtering - ensure at least 10 reads in at least 4 samples
dim(ddsTxi)#original ddsTxi has 40526 features x 8 samples
keep <- rowSums(counts(ddsTxi) >= 10) >= 4
dds <- ddsTxi[keep,]
ddsTxi <- dds
dim(ddsTxi)
#after filtering have 17884 features x 8 samples

#order factors so mNEONneg is used for reference
ddsTxi$mNEON <- factor(ddsTxi$mNEON, levels = c("neg","pos"))

###########################################
#perform DESEQ analysis and save results

#include Diff for batch correction and mNEON as variable of interest
design(ddsTxi) <- formula(~ Diff + mNEON)
ddsTxi <- DESeq(ddsTxi)

saveRDS(ddsTxi, paste(resultsdir, "ddsTxi.rds", sep = "/")) #save for later

res <- results(ddsTxi, alpha = 0.05) #results table with adj. p-value threshold 0.05
res
summary(res) #summary of comparison, tests, number of sig genes etc

#perform shrinkage to remove effect of very low expressed genes
resultsNames(ddsTxi)
resLFC <- lfcShrink(ddsTxi, coef="mNEON_pos_vs_neg", type="apeglm")
resLFC
summary(resLFC) #summary of comparison, tests, number of sig genes etc
saveRDS(resLFC, paste(resultsdir, "resLFC.rds", sep = "/"))

resmNeon <- as.data.frame(resLFC[order(resLFC$padj),])
resmNeon <- cbind(rownames(resmNeon), resmNeon)
colnames(resmNeon)[1] <- "gene_id"
resmNeon <- inner_join(resmNeon, gene_info)

#how many sig
length(unique(dplyr::filter(resmNeon, padj <0.05)$gene_id)) #7644
length(unique(dplyr::filter(resmNeon, padj <0.05 & log2FoldChange > 0)$gene_id))
#3655 up
length(unique(dplyr::filter(resmNeon, padj <0.05 & log2FoldChange < 0)$gene_id))
#3989

write.csv(resmNeon, paste(resultsdir, "mNEONDE_logFCshrink_results.csv", sep = "/")) #save results file
saveRDS(resmNeon, paste(resultsdir, "resmNeon.rds", sep = "/")) #save for later

#Transform data for clustering and plotting
vsd <- vst(ddsTxi, blind=FALSE) #blind false so takes into account grouping of sample 
saveRDS(vsd, paste(resultsdir, "vsd.rds", sep = "/")) #save for later



