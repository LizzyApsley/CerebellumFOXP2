#mneon rnaseq plotting

#libraries
library("DESeq2")
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)

#directories:
resultsdir <- "C:/Users/Lizzy/OneDrive - King's College London/PhD/PhD_paperdrafts/rnaseqanalysis"
vsd <- readRDS(paste(resultsdir, "vsd.rds", sep = "/"))
resmNeon <- readRDS(paste(resultsdir, "resmNeon.rds", sep = "/"))
gene_info <- readRDS(paste(resultsdir, "gene_info.rds", sep = "/"))


#PCA plot - ggplot plotting:
pcaData <- plotPCA(vsd, intgroup=c("mNEON", "Diff"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=mNEON, shape=Diff)) +
  geom_point(size=6) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_classic() +
  theme(axis.text = element_text(size = 16, colour = "black"), 
        axis.title = element_text(size = 18)) +
  scale_colour_manual(values = c( "#666666","#66A61E"), 
                      breaks = c("neg", "pos"))


#volcano plot showing most changed genes
p <- 0.05  #threshold for colour min adjusted p value
fc <- 1.5 #threshold for colour - min absolute log2fc

plot.results <- resmNeon %>%
  dplyr::filter(., baseMean != 0) %>%
  mutate(., colour = ifelse(symbol %in% c("FOXP1" ,"FOXP2", "FOXP4"), "FOXP",
                            ifelse(padj >p, "nonsig",
                                   ifelse(log2FoldChange >fc, "sig.up", 
                                          ifelse(log2FoldChange < -fc,"sig.down", "sig.smallchange")))))

label_genes <- rbind(slice_max(dplyr::filter(plot.results, log2FoldChange >0 & padj <0.05), log2FoldChange, n = 12),
                     slice_min(dplyr::filter(plot.results, log2FoldChange <0 & padj <0.05), log2FoldChange, n = 12),
                     slice_min(dplyr::filter(plot.results, log2FoldChange >0), padj, n = 12),
                     slice_min(dplyr::filter(plot.results, log2FoldChange <0), padj, n = 12),
                     dplyr::filter(plot.results, symbol %in% c("FOXP1" ,"FOXP2", "FOXP4"))) %>%
  unique(.)

ggplot(plot.results, aes(x = log2FoldChange, y = -log10(padj)))+
  geom_point(aes(colour = colour), size = 2, alpha = 0.3) +theme_classic() + 
  geom_hline(yintercept = -log10(p), linetype = "dashed")+
  geom_vline(xintercept = fc, linetype = "dashed")+ 
  geom_vline(xintercept = -fc, linetype = "dashed")+
  geom_text_repel(data = label_genes, aes(label = symbol, colour = colour), nudge_y = 1, 
                  size = 3, max.overlaps = 17)+
  scale_color_manual(breaks = c("sig.smallchange", "sig.up", "sig.down", "nonsig", "FOXP"),
                     values=c("#666666","#66A61E","#1F78B4","#666666", "purple")) +
  theme(axis.text = element_text(size = 12, colour = "black"), 
        axis.title = element_text(size = 12),
        legend.position = "none")

#heatmap showing marker for major cerebellar cell types
GC <- c("PAX6", "NEUROD1", "BARHL1") 
VZ_neuroblast <- c("PTF1A", "KIRREL2")
Purkinje <- c("SKOR2", "ESRRB", "PTGER3") 
progenitor <- c("PAX3", "NOTCH1", "TOP2A")  
oligo <- c("SOX10", "MAG")  
interneuron <- c("KIT", "SLC6A5", "PAX2") 
UBC <- c("EOMES")      
GABA_DN <- c("ZFHX3", "TOX", "SLC24A4")
glut_DN <- c("MEIS2", "NTNG1", "ZNF804A") 
astrocyte <- c("AQP4", "TNC", "EGFR") 

genes <- c(progenitor, VZ_neuroblast, Purkinje, interneuron, GABA_DN, GC, UBC, glut_DN,
           oligo, astrocyte)
my.genes <- dplyr::filter(gene_info, symbol %in% genes)%>%
  dplyr::filter(., gene_id %in% rownames(vsd)) %>%
  mutate(., names = paste(symbol, gene_id, sep = "_")) %>%
  arrange(., match(symbol, genes))

mat <- assay(vsd)[my.genes$gene_id,]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("Diff","mNEON")])
temp <- as.data.frame(mat)
temp$gene_id <- rownames(temp)
temp <- inner_join(temp, gene_info) %>%
  mutate(., names = paste(symbol, gene_id, sep = "_"))
rownames(temp) <- temp$names
toplot <- temp[my.genes$names,1:8]

pheatmap(toplot, annotation_col = dplyr::select(anno, mNEON), 
         annotation_colors = list(mNEON = c(neg = "#666666", 
                                            pos = "#66A61E")),
         show_colnames = F, cluster_rows = F, 
         labels_row = as.expression(lapply(temp$symbol, function(a) bquote(italic(.(a)))))) 

#results for these genes
resultsforplot <- dplyr::filter(resmNeon, gene_id %in% my.genes$gene_id) %>%
  arrange(., match(gene_id, my.genes$gene_id)) %>%
  mutate(., sig = ifelse(log2FoldChange <0 & padj <0.05, "sig.down",
                         ifelse(log2FoldChange >0 & padj <0.05, "sig.up", "non.sig")))

#CAS genes
CAS <- read.csv(paste(resultsdir, "CAS_geneset.csv", sep = "/"))
head(CAS)
my.genes <- dplyr::filter(gene_info, symbol %in% CAS$Gene)%>%
  dplyr::filter(., gene_id %in% rownames(vsd)) %>%
  mutate(., names = paste(symbol, gene_id, sep = "_")) %>%
  arrange(., match(symbol, genes))
mat <- assay(vsd)[my.genes$gene_id,]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("Diff","mNEON")])
temp <- as.data.frame(mat)
temp$gene_id <- rownames(temp)
temp <- inner_join(temp, gene_info) %>%
  mutate(., names = paste(symbol, gene_id, sep = "_"))
rownames(temp) <- temp$names
toplot <- temp[my.genes$names,1:8]

pheatmap(toplot, annotation_col = dplyr::select(anno, mNEON), 
         annotation_colors = list(mNEON = c(neg = "#666666", 
                                            pos = "#66A61E")),
         show_colnames = F, 
         labels_row = as.expression(lapply(temp$symbol, function(a) bquote(italic(.(a)))))) 

resultsforplot <- dplyr::filter(resmNeon, gene_id %in% my.genes$gene_id) %>%
  arrange(., match(gene_id, my.genes$gene_id)) %>%
  mutate(., sig = ifelse(log2FoldChange <0 & padj <0.05, "sig.down",
                         ifelse(log2FoldChange >0 & padj <0.05, "sig.up", "non.sig")))


#ASD genes
sfari <- read.csv(paste(resultsdir, "SFARI-Gene_human-gene-scores_01-23-2023release_03-08-2023export.csv", sep = "/"))
head(sfari)
sfari <- dplyr::select(sfari, gene.symbol, gene.score) %>%
  dplyr::filter(., gene.score == 1)

#plot ASD risk genes differentially expressed in mNeon sorted RNAseq
mneon <- list(unique(dplyr::filter(as.data.frame(resmNeon), padj <0.05 & 
                                     log2FoldChange > 1.5)$symbol), 
              unique(dplyr::filter(as.data.frame(resmNeon),padj <0.05 & 
                                     log2FoldChange < -1.5)$symbol))
names(mneon) <- c("pos", "neg")

genes <- c(sfari$gene.symbol[sfari$gene.symbol %in% mneon$pos], sfari$gene.symbol[sfari$gene.symbol %in% mneon$neg])

my.genes <- dplyr::filter(gene_info, symbol %in% genes)%>%
  dplyr::filter(., gene_id %in% rownames(vsd)) %>%
  mutate(., names = paste(symbol, gene_id, sep = "_")) %>%
  arrange(., match(symbol, genes))

mat <- assay(vsd)[my.genes$gene_id,]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("Diff","mNEON")])
temp <- as.data.frame(mat)
temp$gene_id <- rownames(temp)
temp <- inner_join(temp, gene_info) %>%
  mutate(., names = paste(symbol, gene_id, sep = "_"))
rownames(temp) <- temp$names
toplot <- temp[my.genes$names,1:8]

pheatmap(toplot, annotation_col = dplyr::select(anno, mNEON), 
         annotation_colors = list(mNEON = c(neg = "#666666", 
                                            pos = "#66A61E")),
         show_colnames = F, cluster_rows = F, 
         labels_row = as.expression(lapply(temp$symbol, function(a) bquote(italic(.(a)))))) 
