# Geneset overlap - FOXP2 targets and disorders

#libraries
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(GeneOverlap)
library(RColorBrewer)

#directories
resultsdir <- "C:/Users/Lizzy/OneDrive - King's College London/PhD/PhD_paperdrafts/rnaseqanalysis"
gene_info <- readRDS(paste(resultsdir, "gene_info.rds", sep = "/"))
resmNeon <- readRDS(paste(resultsdir, "resmNeon.rds", sep = "/"))
vsd <- readRDS(paste(resultsdir, "vsd.rds", sep = "/"))

#make lists for mneon significant genes with stricter significance cut off
mneon <- list(unique(dplyr::filter(as.data.frame(resmNeon), padj <0.05 & 
                               log2FoldChange > 1.5)$symbol), 
               unique(dplyr::filter(as.data.frame(resmNeon),padj <0.05 & 
                               log2FoldChange < -1.5)$symbol))
names(mneon) <- c("pos", "neg")
length(mneon$pos) #548
length(mneon$neg) #1211
n = length(unique(as.data.frame(resmNeon)$symbol)) 
n #16110 unique genes (with symbols)

#load foxp2 targets datasets:
spiteri <- read.csv(paste(resultsdir, "Spiteri_2007.csv", sep = "/"))
head(spiteri)
#genes as symbols, expression in BG IFC and lung -> filter this to brain only
spiteri <- spiteri %>% 
  mutate(., symbol = Gene) %>%
  dplyr::filter(., BG == "X"| IFC == "X") 
spiteri_symbols <- spiteri$symbol
length(spiteri_symbols) #258 genes

vernes_2007 <- read.csv(paste(resultsdir, "Vernes_2007_ChIP.csv", sep = "/"))
head(vernes_2007) #Genes as symbols
vernes_symbols <- vernes_2007$Gene
length(vernes_symbols) #302


#all significant genes 
go.obj <- newGeneOverlap(c(mneon$pos, mneon$neg),
                         spiteri_symbols, n)
go.obj <- testGeneOverlap(go.obj)
print(go.obj) #p value 0.68
go.obj <- newGeneOverlap(c(mneon$pos, mneon$neg),
                         vernes_symbols, n)
go.obj <- testGeneOverlap(go.obj)
print(go.obj) #p value 0.79
#mneon pos genes
go.obj <- newGeneOverlap(mneon$pos,
                         spiteri_symbols, n)
go.obj <- testGeneOverlap(go.obj)
print(go.obj) #p value 0.51
go.obj <- newGeneOverlap(mneon$pos,
                         vernes_symbols, n)
go.obj <- testGeneOverlap(go.obj)
print(go.obj) #p value 0.99
#mneon neg genes
go.obj <- newGeneOverlap(mneon$neg,
                         spiteri_symbols, n)
go.obj <- testGeneOverlap(go.obj)
print(go.obj) #p value 0.74
go.obj <- newGeneOverlap(mneon$neg,
                         vernes_symbols, n)
go.obj <- testGeneOverlap(go.obj)
print(go.obj) #p value 0.34
#no significant overlaps

#find sizes of overlaps for venn diagram
#mneonpos
sum(mneon$pos %in% spiteri_symbols) #10
mneon$pos[which(mneon$pos %in% spiteri_symbols)]
#"POU4F2" "PKIA"   "CRH"    "SNCB"   "GABBR1" "MCF2"   "HOXA6"  "DPP6"  "HAS1"   "RAB27B"
sum(mneon$pos[which(mneon$pos %in% spiteri_symbols)] %in% vernes_symbols) #2
mneon$pos[which(mneon$pos[which(mneon$pos %in% spiteri_symbols)] 
                %in% vernes_symbols)] #FOXD3 and PDZRN4
sum(mneon$pos %in% vernes_symbols) #4
mneon$pos[which(mneon$pos %in% vernes_symbols)]
#"PTGER3" "SPOCK3" "GABBR1" "HAS1"  

#mneon_neg
sum(mneon$neg %in% spiteri_symbols) #19
mneon$neg[which(mneon$neg %in% spiteri_symbols)]
sum(mneon$neg[which(mneon$neg %in% spiteri_symbols)] %in% vernes_symbols) #6
mneon$neg[which(mneon$neg[which(mneon$neg %in% spiteri_symbols)] 
                %in% vernes_symbols)] 
#"PTN"    "LRP2"   "COL4A5" "HTRA1"  "LRP4"   "NOTCH2"
sum(mneon$neg %in% vernes_symbols) #25
mneon$neg[which(mneon$neg %in% vernes_symbols)]

sum(spiteri_symbols %in% vernes_symbols) #61


##################################  

#Hickey
Hickey <- read.csv(paste(resultsdir, "Hickey_2019.csv", sep = "/"))
head(Hickey)
#genes as symbols log2FC and padj
Hickey.sig <- list(unique(dplyr::filter(Hickey, padj <0.05 & 
                                       log2FoldChange > 0)$Gene), 
                unique(dplyr::filter(Hickey, padj <0.05 & 
                                       log2FoldChange < 0)$Gene))
names(Hickey.sig) <- c("H.up", "H.down")
length(Hickey.sig$H.up) #288
length(Hickey.sig$H.down) #258

gom.obj <- newGOM(mneon, Hickey.sig,n)
drawHeatmap(gom.obj)
print(gom.obj)

