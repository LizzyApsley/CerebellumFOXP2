#create files for GSEA geneset analysis and plot output results

#libraries
library(reshape2)
library(tidyverse)
library(SingleCellExperiment)
library(scran)

#directories
resultsdir <- "C:/Users/Lizzy/OneDrive - King's College London/PhD/PhD_paperdrafts/rnaseqanalysis"

##########################################

#preparing geneset files for GSEA

#Cerebellar disorder genesets 
#From Aldinger et al. 2021
disorders <- read.csv(paste(resultsdir, "Aldinger_disorder_genesets.csv", sep = "/"))
head(disorders) #gene symbols
length(unique(disorders$Gene)) #510 genes
table(disorders$Condition) #ALZ, ASD, CBLM, ID, JS, SCA
#CAS genes 
#combined from Eising et al., 2019; Hildebrand et al., 2020; Kaspi et al., 2022
CAS <- read.csv(paste(resultsdir, "CAS_geneset.csv", sep = "/"))
head(CAS)
dim(CAS)
#SFARI ASD risk genes
sfari <- read.csv(paste(resultsdir, "SFARI-Gene_human-gene-scores_01-23-2023release_03-08-2023export.csv", sep = "/"))
head(sfari)
dim(sfari)
#combine
temp <- cbind(sfari$gene.symbol, filter(disorders, Condition == "ID")$Gene) %>%
  cbind(., filter(disorders, Condition == "ALZ")$Gene) %>%
  cbind(., filter(disorders, Condition == "CBLM")$Gene )%>%
  cbind(., filter(disorders, Condition == "JOUBERT")$Gene )%>%
  cbind(., filter(disorders, Condition == "SCA")$Gene ) %>%
  cbind(., CAS$Gene) 
colnames(temp) <- c("ASD", "ID", "ALZ", "CBLM", "JOUBERT", "SCA", "CAS")
#save
write.csv(temp, paste(resultsdir, "genesets.csv", sep = "/"))
#duplicates removed and reformatted to .gmx file 

#Cerebellar region marker genes
#data from Aldinger et al., 2021
aldLCM <- read.csv(paste(resultsdir, "aldingerLCMresults.csv", sep= "/"))
head(aldLCM)
RL <- dplyr::filter(aldLCM, log2.FC...RL.vs.Bulk >1.5 & FDR.adj.P.value..RL.vs.Bulk <0.05) %>%
  arrange(., desc(log2.FC...RL.vs.Bulk))
EGL <- dplyr::filter(aldLCM, log2.FC...EGL.vs.Bulk >1.5 & FDR.adj.P.value..EGL.vs.Bulk <0.05) %>%
  arrange(., desc(log2.FC...EGL.vs.Bulk))
PC <- dplyr::filter(aldLCM, log2.FC...PCL.vs.Bulk >1.5 & FDR.adj.P.value..PCL.vs.Bulk <0.05) %>%
  arrange(., desc(log2.FC...PCL.vs.Bulk))
#combine longest first
temp <- cbind(RL$Gene, EGL$Gene) %>%
  cbind(., PC$Gene)
#save
write.csv(temp, paste(resultsdir, "ALD-LCM-GENESETS.CSV", sep = "/"))
#duplicates removed and reformatted to .gmx file 

#NDD disorder DEG from Gandal et al., 2018
#downloaded from http://resource.psychencode.org/
disorders_DEG <- read.csv(paste(resultsdir, "DER-13_Disorder_DEX_Genes.csv", 
                                sep = "/"))
head(disorders_DEG)
table(disorders_DEG$Disorder.DGE_RegulationDirection)
#combine longest first
temp <- cbind(filter(disorders_DEG, Disorder.DGE_RegulationDirection == "SCZ.DGE_up")$Gene_Name,
              filter(disorders_DEG, Disorder.DGE_RegulationDirection == "SCZ.DGE_down")$Gene_Name) %>%
  cbind(., filter(disorders_DEG, Disorder.DGE_RegulationDirection == "ASD.DGE_up")$Gene_Name) %>%
  cbind(., filter(disorders_DEG, Disorder.DGE_RegulationDirection == "ASD.DGE_down")$Gene_Name) %>%
  cbind(., filter(disorders_DEG, Disorder.DGE_RegulationDirection == "BD.DGE_up")$Gene_Name) %>%
  cbind(., filter(disorders_DEG, Disorder.DGE_RegulationDirection == "BD.DGE_down")$Gene_Name)
colnames(temp) <- c("SCZ.up", "SCZ.down", "ASD.up", "ASD.down", "BD.up", "BD.down")  
write.csv(temp, paste(resultsdir, "disorderDEG_genesets.csv", sep = "/"))
#duplicates removed and reformatted to .gmx file 

#Cerebellar cell type markers from Sepp et al., 2023

sce <- readRDS("C:/Users/Lizzy/OneDrive - King's College London/PhD/PhD_paperdrafts/hum_sce_final.rds")
#use scran markers to find markers cell types too
assay(sce, "counts") <- assay(sce, "umi")
sce <- logNormCounts(sce)
sce@assays
markers <- findMarkers(sce, groups = colData(sce)$cell_type, direction = "up",
                       pval.type = "all")
celltypes <- unique(colData(sce)$cell_type)[-1]
gene_info <- read.csv("C:/Users/Lizzy/OneDrive - King's College London/PhD/PhD_paperdrafts/all_genes_meta.csv")
gene_info <- filter(gene_info, species == "HUM")
sepplist <- list()
for(c in celltypes){
  t <-  dplyr::filter(gene_info, gene_id %in% rownames(dplyr::filter(as.data.frame(markers[[c]]),
                                   FDR <0.05 & summary.logFC >0)))$gene_name
  temp <- c(sepplist, list(t))
  sepplist <- temp}
names(sepplist) <- celltypes

lengths <- sapply(sepplist, length)
names(which(lengths == max(lengths)))
maxlength = max(lengths)
#oligo longest (1181)

sepplist.df <- data.frame()
for (c in celltypes){
  temp <- c(sepplist[[c]], rep(NA, maxlength - length(sepplist[[c]])))
  sepplist.df <- rbind(sepplist.df, temp)}
rownames(sepplist.df) <- celltypes
t.sepplist.df <- t(sepplist.df)

write.csv(t.sepplist.df, paste(resultsdir, "sepp-markers-symbols.csv", sep = "/"))
#duplicates removed and reformatted to gmx file outside of R.

