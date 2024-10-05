#Identifying Cell type and PC subtype markers in Sepp mouse snRNAseq data
#data from from Sepp et al., 2023

#libraries
library(SingleCellExperiment)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(scran)
resultsdir <- "C:/Users/Lizzy/OneDrive - King's College London/PhD/PhD_paperdrafts/rnaseqanalysis"

#load sepp mouse data
ms.sce <- readRDS("C:/Users/Lizzy/OneDrive - King's College London/PhD/PhD_paperdrafts/mou_sce_final.rds")
gene_info <- read.csv("C:/Users/Lizzy/OneDrive - King's College London/PhD/PhD_paperdrafts/all_genes_meta.csv")
gene_info <- filter(gene_info, species == "MOU")
names(colData(ms.sce))
rownames(ms.sce)
rowData(ms.sce)$gene_name <- gene_info$gene_name[match(rownames(ms.sce), gene_info$gene_id)]
rowData(ms.sce)

#Purkinje cell subtypes
#subset PC with subtypes only 
pc <- colData(ms.sce) %>% 
  as.data.frame() %>% 
  rownames_to_column('cell_id') %>% 
  as_tibble() %>%
  filter(., cell_type == "Purkinje") %>%
  filter(., !(is.na(subtype)))
PC.data <- ms.sce[,pc$cell_id]

assay(PC.data, "counts") <- assay(PC.data, "umi")
PC.data <- logNormCounts(PC.data)
PC.data@assays
markers <- findMarkers(PC.data, groups = colData(PC.data)$subtype, direction = "up",
                       pval.type = "all") 

Cdh9 <- as.data.frame(markers$Purkinje_defined_Cdh9) %>%
  cbind(rownames(.),.)
colnames(Cdh9)[1] <- "gene_id"
Cdh9 <- inner_join(Cdh9, gene_info)
write.csv(Cdh9, paste(resultsdir, "Cdh9_ms_markers_scran_pall.csv", sep = "/"))

Etv1 <- as.data.frame(markers$Purkinje_defined_Etv1) %>%
  cbind(rownames(.),.)
colnames(Etv1)[1] <- "gene_id"
Etv1 <- inner_join(Etv1, gene_info)
write.csv(Etv1, paste(resultsdir, "Etv1_ms_markers_scran_pall.csv", sep = "/"))

Foxp1 <- as.data.frame(markers$Purkinje_defined_Foxp1) %>%
  cbind(rownames(.),.)
colnames(Foxp1)[1] <- "gene_id"
Foxp1 <- inner_join(Foxp1, gene_info)
write.csv(Foxp1, paste(resultsdir, "Foxp1_ms_markers_scran_pall.csv", sep = "/"))

Rorb <- as.data.frame(markers$Purkinje_defined_Rorb) %>%
  cbind(rownames(.),.)
colnames(Rorb)[1] <- "gene_id"
Rorb <- inner_join(Rorb, gene_info)
write.csv(Rorb, paste(resultsdir, "Rorb_ms_markers_scran_pall.csv", sep = "/"))



