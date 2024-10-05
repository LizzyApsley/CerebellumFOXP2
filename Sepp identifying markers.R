#Identifying Cell type and PC subtype markers in Sepp human snRNAseq data
#data from from Sepp et al., 2023

#libraries
library(SingleCellExperiment)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(scran)
resultsdir <- "C:/Users/Lizzy/OneDrive - King's College London/PhD/PhD_paperdrafts/rnaseqanalysis"

#load and format gene names onto data
sce <- readRDS("C:/Users/Lizzy/OneDrive - King's College London/PhD/PhD_paperdrafts/hum_sce_final.rds")
gene_info <- read.csv("C:/Users/Lizzy/OneDrive - King's College London/PhD/PhD_paperdrafts/all_genes_meta.csv")
gene_info <- filter(gene_info, species == "HUM")
colData(sce)
names(colData(sce))
rownames(sce)
rowData(sce)$gene_name <- gene_info$gene_name[match(rownames(sce), gene_info$gene_id)]
rowData(sce)

#Purkinje cell subtypes
#subset PC with subtypes only 
pc <- colData(sce) %>% 
  as.data.frame() %>% 
  rownames_to_column('cell_id') %>% 
  as_tibble() %>%
  filter(., cell_type == "Purkinje") %>%
  filter(., !(is.na(subtype)))
PC.data <- sce[,pc$cell_id]

assay(PC.data, "counts") <- assay(PC.data, "umi")
PC.data <- logNormCounts(PC.data)
PC.data@assays
markers <- findMarkers(PC.data, groups = colData(PC.data)$subtype, direction = "up",
                       pval.type = "all") 
LB <- as.data.frame(markers$Purkinje_defined_LB) %>%
  cbind(rownames(.),.)
colnames(LB)[1] <- "gene_id"
LB <- inner_join(LB, gene_info)
write.csv(LB, paste(resultsdir, "LB_markers_scran_pall.csv", sep = "/"))

EB <- as.data.frame(markers$Purkinje_defined_EB) %>%
  cbind(rownames(.),.)
colnames(EB)[1] <- "gene_id"
EB <- inner_join(EB, gene_info)

write.csv(EB, paste(resultsdir, "EB_markers_scran_pall.csv", sep = "/"))


#use scran markers to find markers for sepp cell types too
assay(sce, "counts") <- assay(sce, "umi")
sce <- logNormCounts(sce)
sce@assays
markers <- findMarkers(sce, groups = colData(sce)$cell_type, direction = "up",
                       pval.type = "all")
celltypes <- unique(colData(sce)$cell_type)
celltypes <- celltypes[-1]

for (cell in celltypes){
  m <- as.data.frame(markers[[cell]]) %>%
    cbind(rownames(.),.)
  colnames(m)[1] <- "gene_id"
  m <- inner_join(m, gene_info)
  cell <- str_remove(cell, "/")
  write.csv(m, paste(paste(resultsdir, "sepp_markers_pvaltype_all", cell, sep = "/"), "_markers_pall.csv", sep = ""))}

