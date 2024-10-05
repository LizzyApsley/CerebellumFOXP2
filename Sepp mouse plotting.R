#plotting mouse developing cerebellar snRNAseq data from Sepp et al., 2023

#libraries
library(SingleCellExperiment)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

#plotting sepp mouse data
ms.sce <- readRDS("C:/Users/Lizzy/OneDrive - King's College London/PhD/PhD_paperdrafts/mou_sce_final.rds")
gene_info <- read.csv("C:/Users/Lizzy/OneDrive - King's College London/PhD/PhD_paperdrafts/all_genes_meta.csv")
gene_info <- filter(gene_info, species == "MOU")
names(colData(ms.sce))
rownames(ms.sce)
rowData(ms.sce)$gene_name <- gene_info$gene_name[match(rownames(ms.sce), gene_info$gene_id)]
rowData(ms.sce)

# functions from Sepp et al 2023 
# available at https://gitlab.com/kaessmannlab/mammalian-cerebellum
sce.get.gene.exprs <- function(sce, gene, umi=F){
  out <- NULL
  gene.orig <- gene
  
  if(!is.null(ncol(rowData(sce)))){
    if("gene_name" %in% colnames(rowData(sce))){
      if(gene %in% rowData(sce)[,"gene_name"]){
        gene <- rownames(rowData(sce))[which(rowData(sce)[,"gene_name"] == gene)]
      }
    }
  }
  
  if(gene %in% rownames(sce)){
    if(!umi){
      out <- log1p(assay(sce, "umi")[gene,] / colSums(assay(sce, "umi")))
    }else{
      out <- assay(sce, "umi")[gene,]
    }
    
  }
  
  return(out)
}

#calculate cell numbers
cell_age_num <- as.data.frame(table(colData(ms.sce)$Stage, colData(ms.sce)$cell_type))
colnames(cell_age_num)<- c("Stage", "cell_type", "Freq")
cell_num <- as.data.frame(table(colData(ms.sce)$cell_type)) %>%
  arrange(., desc(Freq))
head(cell_num$Var1, 10)

keep <- filter(cell_age_num, Freq > 9) #keep clusters with at least 10 cells
genes <- c("Foxp2", "Foxp1", "Foxp4")
tmp <- lapply(genes, function(g){
  exprs <- sce.get.gene.exprs(ms.sce, g, umi = T)
  if(is.null(exprs)){
    return(NULL)
  }else{
    colData(ms.sce) %>%
      as_tibble() %>%
      add_column(umi = exprs) %>%
      inner_join(., keep) %>%
      mutate(exprs = umi / size_factor * 1e6) %>%
      mutate(expresses = ifelse(umi > 0 , 1, 0)) %>%
      group_by(cell_type, Stage, stage.ord)%>%
      summarise(exprs = mean(exprs),
                expresses = sum(expresses) / n()) %>%
      ungroup() %>%
      add_column(gene = g) %>%
      mutate(max.exp = max(exprs)) %>%
      mutate(exprs = ((exprs ) / max(exprs )))
  }}) %>% do.call(rbind, .)

stages <- c("E10.5", "E11.5", "E12.5", "E13.5", "E14.5", "E15.5", "E17.5", "P0", "P4", "P7", "P14", "Adult")
p1 <- tmp %>% 
  ggplot(aes(x = factor(stage.ord, labels = stages), 
             y = factor(cell_type, 
                        levels = rev(sort(unique(cell_type)))))) + 
  geom_point(aes(color = exprs, size = expresses)) +
  theme_classic() +
  scale_colour_gradientn(name = "Max scaled expression",
                         colours = colorRampPalette(brewer.pal(n=9, name = "Purples"))(100)) +
  labs(x = "Stage",
       y = "Cell type",
       color = "Max. scaled\nexpression",
       size = "Proportion of\npositive cells") + facet_wrap(~gene)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = "italic"))
p1

ms.pc <- colData(ms.sce) %>% 
  as.data.frame() %>% 
  rownames_to_column('cell_id') %>% 
  as_tibble() %>%
  filter(., cell_type == "Purkinje") %>%
  filter(., !(is.na(subtype)))

ms.PC.data <- ms.sce[,ms.pc$cell_id]

genes <- c("Foxp2", "Foxp1", "Foxp4", "Ebf1", "Ebf2", "Cdh9", "Etv1", "Rorb")
tmp <- lapply(genes, function(g){
  exprs <- sce.get.gene.exprs(ms.PC.data, g, umi = T)
  if(is.null(exprs)){
    return(NULL)
  }else{
    colData(ms.PC.data) %>%
      as_tibble() %>%
      add_column(umi = exprs) %>%
      mutate(exprs = umi / size_factor * 1e6) %>%
      mutate(expresses = ifelse(umi > 0 , 1, 0)) %>%
      group_by(subtype) %>%
      summarise(exprs = mean(exprs),
                expresses = sum(expresses) / n()) %>%
      add_column(gene = g) %>%
      mutate(exprs = ((exprs ) / max(exprs )))
  }
}) %>% do.call(rbind, .)
p2 <- tmp %>% 
  ggplot(aes(x = factor(subtype, labels = c("Cdh9", "Etv1", "Foxp1", "Rorb")),
             y = factor(gene, levels = rev(genes)))) + 
  geom_point(aes(color = exprs, size = expresses)) +
  theme_classic() +
  scale_colour_gradientn(name = "Max scaled expression",
                         colours = colorRampPalette(brewer.pal(n=9, name = "Purples"))(100)) +
  #scale_colour_viridis_c(option = "B", direction = -1)  +
  labs(
    y = "Gene",
    x = "Purkinje cell subtype",
    color = "Max. scaled\nexpression",
    size = "Proportion of\npositive cells")+
  theme(axis.text = element_text(colour = "black", size = 12),
        axis.text.y = element_text(face = "italic"),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))+
  guides(colour = guide_colorbar(order = 1),
         size = guide_legend(order = 2))
p2
