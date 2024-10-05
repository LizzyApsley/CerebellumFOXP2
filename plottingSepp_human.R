#plotting human developing cerebellar snRNAseq data from Sepp et al., 2023
#expression of FOXP genes, PC subtype markers and cell type markers 

#libraries
library(SingleCellExperiment)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

#data
sce <- readRDS("C:/Users/Lizzy/OneDrive - King's College London/PhD/PhD_paperdrafts/hum_sce_final.rds")
gene_info <- read.csv("C:/Users/Lizzy/OneDrive - King's College London/PhD/PhD_paperdrafts/all_genes_meta.csv")
gene_info <- filter(gene_info, species == "HUM")
colData(sce)
names(colData(sce))
rownames(sce)
rowData(sce)$gene_name <- gene_info$gene_name[match(rownames(sce), gene_info$gene_id)]
rowData(sce)

#calculating cell numbers
cell_age_num <- as.data.frame(table(colData(sce)$Stage, colData(sce)$cell_type))
colnames(cell_age_num)<- c("Stage", "cell_type", "Freq")
cell_num <- as.data.frame(table(colData(sce)$cell_type)) %>%
  arrange(., desc(Freq))
head(cell_num$Var1, 10)

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

#modify sce.dotplot.area this to split cell type (y) vs stage (x)
keep <- filter(cell_age_num, Freq > 9) #keep clusters with at least 10 cells
genes <- c("FOXP2", "FOXP1", "FOXP4")
tmp <- lapply(genes, function(g){
  exprs <- sce.get.gene.exprs(sce, g, umi = T)
  if(is.null(exprs)){
    return(NULL)
  }else{
    colData(sce) %>%
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
tmp <- filter(tmp, !(stage.ord %in% c("13_adult", "11_infant", "12_toddler")))

stages <- c("7 wpc", "8 wpc", "9 wpc", "11 wpc", "17 wpc", "20 wpc", "newborn")
p1 <- tmp %>% 
  ggplot(aes(x = factor(stage.ord, labels = stages), 
             y = factor(cell_type, 
                      levels = rev(sort(unique(cell_type)))))) + 
  geom_point(aes(color = exprs, size = expresses)) +
  theme_classic() +
  #scale_colour_viridis_c(option = "B", direction = -1)  +
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

#######################################################

#plot FOXPs and EBF1/2 in PC subtypes where available

#subset PC with subtypes only 
pc <- colData(sce) %>% 
  as.data.frame() %>% 
  rownames_to_column('cell_id') %>% 
  as_tibble() %>%
  filter(., cell_type == "Purkinje") %>%
  filter(., !(is.na(subtype)))
PC.data <- sce[,pc$cell_id]

genes <- c("EBF1", "EBF2", "FOXP1", "FOXP2", "FOXP4")

tmp <- lapply(genes, function(g){
  exprs <- sce.get.gene.exprs(PC.data, g, umi = T)
  if(is.null(exprs)){
    return(NULL)
  }else{
    colData(PC.data) %>%
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
  ggplot(aes(x = factor(subtype, labels = c("early-born", "late-born")),
             y = factor(gene, levels = rev(genes)))) + 
  geom_point(aes(color = exprs, size = expresses)) +
  theme_classic() +
  scale_colour_gradientn(name = "Max scaled expression",
                         colours = colorRampPalette(brewer.pal(n=9, name = "Purples"))(100)) +
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

####################################################

#plot cell type marker examples in Sepp data

#major cell types and Sepp markers
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

tmp <- lapply(genes, function(g){
  exprs <- sce.get.gene.exprs(sce, g, umi = T)
  if(is.null(exprs)){
    return(NULL)
  }else{
    colData(sce) %>%
      as_tibble() %>%
      add_column(umi = exprs) %>%
      filter(!(is.na(cell_type))) %>%
      mutate(exprs = umi / size_factor * 1e6) %>%
      mutate(expresses = ifelse(umi > 0 , 1, 0)) %>%
      group_by(cell_type) %>%
      summarise(exprs = mean(exprs),
                expresses = sum(expresses) / n()) %>%
      add_column(gene = g) %>%
      mutate(exprs = ((exprs ) / max(exprs )))
  }
}) %>% do.call(rbind, .)
p3 <- tmp %>% 
  ggplot(aes(x = cell_type,
             y = factor(gene, levels = rev(genes)))) + 
  geom_point(aes(color = exprs, size = expresses)) +
  theme_classic() +
  scale_colour_gradientn(name = "Max scaled expression",
                         colours = colorRampPalette(brewer.pal(n=9, name = "Purples"))(100)) +
  #scale_colour_viridis_c(option = "B", direction = -1)  +
  labs(
    y = "Gene",
    x = "Celltype",
    color = "Max. scaled\nexpression",
    size = "Proportion of\npositive cells")+
  theme(axis.text = element_text(colour = "black", size = 12),
        axis.text.y = element_text(face = "italic"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))+
  guides(colour = guide_colorbar(order = 1),
         size = guide_legend(order = 2))
p3
