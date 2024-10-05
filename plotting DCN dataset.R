#Plotting human DCN RNAseq from Kebschull et al., 2020

library(Seurat)
library(tidyverse)
library(cowplot)

#data avaliable at 
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4873766
#integrated inh and ex human data following https://github.com/justuskebschull/CNcode_final

neurons2 <- readRDS("C:/Users/Lizzy/OneDrive - King's College London/PhD/PhD_paperdrafts/DCN_neurons2.rds")
FeaturePlot(neurons2, "FOXP2")
FeaturePlot(neurons2, "SLC17A6")
FeaturePlot(neurons2, "GAD1")
FeaturePlot(neurons2, "NR4A2")
FeaturePlot(neurons2, "ZFHX4")
FeaturePlot(neurons2, "MEIS2")