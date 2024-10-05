#plotting GSEA results

#libraries
library(reshape2)
library(tidyverse)

#directories
resultsdir <- "C:/Users/Lizzy/OneDrive - King's College London/PhD/PhD_paperdrafts/rnaseqanalysis"


#Cerebellar regions (Aldinger et al., 2021 - LCM RNAseq)
LCM_pos <- read.table(file = paste(resultsdir, "GSEAwithmyowngenesets/599396-GSEA-AldLCM/599396/gsea_report_for_pos_1724189826370.tsv", sep = "/" ), 
                      sep = "\t", header = T)
LCM_neg <- read.table(file = paste(resultsdir, "GSEAwithmyowngenesets/599396-GSEA-AldLCM/599396/gsea_report_for_neg_1724189826370.tsv", sep = "/" ), 
                      sep = "\t", header = T)
LCM <- rbind(LCM_neg, LCM_pos) %>%
  mutate(., sig = ifelse(FDR.q.val <0.05, ifelse(NES > 0, "pos", "neg"), "nonsig"))

ggplot(LCM, aes(x = NES, y = factor(NAME, levels = c("RL", "EGL", "PC")),
                fill = sig))+
  geom_col(colour = "black")+
  theme_classic()+
  ylab("Cerebellar region") + xlab ("Normalised enrichment score")+
  theme(axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(colour = "black", size = 12))+
  scale_fill_manual(values = c( "#666666","#66A61E"), 
                    breaks = c("neg", "pos"))

#sepp cell type markers 
sepp_pos <- read.table(file = paste(resultsdir, "GSEAwithmyowngenesets/600165-GSEA-Seppgenesets/600165/gsea_report_for_pos_1724492674271.tsv", sep = "/" ), 
                       sep = "\t", header = T)
sepp_neg <- read.table(file = paste(resultsdir, "GSEAwithmyowngenesets/600165-GSEA-Seppgenesets/600165/gsea_report_for_neg_1724492674271.tsv", sep = "/" ), 
                       sep = "\t", header = T)
sepp <- rbind(sepp_neg, sepp_pos) %>%
  arrange(., NES) %>%
  mutate(., sig = ifelse(FDR.q.val <0.05, ifelse(NES > 0, "pos", "neg"), "nonsig"))

ggplot(sepp, aes(x = NES, y = factor(NAME, levels = NAME),
                 fill = sig))+
  geom_col(colour = "black")+
  theme_classic()+
  ylab("Cell type") + xlab ("Normalised enrichment score")+
  theme(axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(colour = "black", size = 12))+
  scale_fill_manual(values = c( "#666666","#66A61E", "white"), 
                    breaks = c("neg", "pos", "nonsig"))


#disease risk gene sets
disease_pos <- read.table(file = paste(resultsdir, "GSEAwithmyowngenesets/599391-GSEA-disordergenesets/599391/gsea_report_for_pos_1724188605955.tsv", sep = "/" ), 
                          sep = "\t", header = T)
disease_neg <- read.table(file = paste(resultsdir, "GSEAwithmyowngenesets/599391-GSEA-disordergenesets/599391/gsea_report_for_neg_1724188605955.tsv", sep = "/" ), 
                          sep = "\t", header = T)
disease <- rbind(disease_neg, disease_pos) %>%
  arrange(., NES) %>%
  mutate(., sig = ifelse(FDR.q.val <0.05, ifelse(NES > 0, "pos", "neg"), "nonsig"))

ggplot(disease, aes(x = NES, y = factor(NAME, levels = NAME),
                    fill = sig))+
  geom_col(colour = "black")+
  theme_classic()+
  ylab("Disorder") + xlab ("Normalised enrichment score")+
  theme(axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(colour = "black", size = 12))+
  scale_fill_manual(values = c( "#666666","#66A61E", "white"), 
                    breaks = c("neg", "pos", "nonsig"))

#disease DEG gene sets
diseaseDEG_pos <- read.table(file = paste(resultsdir, "GSEAwithmyowngenesets/599395-GSEA-disorderDEG/599395/gsea_report_for_pos_1724189579852.tsv", sep = "/" ), 
                             sep = "\t", header = T)
diseaseDEG_neg <- read.table(file = paste(resultsdir, "GSEAwithmyowngenesets/599395-GSEA-disorderDEG/599395/gsea_report_for_neg_1724189579852.tsv", sep = "/" ), 
                             sep = "\t", header = T)
diseaseDEG <- rbind(diseaseDEG_neg, diseaseDEG_pos) %>%
  arrange(., NAME) %>%
  mutate(., sig = ifelse(FDR.q.val <0.05, ifelse(NES > 0, "pos", "neg"), "nonsig"))

ggplot(diseaseDEG, aes(x = NES, y = factor(NAME, levels = rev(NAME)),
                       fill = sig))+
  geom_col(colour = "black")+
  theme_classic()+
  ylab("Disorder") + xlab ("Normalised enrichment score")+
  theme(axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(colour = "black", size = 12))+
  scale_fill_manual(values = c( "#666666","#66A61E", "white"), 
                    breaks = c("neg", "pos", "nonsig"))

#GOBP pathways
GOBP_neg <- read.table(file = "C:/Users/Lizzy/OneDrive - Nexus365/PhD/mNeon sorted CbO RNAseq/GSEA/497669-GOBP/497669/gsea_report_for_neg_1678197364749.tsv", 
                       sep = "\t", header = T)
GOBP_neg <- filter(GOBP_neg, FDR.q.val <0.05) %>%
  arrange(., NES) %>%
  head(., 10)%>%
  mutate(., sig = ifelse(FDR.q.val <0.05, ifelse(NES > 0, "pos", "neg"), "nonsig")) %>%
  mutate(., temp = NAME) %>%
  separate(., temp, into = c("Gotype", "t", "Goname"), sep = c(4,5)) %>%
  mutate(., name = paste(substring(Goname, 0,1), tolower(substring(Goname, 2)), sep = "")) %>%
  mutate(., newname = str_replace_all(name, pattern = "_", replacement =" "))


GOBP_pos <- read.table(file = "C:/Users/Lizzy/OneDrive - Nexus365/PhD/mNeon sorted CbO RNAseq/GSEA/497669-GOBP/497669/gsea_report_for_pos_1678197364749.tsv", 
                       sep = "\t", header = T)
GOBP_pos <- filter(GOBP_pos, FDR.q.val <0.05) %>%
  arrange(., desc(NES)) %>%
  head(., 10)%>%
  mutate(., sig = ifelse(FDR.q.val <0.05, ifelse(NES > 0, "pos", "neg"), "nonsig")) %>%
  mutate(., temp = NAME) %>%
  separate(., temp, into = c("Gotype", "t", "Goname"), sep = c(4,5)) %>%
  mutate(., name = paste(substring(Goname, 0,1), tolower(substring(Goname, 2)), sep = "")) %>%
  mutate(., newname = str_replace_all(name, pattern = "_", replacement =" "))


ggplot(GOBP_neg, aes(x = NES, y = factor(newname, levels = rev(newname)),
                     fill = sig))+
  geom_col(colour = "black")+
  theme_classic()+
  ylab("GO biological process") + xlab ("Normalised enrichment score")+
  theme(axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(colour = "black", size = 12))+
  scale_fill_manual(values = c( "#666666","#66A61E", "white"), 
                    breaks = c("neg", "pos", "nonsig"))

ggplot(GOBP_pos, aes(x = NES, y = factor(newname, levels = rev(newname)),
                     fill = sig))+
  geom_col(colour = "black")+
  theme_classic()+
  ylab("GO biological process") + xlab ("Normalised enrichment score")+
  theme(axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(colour = "black", size = 12))+
  scale_fill_manual(values = c( "#666666","#66A61E", "white"), 
                    breaks = c("neg", "pos", "nonsig"))

