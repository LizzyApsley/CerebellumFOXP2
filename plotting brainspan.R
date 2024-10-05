#plotting brainspan bulk RNAseq data - FOXP2

#load programs
require(tidyverse) 
require(ggplot2) 
require(reshape2) 
library(RColorBrewer)

#read in and arrange data
dir <- "C:/Users/Lizzy/OneDrive - King's College London/brainspan/genes_matrix_csv"
#data available for download from https://www.brainspan.org/static/download.html

coldata <- read.csv(paste(dir, "columns_metadata.csv", sep = "/"))
rowdata <- read.csv(paste(dir, "rows_metadata.csv", sep = "/"))                    
expdata <- read.csv(paste(dir, "expression_matrix.csv", sep = "/"), header= F)

geneposition <- rowdata[rowdata$gene_symbol == "FOXP2",]$row_num
exp.FOXP2 <- expdata[geneposition,]
data <- cbind(coldata, t(exp.FOXP2[,-1]))
colnames(data)[9] <- "FOXP2"
data <- mutate(data, exp = log2(FOXP2+0.5)) %>%
  select(., age, donor_id, structure_id, structure_acronym, exp) %>%
  mutate(., sample.info = paste(age, structure_acronym))

ages.key <- data.frame(
  age = unique(data$age),
  ref_age = seq(1:31))
ages.key <- arrange(ages.key, ref_age)
data <- inner_join(data, ages.key)
age.labels <- unique(data$age)

summary.data <- data %>%
  group_by(., structure_acronym, ref_age) %>%
  summarise(mean.exp = mean (exp))

max.exp <- group_by(data, structure_acronym) %>%
  summarise(max.exp = max(exp)) %>%
  arrange(., desc(max.exp))

colour.scheme <- cbind(max.exp, 
                       c(brewer.pal(8, "Dark2"), 
                         rep("grey", length(max.exp$structure_acronym)-8)))
colnames(colour.scheme)[3] <- "colour"

ggplot() + 
  geom_point(data = data, aes (x = factor(ref_age, labels = age.labels), y = exp, 
                               colour = structure_acronym), alpha = 0.2)+
  geom_point(data = summary.data, aes (x = factor(ref_age, labels = age.labels), 
                                       y = mean.exp, colour = structure_acronym))+
  geom_line(data = summary.data, aes (x = factor(ref_age, labels = age.labels), 
                                      y = mean.exp, group = structure_acronym, 
                                      colour = structure_acronym), alpha = 0.6)+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, colour = "black", size = 10, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 10))+
  xlab("Age") + ylab("Expression (log2 RPKM)") +
  scale_color_manual(breaks = colour.scheme$structure_acronym,
                     values= colour.scheme$colour) + 
  labs(colour = "Brain region")
