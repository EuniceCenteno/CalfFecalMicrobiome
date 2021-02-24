#             R-Script for Calf Fecal microbiota Taxonomy Barplots
#                     Author: Ruth Eunice Centeno Martinez, 2020
#          Written at Dr. Tim Johnson Lab, Dept. of Animal Sciences,       #
#                             Purdue University, 2020
#                               rcenteno@purdue.edu

#install.packages("ggpubr")
library(tidyverse)
library(ggpubr)
library(rstatix)
library(afex)
library(emmeans)
library(lme4)
library(lattice)
library("latticeExtra")
library(dplyr)
library(car)

setwd("~/Desktop/eunice/CalfFecal/NewQiime/TrimmedFecal/exported/") #sets new working directory for Windows systems (remember to replace … with your filepath)
# Importing the metadata
meta <- read.table("CalfFecal_meta.txt", header=TRUE, row.names=2, sep="\t")

# All the alpha diversity metrics were obtained from Qiime2
otu_table <- read.table("observed_otus.tsv", header=TRUE, row.names=1, sep="\t")
shannon <-read.table("shannon.tsv", header=TRUE, row.names=1, sep="\t")
evenness <- read.table("evenness.tsv", header=TRUE, row.names=1, sep="\t")
alpha_diversity <- merge(otu_table, shannon, by.x = 0, by.y = 0)
alpha_diversity <- merge(alpha_diversity, evenness, by.x = "Row.names", by.y = 0)
meta <- merge(meta, alpha_diversity, by.x = 0, by.y = "Row.names")
row.names(meta) <- meta$Row.names
meta <- meta[,-1]
meta$sample <- rownames(meta)

## Checking data
str(meta)
## we need to make calf (this is the ID variable) and day a factor 
meta$d <- factor(meta$d)
meta$calf <- factor(meta$calf)
meta$batch <- factor(meta$batch)
meta$block <- factor(meta$block)
meta$treatment <- factor(meta$treatment)
meta$LPS <- factor(meta$LPS)
levels(meta$LPS)
levels(meta$LPS) <- list("X"="No", "Yes"="Yes_56","No"="N_56", "P-LPS84"="P-LPS84","No_84"="No_84", "P-LPS112"="P-LPS112","No_112"="No_112")
str(meta)
levels(meta$d)
levels(meta$d) <- list("0"="0d", "28"="28d", "56"="56d", "84"="84d", "112"="112d")
levels(meta$treatment)
levels(meta$treatment) <- list("CON"="Pre-CON", "SCFP"="Pre-SCFP", "CON"="CON","SCFP"="SCFP")

# data visualization

my_colors <- c("black", "gray35", "gray48", "gray83", "white")

Ch <- ggplot(meta, aes(x=d, y=observed_otus, fill=d)) + 
  geom_boxplot() + 
  scale_fill_manual(values = my_colors) +
  guides(fill=FALSE) +
  labs(fill= "Treatment") +
  theme_bw()+
  ylab("Observed OTUs (Chao)") +xlab ("d") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)) +
  theme(text=element_text(family="Times New Roman"))

EV <- ggplot(meta, aes(x=d, y=pielou_e, fill=d)) + 
  geom_boxplot() + 
  scale_fill_manual(values = my_colors) +
  guides(fill=FALSE) +
  labs(fill= "Treatment") +
  theme_bw()+
  ylab("Evenness (Pielou)") +xlab ("d") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)) +
  theme(text=element_text(family="Times New Roman"))

ggarrange(Ch, EV, labels = c("A", "B"),
          ncol = 2, font.label = list(family = "Times New Roman"))

my_colors2 <- c("black", "gray")
CH <- ggplot(meta, aes(x=d, y=observed_otus, fill=treatment)) + 
  geom_boxplot() + 
  scale_fill_manual(values = my_colors2) +
  guides(fill=FALSE) +
  labs(fill= "Treatment") +
  theme_bw()+
  ylab("Observed OTUs (Chao)") +xlab ("d") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)) +
  theme(text=element_text(family="Times New Roman"))

Ev <- ggplot(meta, aes(x=d, y=pielou_e, fill=treatment)) + 
  geom_boxplot() + 
  scale_fill_manual(values = my_colors2) +
  #guides(fill=FALSE) +
  labs(fill= "Treatment") +
  theme_bw()+
  ylab("Evenness (Pielou)") +xlab ("d") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)) +
  theme(text=element_text(family="Times New Roman"))

ggarrange(CH, Ev, labels = c("A", "B"),
          ncol = 2, font.label = list(family = "Times New Roman"))

##LPS figure
levels(meta$LPS)
LPS <- subset(meta, subset=LPS %in% c("No", "Yes")) #Samples on day 56
str(LPS)

lps <- ggplot(LPS, aes(x=LPS, y=pielou_e, fill=LPS)) + 
  geom_boxplot() + 
  scale_fill_manual(values = my_colors) +
  guides(fill=FALSE) +
  #labs(fill= "Treatment") +
  theme_bw()+
  ylab("Evenness (Pielou)") +xlab ("LPS Challenge") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)) +
  theme(text=element_text(family="Times New Roman"))

lps2 <- ggplot(LPS, aes(x=LPS, y=observed_otus, fill=LPS)) + 
  geom_boxplot() + 
  scale_fill_manual(values = my_colors) +
  guides(fill=FALSE) +
  #labs(fill= "Treatment") +
  theme_bw()+
  ylab("Observed OTUs (Chao)") +xlab ("LPS Challenge") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)) +
  theme(text=element_text(family="Times New Roman"))


ggarrange(lps2, lps, labels = c("A", "B"),
          ncol = 2, font.label = list(family = "Times New Roman"))
