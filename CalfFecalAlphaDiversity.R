#             R-Script for Calf Fecal microbiota Taxonomy Barplots
#                     Author: Ruth Eunice Centeno Martinez, 2020
#          Written at Dr. Tim Johnson Lab, Dept. of Animal Sciences,       #
#                             Purdue University, 2020
#                               rcenteno@purdue.edu

#install.packages("car")
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

setwd("~/Desktop/eunice/CalfFecal/CalfFecal01_2021/NewQiime/TrimmedFecal/exported/") #sets new working directory for Windows systems (remember to replace … with your filepath)
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
  ylab("Observed ASVs") +xlab ("d") +
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
  ylab("Observed ASVs") +xlab ("d") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=16, face="bold"), axis.title.y = element_text(color="black", size=16, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 16), axis.text.y = element_text(color = "black", size = 16)) +
  theme(text=element_text(family="Times New Roman"))

Ev <- ggplot(meta, aes(x=d, y=pielou_e, fill=treatment)) + 
  geom_boxplot() + 
  scale_fill_manual(values = my_colors2) +
  #guides(fill=FALSE) +
  labs(fill= "Treatment") +
  theme_bw()+
  ylab("Evenness (Pielou)") +xlab ("d") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=16)) +
  theme(legend.title = element_text(size = 16, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=16, face="bold"), axis.title.y = element_text(color="black", size=16, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 16), axis.text.y = element_text(color = "black", size = 16)) +
  theme(text=element_text(family="Times New Roman"))

ggarrange(CH, Ev, labels = c("A", "B"),
          ncol = 2, font.label = list(family = "Times New Roman"))

##LPS figure
levels(meta$LPS)
LPS <- subset(meta, subset=LPS %in% c("No", "Yes")) #Samples on day 56
str(LPS)

#Treatment effect
LPSml <- lm(observed_otus ~ treatment, data = LPS)
#Running ANOVA
options(contrasts=c("contr.sum", "contr.poly"))

#Anova type III SS
ANOV_LPSml <- Anova(LPSml, type = 3)
ANOV_LPSml

mean_LPSml<-emmeans(LPSml, ~ treatment)
mean_LPSml <- data.frame(mean_LPSml)
mean_LPSml
day <- "56 d"
mean_LPSml$day <- day



LPSml2 <- lm(pielou_e ~ treatment, data = LPS)
#Running ANOVA
options(contrasts=c("contr.sum", "contr.poly"))

#Anova type III SS
ANOV_LPSml2 <- Anova(LPSml2, type = 3)
ANOV_LPSml2

mean_LPSml2<-emmeans(LPSml2, ~ treatment)
mean_LPSml2 <- data.frame(mean_LPSml2)
mean_LPSml2
day <- "56 d"
mean_LPSml2$day <- day

lps <- ggplot(mean_LPSml, aes(x=treatment, y=emmean, fill=treatment)) + 
  geom_bar(stat = "identity", alpha=0.8, colour="black") +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2, position = position_dodge(width = 1)) + 
  scale_fill_manual(values = my_colors) +
  guides(fill=FALSE) +
  labs(fill= "Treatment") +
  facet_grid(.~day) +
  theme_bw()+
  ylab("Observed ASVs") +xlab ("Treatment") +
  theme(strip.text = element_text(size = 16, face = "bold")) +
  theme(legend.text = element_text(size=16)) +
  theme(legend.title = element_text(size = 16, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=16, face="bold"), axis.title.y = element_text(color="black", size=16, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 16), axis.text.y = element_text(color = "black", size = 16)) +
  theme(text=element_text(family="Times New Roman"))

lps2 <- ggplot(mean_LPSml2, aes(x=treatment, y=emmean, fill=treatment)) + 
  geom_bar(stat = "identity", alpha=0.8, colour="black") +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2, position = position_dodge(width = 1)) + 
  scale_fill_manual(values = my_colors) +
  guides(fill=FALSE) +
  labs(fill= "Treatment") +
  facet_grid(.~day) +
  theme_bw()+
  ylab("Evenness (Pielou)") +xlab ("Treatment") +
  theme(strip.text = element_text(size = 16, face = "bold")) +
  theme(legend.text = element_text(size=16)) +
  theme(legend.title = element_text(size = 16, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=16, face="bold"), axis.title.y = element_text(color="black", size=16, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 16), axis.text.y = element_text(color = "black", size = 16)) +
  theme(text=element_text(family="Times New Roman"))


ggarrange(lps, lps2, labels = c("A", "B"),
          ncol = 2, font.label = list(family = "Times New Roman"))

#LPS effect
str(LPS)
LPS %>% tally()
LPS %>% count(LPS)
#LPS  n
#1 Yes 27
#2  No 27
levels(LPS$LPS)
LPS$LPS <- factor(LPS$LPS)
LPS$sample <- factor(LPS$sample)
LPSe <- lm(observed_otus ~ LPS, data = LPS)
#Running ANOVA
ANOV_LPSe <- Anova(LPSe)
ANOV_LPSe

mean_LPSe<-emmeans(LPSe, ~ LPS)
mean_LPSe <- data.frame(mean_LPSe)
mean_LPSe
day <- "56 d"
mean_LPSe$day <- day

LPSp <- lm(pielou_e ~ LPS, data = LPS)
#Running ANOVA
ANOV_LPSp <- Anova(LPSp)
ANOV_LPSp

mean_LPSp<-emmeans(LPSp, ~ LPS)
mean_LPSp <- data.frame(mean_LPSp)
mean_LPSp
day <- "56 d"
mean_LPSp$day <- day

LPSCh <- ggplot(mean_LPSe, aes(x=LPS, y=emmean, fill=LPS)) + 
  geom_bar(stat = "identity", alpha=0.8, colour="black") +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2, position = position_dodge(width = 1)) + 
  scale_fill_manual(values = my_colors) +
  guides(fill=FALSE) +
  labs(fill= "LPS") +
  facet_grid(.~day) +
  theme_bw()+
  ylab("Observed ASVs") +xlab ("LPS") +
  theme(strip.text = element_text(size = 16, face = "bold")) +
  theme(legend.text = element_text(size=16)) +
  theme(legend.title = element_text(size = 16, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=16, face="bold"), axis.title.y = element_text(color="black", size=16, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 16), axis.text.y = element_text(color = "black", size = 16)) +
  theme(text=element_text(family="Times New Roman"))

LPSPi <- ggplot(mean_LPSp, aes(x=LPS, y=emmean, fill=LPS)) + 
  geom_bar(stat = "identity", alpha=0.8, colour="black") +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2, position = position_dodge(width = 1)) + 
  scale_fill_manual(values = my_colors) +
  guides(fill=FALSE) +
  labs(fill= "LPS") +
  facet_grid(.~day) +
  theme_bw()+
  ylab("Evenness (Pielou)") +xlab ("LPS") +
  theme(strip.text = element_text(size = 16, face = "bold")) +
  theme(legend.text = element_text(size=16)) +
  theme(legend.title = element_text(size = 16, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=16, face="bold"), axis.title.y = element_text(color="black", size=16, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 16), axis.text.y = element_text(color = "black", size = 16)) +
  theme(text=element_text(family="Times New Roman"))

ggarrange(LPSCh, LPSPi, labels = c("A", "B"),
          ncol = 2, font.label = list(family = "Times New Roman"))
