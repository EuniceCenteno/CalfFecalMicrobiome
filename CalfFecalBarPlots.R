#             R-Script for Calf Fecal microbiota Taxonomy Barplots
#                     Author: Ruth Eunice Centeno Martinez, 2020
#          Written at Dr. Tim Johnson Lab, Dept. of Animal Sciences,       #
#                             Purdue University, 2020
#                               rcenteno@purdue.edu
library(ggplot2)
library(tidyr) #separate function
library(reshape2) #melt function
library(dplyr)
library(naniar) # for replace_with_na_all function
library(data.table)

setwd("~/Desktop/eunice/CalfFecal/NewQiime/TrimmedFecal/exported/")
##

##### Taxonomy barplot
#The OTU table as exported from qiime has a pound sign before the header row. You need to delete that pound sign in a text editor.
OTU = read.csv("rarified-table.tsv", header=TRUE, sep="\t", stringsAsFactors = F)
str(OTU)
row.names(OTU) = OTU[,1]
OTU = OTU[,-1]
OTU2 <- transpose(OTU)
rownames(OTU2) <- colnames(OTU)
colnames(OTU2) <- rownames(OTU)
OTU.clean <- as.data.frame(OTU2)
OTU.clean

#Importing the metadata file
meta = read.table("CalfFecal_meta.txt", header=TRUE, row.names=1, sep="\t")
meta$d <- factor(meta$d)
meta$calf <- factor(meta$calf)
meta$batch <- factor(meta$batch)
meta$block <- factor(meta$block)
meta$treatment <- factor(meta$treatment)
meta$sample <- rownames(meta)
order_groups <- meta$sample
meta$sample<- factor(meta$sample)
levels(meta$treatment) <- list("CON"="Pre-CON", "SCFP"="Pre-SCFP", "CON"="CON","SCFP"="SCFP")

### CALCULATION OF THE ABUNDANCE OF EACH OTU  
otu.summary <- prop.table(as.matrix(OTU.clean), 1) 
str(otu.summary)
otu_abund <- colSums(otu.summary)
otu_abund2 <- as.data.frame(otu_abund)
otu.summary <- rbind(otu_abund, otu.summary)
str(otu.summary)
otu.summary_sorted <- otu.summary[,order(otu.summary[1,], decreasing = TRUE)]
str(otu.summary_sorted)
melt_otu <- reshape2::melt(otu.summary_sorted[, c(1:6212)]) ###TOTAL NUMBER OF OTUS
str(melt_otu)
colnames(melt_otu) <- c("Sample", "OTU", "Abundance")
str(melt_otu)
levels(melt_otu$Sample)

#merging the abundance of each OTU with the metadata and the taxonomy file
taxonomy = read.table("taxaCalfFecal.txt", header=TRUE, row.names=1, sep="\t")
meta_otu <- merge(meta, melt_otu, by.x = 0, by.y = "Sample")
meta_otu <- meta_otu[-c(2:3, 8:9)]
meta_otu_tax <- merge(meta_otu, taxonomy, by.x = "OTU", by.y = 0)
str(meta_otu_tax)
levels(meta_otu_tax$d)
meta_otu_tax$Row.names <- factor(meta_otu_tax$Row.names, levels = order_groups)
summary(meta_otu_tax$Row.names) ###to check that all the samples have the same number of OTUs (6199 total, same value from the taxonomy file) 
meta_otu_tax$treatment <- factor(meta_otu_tax$treatment)
str(meta_otu_tax)
meta_otu_tax$d <- factor(meta_otu_tax$d)
meta_otu_tax$Family <- factor(meta_otu_tax$Family)
meta_otu_tax$Genus <- factor(meta_otu_tax$Genus)
meta_otu_tax$Phylum <- factor(meta_otu_tax$Phylum)
meta_otu_tax$OTU <- factor(meta_otu_tax$OTU)
str(meta_otu_tax)


## PHYLUM LEVEL
num_genera <- 1500 # we need 100 OTUs in order to get the 25 most abundant Genus

melt_otu1 <- reshape2::melt(otu.summary_sorted[, c(1:num_genera)])
colnames(melt_otu1) <- c("Sample", "OTU", "Abundance")
tail(melt_otu1)

#Putting it all together: merge melt_otu, metadata, taxonomy tables
meta_otu1 <- merge(meta, melt_otu1, by.x = 0, by.y = "Sample")
meta_otu1 <- meta_otu1[-c(2:3, 8:9)]
meta_otu_tax1 <- merge(meta_otu1, taxonomy, by.x = "OTU", by.y = 0)
str(meta_otu_tax1)
levels(meta_otu_tax1$d)
meta_otu_tax1$Row.names <- factor(meta_otu_tax1$Row.names, levels = order_groups)
summary(meta_otu_tax1$Row.names) ###to check that all the samples have the same number of OTUs (346 total) 
meta_otu_tax1$treatment <- factor(meta_otu_tax1$treatment)
str(meta_otu_tax1)
meta_otu_tax1$d <- factor(meta_otu_tax1$d)
meta_otu_tax1$Family <- factor(meta_otu_tax1$Family)
meta_otu_tax1$Phylum <- factor(meta_otu_tax1$Phylum)
meta_otu_tax1$Genus <- factor(meta_otu_tax1$Genus)
meta_otu_tax1$Family <- factor(meta_otu_tax1$Family)
levels(meta_otu_tax1$Phylum)


##Whole phylum abuundance 
Phylum <- meta_otu_tax %>% 
  group_by(Row.names, Phylum) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(Phylum) %>%
  summarise(taxa.average = mean(taxa.sum))
attach(Phylum)
Phylum <- Phylum[order(-taxa.average),]### relative abundance 

### Calculation of the Phylum relative abundance for each time and treatment
PhylumAB <- meta_otu_tax1 %>% 
  group_by(d, Row.names, treatment, Phylum) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(d, treatment, Phylum) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance
str(PhylumAB)
PhylumAB$Phylum <- factor(PhylumAB$Phylum)
levels(PhylumAB$Phylum)

my_colors <- c(
  '#a6cee3','#1f78b4','#b3df8a','#33a03c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab3d6','#6a3d9a','#ffff99','#b15938', 
  "#CBD588", "#5F7FC7", "orange","#DA5734", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14385", "#653936", "#C84348", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black"
)

#Plot the graph 
ggplot(PhylumAB, aes(x = treatment, y = taxa.average, fill =Phylum)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  facet_wrap(vars(d), scales = "free") +
  ylim(c(0,1)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  theme(strip.text = element_text(size = 13, face = "bold")) +
  theme(legend.text = element_text(size=13)) +
  theme(legend.text = element_text(margin = margin(t = 8))) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=13, face="bold"), axis.title.y = element_text(color="black", size=13, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 13), axis.text.y = element_text(color = "black", size = 13)) +
  theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Relative Abundance Phylum (Top 15)")) +  labs(x='Treatment')

### FAMILY LEVEL
num_genera <- 120 # we need 100 OTUs in order to get the 25 most abundant Genus

melt_otu3 <- reshape2::melt(otu.summary_sorted[, c(1:num_genera)])
colnames(melt_otu3) <- c("Sample", "OTU", "Abundance")
tail(melt_otu3)

#Putting it all together: merge melt_otu, metadata, taxonomy tables
meta_otu3 <- merge(meta, melt_otu3, by.x = 0, by.y = "Sample")
meta_otu3 <- meta_otu3[-c(2:3, 8:9)]
meta_otu_tax3 <- merge(meta_otu3, taxonomy, by.x = "OTU", by.y = 0)
str(meta_otu_tax3)
levels(meta_otu_tax3$d)
meta_otu_tax3$Row.names <- factor(meta_otu_tax3$Row.names, levels = order_groups)
summary(meta_otu_tax3$Row.names) ###to check that all the samples have the same number of OTUs (120) 
meta_otu_tax3$treatment <- factor(meta_otu_tax3$treatment)
str(meta_otu_tax3)
meta_otu_tax3$d <- factor(meta_otu_tax3$d)
meta_otu_tax3$Family <- factor(meta_otu_tax3$Family)
meta_otu_tax3$Genus <- factor(meta_otu_tax3$Genus)
meta_otu_tax3$Family <- factor(meta_otu_tax3$Family)
levels(meta_otu_tax3$Family)

##Whole family abuundance 
Family <- meta_otu_tax %>% 
  group_by(Row.names, Family) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(Family) %>%
  summarise(taxa.average = mean(taxa.sum))
attach(Family)
Family <- Family[order(-taxa.average),]### relative abundance


FamilyAB <- meta_otu_tax3 %>% 
  group_by(d, Row.names, treatment, Family) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(d, treatment, Family) %>%
  summarise(taxa.average = mean(taxa.sum)) 
str(FamilyAB)
FamilyAB$Family <- factor(FamilyAB$Family)
levels(FamilyAB$Family)

#Plot the graph 
ggplot(FamilyAB, aes(x = treatment, y = taxa.average, fill =Family)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  facet_wrap(vars(d), scales = "free") +
  ylim(c(0,1)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  #theme(legend.position="bottom") +
  #theme(axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
  theme(strip.text = element_text(size = 13, face = "bold")) +
  theme(legend.text = element_text(size=13)) +
  theme(legend.text = element_text(margin = margin(t = 8))) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=13, face="bold"), axis.title.y = element_text(color="black", size=13, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 13), axis.text.y = element_text(color = "black", size = 13)) +
  theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Average Relative Abundance Family (Top 25)")) +  labs(x='Treatment')


###----- TAXONOMY CALCULATION FOR TIME  
## The abundance at a family level grouped by time
num_genera <- 53 # we need 100 OTUs in order to get the 25 most abundant Genus

melt_otu4 <- reshape2::melt(otu.summary_sorted[, c(1:num_genera)])
colnames(melt_otu4) <- c("Sample", "OTU", "Abundance")
tail(melt_otu4)


#Putting it all together: merge melt_otu, metadata, taxonomy tables
meta_otu4 <- merge(meta, melt_otu4, by.x = 0, by.y = "Sample")
meta_otu4 <- meta_otu4[-c(3:4, 8:9)]
meta_otu_tax4 <- merge(meta_otu4, taxonomy, by.x = "OTU", by.y = 0)
str(meta_otu_tax4)
levels(meta_otu_tax4$d)
meta_otu_tax4$Row.names <- factor(meta_otu_tax4$Row.names, levels = order_groups)
summary(meta_otu_tax4$Row.names) ###to check that all the samples have the same number of OTUs (53) 
meta_otu_tax4$treatment <- factor(meta_otu_tax4$treatment)
str(meta_otu_tax4)
meta_otu_tax4$d <- factor(meta_otu_tax4$d)
meta_otu_tax4$Family <- factor(meta_otu_tax4$Family)
meta_otu_tax4$Genus <- factor(meta_otu_tax4$Genus)
meta_otu_tax4$Phylum <- factor(meta_otu_tax4$Phylum)
levels(meta_otu_tax4$Genus)

##Whole family abuundance 
Genus <- meta_otu_tax %>% 
  group_by(Row.names, Genus) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(Genus) %>%
  summarise(taxa.average = mean(taxa.sum))
attach(Genus)
Genus <- Genus[order(-taxa.average),]### relative abundance


GenusAB <- meta_otu_tax4 %>% 
  group_by(d, sample, treatment, Genus) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(d, treatment, Genus) %>%
  summarise(taxa.average = mean(taxa.sum)) 
str(PhylumAB)
GenusAB$Genus <- factor(GenusAB$Genus)
levels(GenusAB$Genus) #25 Genus total

# PLOT FOR THE FIRST 25 GENUS
my_colors <- c(
  '#a6cee3','#1f78b4','#b3df8a','#33a03c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab3d6','#6a3d9a','#ffff99','#b15938', 
  "#CBD588", "#5F7FC7", "orange","#DA5734", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14385", "#653936", "#C84348", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black"
)

ggplot(GenusAB, aes(x = treatment, y = taxa.average, fill =Genus)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  facet_wrap(vars(d), scales = "free") +
  ylim(c(0,1)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  #theme(legend.text=element_text(size=8)) +
  #theme(legend.position="bottom") +
  #theme(axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
  theme(strip.text = element_text(size = 13, face = "bold")) +
  theme(legend.text = element_text(size=13)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  #theme(legend.key.size = unit(1.5, "cm")) +
  theme(legend.text = element_text(margin = margin(t = 8))) +
  theme(axis.title.x = element_text(color="black", size=13, face="bold"), axis.title.y = element_text(color="black", size=13, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 13), axis.text.y = element_text(color = "black", size = 13)) +
  theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Average Relative Abundance Genus (Top 25)")) +  labs(x='Treatment')


