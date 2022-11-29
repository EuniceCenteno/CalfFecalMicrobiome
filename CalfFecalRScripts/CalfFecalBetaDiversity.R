#                 R-Script for Calf Fecal microbiota Alpha Diversity
#                     Author: Ruth Eunice Centeno Martinez, 2020
#          Written at Dr. Tim Johnson Lab, Dept. of Animal Sciences,       #
#                             Purdue University, 2020

###---------------- BETA DIVERSITY---------------------------
## All the information related to PCoA (coordinates) and rarefied table coordinates were obtained from Qiime2
#install.packages("car")
library(vegan) 
library(ggplot2)
library(ggpubr)
library(dplyr) #tally function
library(data.table) #transpose function
library(car) ##for type III SS 

setwd("~/Desktop/eunice/MSc/CalfFecal/CalfFecal01_2021/NewQiime/TrimmedFecal/exported/") #sets new working directory for Windows systems (remember to replace … with your filepath)
calffecal <- read.table("CalfFecal_meta.txt", header=TRUE, row.names=1, sep="\t")
OTU = read.csv("rarified-table.tsv", header=TRUE, sep="\t", stringsAsFactors = F)
str(OTU)
row.names(OTU) = OTU[,1]
OTU = OTU[,-1]
OTU2 <- transpose(OTU)
rownames(OTU2) <- colnames(OTU)
colnames(OTU2) <- rownames(OTU)
OTU.clean <- as.data.frame(OTU2)
OTU.clean

str(calffecal)
calffecal$d <- factor(calffecal$d)
calffecal$calf <- factor(calffecal$calf)
calffecal$batch <- factor(calffecal$batch)
calffecal$block <- factor(calffecal$block)
calffecal$treatment <- factor(calffecal$treatment)
calffecal$LPS <- factor(calffecal$LPS)
str(calffecal)
levels(calffecal$LPS)
levels(calffecal$LPS) <- list("X"="No", "Yes"="Yes_56","No"="N_56", "P-LPS84"="P-LPS84","No_84"="No_84", "P-LPS112"="P-LPS112","No_112"="No_112")
levels(calffecal$d)
levels(calffecal$d) <- list("0"="0d", "28"="28d", "56"="56d", "84"="84d", "112"="112d")
levels(calffecal$treatment)
levels(calffecal$treatment) <- list("CON"="Pre-CON", "SCFP"="Pre-SCFP", "CON"="CON","SCFP"="SCFP")





##Plotting the Time effect
library(ellipse)
centroids <- aggregate(calffecal[,2:3], list(Group=calffecal$d), mean)
colnames(centroids) <- c('d','centroidX', 'centroidY')
conf.rgn  <- do.call(rbind,lapply(unique(calffecal$d),function(t)
  data.frame(d=as.character(t),
             ellipse(cov(calffecal[calffecal$d==t,2:3]),
                     centre=as.matrix(centroids[t,2:3]),
                     level=0.95),
             stringsAsFactors=FALSE)))
colnames(conf.rgn) <- c('d','ellipseX', 'ellipseY')

#merge the centroinds with the ellipse coordinates
TimeE <- merge(centroids, conf.rgn, by.x = "d", by.y= "d")

#my_colors <- c("gray32","black", "gray12", "gray48","gray0")
my_colors <- c(
 "#fb9a99","steelblue2",'palegreen3',"gray12","#D14385",'#e31a1c',
  '#fdbf6f','#ff7f00','#cab3d6','#6a3d9a','#ffff99','#b15938', 
  "#CBD588", "#5F7FC7", "orange","#DA5734", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14385", "#653936", "#C84348", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black"
)

ggplot(calffecal, aes(x=WUAxis1, y=WUAxis2)) +
  theme_bw()+
  geom_point(data = calffecal, aes(x=WUAxis1, y= WUAxis2,color=d), size=2.5) +
  stat_ellipse(geom = "polygon", data = TimeE, aes(x=ellipseX, y=ellipseY, group=d, fill= d), alpha = 0.1) + 
  scale_color_manual(values = c(my_colors)) +
  scale_fill_manual(values = c(my_colors)) +
  scale_shape_manual(values = c(1, 11, 15, 17, 4)) +
  labs(x='Axis 1 (44.8%)', y= 'Axis 2 (11.6%)') +
  labs(color= "d") +
  #facet_grid(id~Sickness) +
  guides(fill=FALSE) +
  #guides(color=FALSE) +
  #guides(color=FALSE) +
  theme(strip.text = element_text(size = 15, face = "bold")) +
  theme(legend.text = element_text(size=15)) +
  theme(legend.title = element_text(size = 15, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=15, face="bold"), axis.title.y = element_text(color="black", size=15, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 15), axis.text.y = element_text(color = "black", size = 15)) +
  theme(text=element_text(family="Times New Roman"))

#### Treatment comparison in each day
my_colors2 <- c("black", "gray")

day0 <- subset(calffecal,
               subset=d %in% c("0"))

A <- ggplot(day0, aes(x=WUAxis1, y=WUAxis2)) +
  theme_bw()+
  geom_point(data = day0, aes(x=WUAxis1, y=WUAxis2,color=treatment)) +
  scale_color_manual(values = my_colors2) +
  labs(x='Axis 1 (44.8%)', y= 'Axis 2 (11.6%)') +
  labs(color= "Treatment") +
  facet_grid(.~d) +
  guides(size=FALSE) +
  theme(strip.text = element_text(size = 14, face = "bold")) +
  theme(legend.text = element_text(size=14)) +
  theme(legend.title = element_text(size = 14, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold"), axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 14), axis.text.y = element_text(color = "black", size = 14)) +
  theme(text=element_text(family="Times New Roman"))

day28 <- subset(calffecal,
                subset=d %in% c("28"))

B <- ggplot(day28, aes(x=WUAxis1, y=WUAxis2)) +
  theme_bw()+
  geom_point(data = day28, aes(x=WUAxis1, y=WUAxis2,color=treatment)) +
  scale_color_manual(values = my_colors2) +
  labs(x='Axis 1 (44.8%)', y= 'Axis 2 (11.6%)') +
  labs(color= "Treatment") +
  facet_grid(.~d) +
  guides(size=FALSE) +
  theme(strip.text = element_text(size = 14, face = "bold")) +
  theme(legend.text = element_text(size=14)) +
  theme(legend.title = element_text(size = 14, face= "bold")) +
  theme(legend.key.size = unit(14, "point")) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold"), axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 14), axis.text.y = element_text(color = "black", size = 14)) +
  theme(text=element_text(family="Times New Roman"))

day56 <- subset(calffecal,
                subset=d %in% c("56"))

C <- ggplot(day56, aes(x=WUAxis1, y=WUAxis2)) +
  theme_bw()+
  geom_point(data = day56, aes(x=WUAxis1, y=WUAxis2,color=treatment)) +
  scale_color_manual(values = my_colors2) +
  labs(x='Axis 1 (44.8%)', y= 'Axis 2 (11.6%)') +
  labs(color= "Treatment") +
  facet_grid(.~d) +
  guides(size=FALSE) +
  theme(strip.text = element_text(size = 14, face = "bold")) +
  theme(legend.text = element_text(size=14)) +
  theme(legend.title = element_text(size = 14, face= "bold")) +
  theme(legend.key.size = unit(14, "point")) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold"), axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 14), axis.text.y = element_text(color = "black", size = 14)) +
  theme(text=element_text(family="Times New Roman"))

day84 <- subset(calffecal,
                subset=d %in% c("84"))

D <- ggplot(day84, aes(x=WUAxis1, y=WUAxis2)) +
  theme_bw()+
  geom_point(data = day84, aes(x=WUAxis1, y=WUAxis2,color=treatment)) +
  scale_color_manual(values = my_colors2) +
  labs(x='Axis 1 (44.8%)', y= 'Axis 2 (11.6%)') +
  labs(color= "Treatment") +
  facet_grid(.~d) +
  guides(size=FALSE) +
  theme(strip.text = element_text(size = 14, face = "bold")) +
  theme(legend.text = element_text(size=14)) +
  theme(legend.title = element_text(size = 14, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold"), axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 14), axis.text.y = element_text(color = "black", size = 14)) +
  theme(text=element_text(family="Times New Roman"))

day112 <- subset(calffecal,
                 subset=d %in% c("112"))

E <- ggplot(day112, aes(x=WUAxis1, y=WUAxis2)) +
  theme_bw()+
  geom_point(data = day112, aes(x=WUAxis1, y=WUAxis2,color=treatment)) +
  scale_color_manual(values = my_colors2) +
  labs(x='Axis 1 (44.8%)', y= 'Axis 2 (11.6%)') +
  labs(color= "Treatment") +
  facet_grid(.~d) +
  guides(size=FALSE) +
  theme(strip.text = element_text(size = 14, face = "bold")) +
  theme(legend.text = element_text(size=14)) +
  theme(legend.title = element_text(size = 14, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold"), axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 14), axis.text.y = element_text(color = "black", size = 14)) +
  theme(text=element_text(family="Times New Roman"))

ggarrange(A, B, C, D, E, labels = c("A", "B", "C", "D", "E"),
          ncol = 3, nrow=2, font.label = list(family = "Times New Roman"))

#LPS challenge
levels(calffecal$LPS)
LPS <- subset(calffecal, subset=LPS %in% c("No", "Yes")) #Samples on day 56
str(LPS)
LPS %>% tally()
LPS %>% count(LPS)
#LPS  n
#1 Yes 27
#2  No 27
levels(LPS$LPS)
LPS$LPS <- factor(LPS$LPS)
LPS$sample <- factor(LPS$sample)


my_colors <- c("black", "gray")
ggplot(LPS, aes(x=WUAxis1, y=WUAxis2)) +
  theme_bw()+
  geom_point(data = LPS, aes(x=WUAxis1, y=WUAxis2,color=LPS, shape=treatment), size=2) +
  scale_color_manual(values = my_colors) +
  scale_shape_manual(values = c(15, 11,4)) +
  labs(x='Axis 1 (44.8%)', y= 'Axis 2 (11.6%)') +
  labs(color= "LPS") +
  labs(shape= "Treatment") +
  guides(size=FALSE) +
  theme(strip.text = element_text(size = 14, face = "bold")) +
  theme(legend.text = element_text(size=14)) +
  theme(legend.title = element_text(size = 14, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold"), axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 14), axis.text.y = element_text(color = "black", size = 14)) +
  theme(text=element_text(family="Times New Roman"))

LPSotu <- OTU.clean[rownames(OTU.clean) %in% rownames(LPS),]
LPS
LPSp <- adonis(LPSotu ~ LPS$treatment, permutations = 999)
LPSp
LPSe <- adonis(LPSotu ~ LPS$LPS, permutations = 999)
LPSe

#####-----STATISTICAL ANALYSIS---------
## PERMANOVA using adonis function

#subet beta diversity for day 13 
calffecal %>% tally()
calffecal %>% count(Time)
calffecal %>% count(treatment)


## We need the metadata and the OTU table
str(calffecal)
str(OTU.clean)

#PERMANOVA 
Time <- adonis(OTU.clean ~ calffecal$d, strata=calffecal$calf, permutations = 999)
Time

#we can also performe a pairwise comparison with the function 
# Pairwise Adonis funtion by edro Martinez Arbizu & Sylvain Monteux
#https://github.com/pmartinezarbizu/pairwiseAdonis/blob/master/pairwiseAdonis/R/pairwise.adonis.R

pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {
  
  #describe parent call function 
  ststri <- ifelse(is.null(strata),'Null',strata)
  fostri <- as.character(x)
  #list to store results
  
  #copy model formula
  x1 <- x
  # extract left hand side of formula
  lhs <- x1[[2]]
  # extract factors on right hand side of formula 
  rhs <- x1[[3]]
  # create model.frame matrix  
  x1[[2]] <- NULL   
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE) 
  
  # create unique pairwise combination of factors 
  co <- combn(unique(as.character(rhs.frame[,1])),2)
  
  # create names vector   
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
  }
  #create results list  
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres
  
  #add parent call to res 
  res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri, ', permutations',nperm ))
  
  
  #start iteration trough pairwise combination of factors  
  for(elem in 1:ncol(co)){
    
    #reduce model elements  
    if(inherits(eval(lhs),'dist')){	
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                           rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
    }else{
      xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    }
    
    mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),] 
    
    # redefine formula
    if(length(rhs) == 1){
      xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))	
    }else{
      xnew <- as.formula(paste('xred' , 
                               paste(rhs[-1],collapse= as.character(rhs[1])),
                               sep='~'))}
    
    #pass new formula to adonis
    if(is.null(strata)){
      ad <- adonis2(xnew,data=mdat1, ... )
    }else{
      perm <- how(nperm = nperm)
      setBlocks(perm) <- with(mdat1, mdat1[,ststri])
      ad <- adonis2(xnew,data=mdat1,permutations = perm, ... )}
    
    res[nameres[elem+1]] <- list(ad[1:5])
  }
  #names(res) <- names  
  class(res) <- c("pwadstrata", "list")
  return(res)
} 


### Method summary


# Pairwise comparison
TimePair <- pairwise.adonis2(OTU.clean ~ d,  data=calffecal, strata="calf")
TimePair

###------------------------------------ Test the result--- Separation of each data into day
# We want to measure the effect of each treatment, bacth and block in each timepoint 0d, 28d, 56d, 84d, and 112
#Permanova d0

str(day0)
day0$treatment <- factor(day0$treatment)

##Preparing OTU table
# this makes sure there are no samples in the OTU table that are not in our metadata
str(OTU.clean)
str(day0)
otu_sub0 <- OTU.clean[rownames(OTU.clean) %in% rownames(day0),]

Perd0 <- adonis(otu_sub0 ~ day0$treatment + day0$batch + day0$block + day0$treatment:day0$batch + day0$treatment:day0$block, permutations = 999)
Perd0 #Batch sign P =0.003

#Permanova 28 d
day28 <- subset(calffecal,
               subset=d %in% c("28"))
str(day28)
day28$treatment <- factor(day28$treatment)

##Preparing OTU table
# this makes sure there are no samples in the OTU table that are not in our metadata
str(OTU.clean)
str(day28)
otu_sub28 <- OTU.clean[rownames(OTU.clean) %in% rownames(day28),]

Perd28 <- adonis(otu_sub28 ~ day28$treatment + day28$batch + day28$block + day28$treatment:day28$batch + day28$treatment:day28$block, permutations = 999)
Perd28 #batch P= 0.001

#Permanova 56 d
day56 <- subset(calffecal,
                subset=d %in% c("56"))
str(day56)
day56$treatment <- factor(day56$treatment)

##Preparing OTU table
# this makes sure there are no samples in the OTU table that are not in our metadata
str(OTU.clean)
str(day56)
otu_sub56 <- OTU.clean[rownames(OTU.clean) %in% rownames(day56),]

Perd56 <- adonis(otu_sub56 ~ day56$treatment + day56$batch + day56$block + day56$treatment:day56$batch + day56$treatment:day56$block, permutations = 999)
Perd56 #batch is sig P=0.001

#Permanova 84 d
day84 <- subset(calffecal,
                subset=d %in% c("84"))
str(day84)
day84$treatment <- factor(day84$treatment)

##Preparing OTU table
# this makes sure there are no samples in the OTU table that are not in our metadata
str(OTU.clean)
str(day84)
otu_sub84 <- OTU.clean[rownames(OTU.clean) %in% rownames(day84),]

Perd84 <- adonis(otu_sub84 ~ day84$treatment + day84$batch + day84$block + day84$treatment:day84$batch + day84$treatment:day84$block, permutations = 999)
Perd84 #batch P= 0.001 and block P = 0.018

BPair84 <- pairwise.adonis2(otu_sub84 ~ block,  data=day84)
BPair84 #no difference in any comparison

#Permanova 112 d
day112 <- subset(calffecal,
                subset=d %in% c("112"))
str(day112)
day112$treatment <- factor(day112$treatment)

##Preparing OTU table
# this makes sure there are no samples in the OTU table that are not in our metadata
str(OTU.clean)
str(day112)
otu_sub112 <- OTU.clean[rownames(OTU.clean) %in% rownames(day112),]

Perd112 <- adonis(otu_sub112 ~ day112$treatment + day112$batch + day112$block + day112$treatment:day112$batch + day112$treatment:day112$block, permutations = 999)
Perd112 # treatment, P 0.043 batch P= 0.001 and block P=0.017

BPair112 <- pairwise.adonis2(otu_sub112 ~ block,  data=day112)
BPair112 #no difference

####---------------------------- DISPERSION TEST -----------------------------
### Using vegan's betadisper() function ###
DistWU = read.csv("weighted-distance-matrix.tsv", header=TRUE, sep="\t", stringsAsFactors = F)
row.names(DistWU) <- DistWU[,1]
DistWU = DistWU[,-c(1)]
distWU <- as.dist(DistWU)
#check homogeneity of variance

group <- (calffecal$d*calffecal$calf)

timeD <- betadisper(distWU, type = c("centroid"), group = calffecal$d)
timeD
boxplot(timeD)

anova(timeD)
ptimeD<- permutest(timeD, permutations = 999)
ptimeD# significant 
TukeyHSD(timeD)

#------Day 28
Dist28 = read.csv("dist28.txt", header=TRUE, sep="\t", stringsAsFactors = F)
row.names(Dist28) <- Dist28[,1]
Dist28 = Dist28[,-c(1)]
dist28 <- as.dist(Dist28)

#treatment
D28 <- betadisper(dist28, type = c("centroid"), group = day28$treatment)
D28
boxplot(D28)


pD28<- permutest(D28, permutations = 999)#, pairwise = TRUE)
pD28

# Batch
D28 <- betadisper(dist28, type = c("centroid"), group = day28$batch)
D28
boxplot(D28)

pD28<- permutest(D28, permutations = 999)#, pairwise = TRUE)
pD28

# Batch
D28 <- betadisper(dist28, type = c("centroid"), group = day28$block)
D28
boxplot(D28)

pD28<- permutest(D28, permutations = 999)#, pairwise = TRUE)
pD28

#------D 56
Dist56 = read.csv("dist56.txt", header=TRUE, sep="\t", stringsAsFactors = F)
row.names(Dist56) <- Dist56[,1]
Dist56 = Dist56[,-c(1)]
dist56 <- as.dist(Dist56)

#Treatment
D56 <- betadisper(dist56, type = c("centroid"), group = day56$treatment)
D56
boxplot(D56)

pD56<- permutest(D56, permutations = 999)#, pairwise = TRUE)
pD56

#Batch
D56 <- betadisper(dist56, type = c("centroid"), group = day56$batch)
D56
boxplot(D56)

pD56<- permutest(D56, permutations = 999)#, pairwise = TRUE)
pD56 # <-- batch is significant

#Block 
D56 <- betadisper(dist56, type = c("centroid"), group = day56$block)
D56
boxplot(D56)

pD56<- permutest(D56, permutations = 999)#, pairwise = TRUE)
pD56 

#------D 84
Dist84 = read.csv("dist84.txt", header=TRUE, sep="\t", stringsAsFactors = F)
row.names(Dist84) <- Dist84[,1]
Dist84 = Dist84[,-c(1)]
dist84 <- as.dist(Dist84)

#Treatment
D84 <- betadisper(dist84, type = c("centroid"), group = day84$treatment)
D84
boxplot(D84)


pD84<- permutest(D84, permutations = 999)#, pairwise = TRUE)
pD84

#Batch
D84 <- betadisper(dist84, type = c("centroid"), group = day84$batch)
D84
boxplot(D84)


pD84<- permutest(D84, permutations = 999)#, pairwise = TRUE)
pD84 # 

#Block 
D84 <- betadisper(dist84, type = c("centroid"), group = day84$block)
D84
boxplot(D84)

pD84<- permutest(D84, permutations = 999)#, pairwise = TRUE)
pD84 

#------D 112
Dist112 = read.csv("dist112.txt", header=TRUE, sep="\t", stringsAsFactors = F)
row.names(Dist112) <- Dist112[,1]
Dist112 = Dist112[,-c(1)]
dist112 <- as.dist(Dist112)

#Treatment
D112 <- betadisper(dist112, type = c("centroid"), group = day112$treatment)
D112
boxplot(D112)

pD112<- permutest(D112, permutations = 999)#, pairwise = TRUE)
pD112

#Batch
D112 <- betadisper(dist112, type = c("centroid"), group = day112$batch)
D112
boxplot(D112)

pD112<- permutest(D112, permutations = 999)#, pairwise = TRUE)
pD112 #-- Significant

#Block 
D112 <- betadisper(dist112, type = c("centroid"), group = day112$block)
D112
boxplot(D112)


pD112<- permutest(D112, permutations = 999)#, pairwise = TRUE)
pD112 
#------D 0
Dist0 = read.csv("distd0.txt", header=TRUE, sep="\t", stringsAsFactors = F)
row.names(Dist0) <- Dist0[,1]
Dist0 = Dist0[,-c(1)]
dist0 <- as.dist(Dist0)

#Treatment
D0 <- betadisper(dist0, type = c("centroid"), group = day0$treatment)
D0
boxplot(D0)

disperD0D <- data.frame(Distance_to_centroid=D0$distances,Group=D0$group)

pD0<- permutest(D0, permutations = 999)#, pairwise = TRUE)
pD0

#Batch
D0 <- betadisper(dist0, type = c("centroid"), group = day0$batch)
D0
boxplot(D0)

disperD0D <- data.frame(Distance_to_centroid=D0$distances,Group=D0$group)

pD0<- permutest(D0, permutations = 999)#, pairwise = TRUE)
pD0 #-- Significant

#Block 
D0 <- betadisper(dist0, type = c("centroid"), group = day0$block)
D0
boxplot(D0)

disperD0D <- data.frame(Distance_to_centroid=D0$distances,Group=D0$group)

pD0<- permutest(D0, permutations = 999)#, pairwise = TRUE)
pD0 

