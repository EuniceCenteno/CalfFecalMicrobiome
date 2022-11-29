## Regression Analysis
## Random Forest prediction using regression
# the prediction is based on the average prediction accross all tree ~ regression

library(ggplot2)
library(vegan)
library(dplyr)
library(magrittr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(randomForest)
library(knitr)
library(qiime2R)
library(tidyr) #for separate function
library(naniar)# for replace all function
library(ggpubr)
library(pheatmap)
library(zoo)

setwd("~/Desktop/eunice/MSc/CalfFecal/CalfFecal01_2021/NewQiime/TrimmedFecal/exported/")
#We mainly have two varieties, planted in two locations. Here we will first demonstrate with A50 modeling and IR24 verification scheme. 
#This experiment is more complicated , and there are many combinations of specific methods.

## we need to remove day 0 because these samples we collected before receiving the treatment
metadata = read.table("CalfFecal_metaNod0.txt",header = T, row.names = 1)
metadata$calf = as.factor(metadata$calf) #214
str(metadata)
metadata$d <- as.numeric(metadata$d)
metadata$treatment <- as.factor(metadata$treatment)
levels(metadata$treatment) <- list("CON"="Pre-CON", "SCFP"="Pre-SCFP", "CON"="CON","SCFP"="SCFP")
metadata$id <- rownames(metadata)
str(metadata)

#otu_table
ASVs <- read_qza("table.qza")
ASV_s <- as.data.frame(ASVs$data)
ASV_table <- as.data.frame(ASVs$data) #7529 ASVs
ASV_table$ASVnos <- paste0("ASV", 1:nrow(ASV_table))
ASV_table$ASVstring <- rownames(ASV_table)
rownames(ASV_table) <- ASV_table$ASVnos ##We change the ASV name created in Qiime to ASVn
ASVkey <- ASV_table[, (ncol(ASV_table)-1):ncol(ASV_table)] #the key withe the names
ASV_table <- ASV_table[,-(ncol(ASV_table)-1):-ncol(ASV_table)]
ASV_table <- t(ASV_table) #calf id as row names

#change the rownames in the ASV table
ASV_table <- merge(metadata, ASV_table, by.x = "sample", 0)

rownames(ASV_table) <- ASV_table$id
ASV_table <- ASV_table[,-c(1:11)]
ASV_table <- as.matrix(ASV_table) #7529 we change the sample name 

#adding taxonomy
#Taxonomy of each OTU
tax <- read_qza("taxonomy.qza")
tax <- as.data.frame(tax$data)
tax2 = separate(tax, Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
#This warning means that some cells are empty and that R is replacing empty cells with NA. Now there are other cells that are unclassified that  say, for example `s__` or `g__`. All these need to be removed and replaced with `NA`. 
#All this is OK except that in future use of the taxonomy table, these ASVs will be ignored because they are not classified. Why are ASVs not classified? Its because there is not a close enough match in the database. Just because there is not a good match in the database does not mean they don‚Äôt exist, so I wanted to make sure this data was not lost. So in my new code, from lines 300 ‚Äì 334 I make it so that ASVs that are unclassified at any level are classified as the lowest taxonomic level for which there is a classification.

#All the strings that need to be removed and replaced with NA
na_strings <- c(" s__", " g__", " f__", " o__", " c__")

tax3 = replace_with_na_all(tax2, condition = ~.x %in% na_strings)

#This code is great because an ASV that is unclassified at a certain level are all listed as `NA`.
#Unfortunately this command changed ou Feature.ID names

#Next, all these `NA` classifications with the last level that was classified
tax3[] <- t(apply(tax3, 1, zoo::na.locf))
tax3 <- as.data.frame(tax3)
row.names(tax3) <- tax3[,1]
tax3 = tax3[,-c(1:2)]
tax.clean <- as.data.frame(tax3)
tax.clean$OTUs <- rownames(tax.clean)
#Would be good to check here to make sure the order of the two data frames was the same. You should do this on your own.

###Remove all the OTUs that don't occur in our OTU.clean data set
tax.final = tax.clean[row.names(tax.clean) %in% row.names(ASV_s),]

##Remove unneccesary information from the taxonomy names
tax.final$Phylum <- sub(" p__*", "", tax.final[,1])
tax.final$Phylum <- sub(" k__*", "", tax.final[,1])
tax.final$Class <- sub(" p__*", "", tax.final[,2])
tax.final$Class <- sub(" c__*", "", tax.final[,2])
tax.final$Class <- sub(" o__*", "", tax.final[,2])
tax.final$Class <- sub(" k__*", "", tax.final[,2])
tax.final$Order <- sub(" p__*", "", tax.final[,3])
tax.final$Order <- sub(" c__*", "", tax.final[,3])
tax.final$Order <- sub(" o__*", "", tax.final[,3])
tax.final$Order <- sub(" o__*", "", tax.final[,3])
tax.final$Order <- sub(" k__*", "", tax.final[,3])
tax.final$Family <- sub(" p__*", "", tax.final[,4])
tax.final$Family <- sub(" c__*", "", tax.final[,4])
tax.final$Family <- sub(" o__*", "", tax.final[,4])
tax.final$Family <- sub(" f__*", "", tax.final[,4])
tax.final$Family <- sub(" k__*", "", tax.final[,4])
tax.final$Genus <- sub(" p__*", "", tax.final[,5])
tax.final$Genus <- sub(" c__*", "", tax.final[,5])
tax.final$Genus <- sub(" o__*", "", tax.final[,5])
tax.final$Genus <- sub(" f__*", "", tax.final[,5])
tax.final$Genus <- sub(" g__*", "", tax.final[,5])
tax.final$Genus <- sub(" k__*", "", tax.final[,5])
tax.final$Species <- sub(" p__*", "", tax.final[,6])
tax.final$Species <- sub(" c__*", "", tax.final[,6])
tax.final$Species <- sub(" o__*", "", tax.final[,6])
tax.final$Species <- sub(" f__*", "", tax.final[,6])
tax.final$Species <- sub(" g__*", "", tax.final[,6])
tax.final$Species <- sub(" s__*", "", tax.final[,6])
tax.final$Species <- sub(" k__*", "", tax.final[,6])

TaxASV <- merge(tax.final, ASVkey, by.x = 0, by.y = "ASVstring")
# some ASVs where not correctly taxonomically assigned
row.names(TaxASV) <- TaxASV[,10]
TaxW <- c("ASV6720", "ASV7520", "ASV7399", "ASV5508", "ASV1017","ASV5132", "ASV6494", "ASV6982",
          "ASV4282", "ASV6599", "ASV6386", "ASV4611", "ASV2772", "ASV4592", "ASV6115", "ASV6146", 
          "ASV3450", "ASV5683", "ASV6609")
TaxASV<- TaxASV[!(row.names(TaxASV) %in% c(TaxW)), ]
TaxASV = TaxASV[,-c(1,10)]

#removing 
ASV_table <- t(ASV_table)
ASV_table<- ASV_table[!(row.names(ASV_table) %in% c(TaxW)), ]
ASV_table <- t(ASV_table) # 7510 ASVs and 214 samples

##subsetting for treatment 
#separating the treatment 
CON <- subset(metadata, treatment == "CON")
CON <-CON[-c(5),]
SCFP <- subset(metadata, treatment == "SCFP")

#New to generate the OTU for the healthy period
C_OTU <- ASV_table[rownames(ASV_table) %in% rownames(CON),] #105
S_OTU <- ASV_table[rownames(ASV_table) %in% rownames(SCFP),] #108

#---------------MIKING THE MODEL CON 
#generating a dataframe with all the response (days) from the samples
# Make one column for our outcome/response variable 
OTU.physeq = otu_table(as.matrix(ASV_table), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(TaxASV))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(metadata)

#We then merge these into an object of class phyloseq.
physeq_deseq = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_deseq # [ 7510 taxa and 214 samples ]
colnames(tax_table(physeq_deseq))
## Filter any non-baxteria, chloroplast and mitochondria
physeq_deseq %>%
  subset_taxa(Family != "Mitochondria" & 
                Genus != "Mitochondria" &
                Species != "Mitochondria" &
                Order != "Chloroplast" &
                Family != "Chloroplast" &
                Genus != "Chloroplast" &
                Species != "Chloroplast") -> physeq_deseq
physeq_deseq ##[ 7510 taxa and 214 samples ]
## Random forest, we want to make the model so it classify samples based on the age
ntaxa(physeq_deseq) #total of 7510 ASVs

#making the data frame with all the combine responses
response <- sample_data(physeq_deseq)$d ##classification
class(response)
predictors <- (otu_table(physeq_deseq))
dim(predictors) # we have 214 7510

# Combine them into 1 data frame
rf.data <- data.frame(response, predictors) ## complete dataframe with the sample ID, classification and OTU table
str(rf.data)
tail(rf.data$response)
class(rf.data)

##Running the model first for the CON group that will be used to predict age in the SCFP animals 
CON <- rf.data[rownames(rf.data) %in% rownames(CON),]
class (CON)
CON$Rownames <- rownames(CON)
tail(CON)
CON <- CON[,-c(2:7511)]

##First we need to specify the taxa with the highes MSE that will be used to model the equation]
#Training set
set.seed(315)
rf = randomForest(C_OTU, CON$response, importance=TRUE, proximity=TRUE)
print(rf)
plot(rf)
which.min(rf$mse) ## we need 238 trees to get the lowest MSE

#tying 238 tress
rf1 = randomForest(C_OTU, CON$response, importance=TRUE, proximity=TRUE, ntree=238)
print(rf1)
plot(rf1)
which.min(rf1$mse)

#rfcv：Random Forest Cross Validation
result = rfcv(C_OTU, CON$response, cv.fold=10)
result$error.cv #corresponding vector of error rates or MSEs at each step
#the mean squared error (MSE) or mean squared deviation (MSD) of an estimator 
#(of a procedure for estimating an unobserved quantity) measures the average of the squares of the errors—that is, the average squared difference between the estimated values and the actual value.
#plot
#10% of the data is used to test the function
# Draw validation, results
with(result, plot(n.var, error.cv, log="x", type="o", lwd=2))

# Export feature importance
imp = as.data.frame(rf1$importance)
imp = imp[order(imp[,1],decreasing = T),]
head(imp)

# Simple visualization
varImpPlot(rf, main = "Top 40-Feature OTU importance",n.var = 40, bg = par("bg"),
           color = par("fg"), gcolor = par("fg"), lcolor = "gray")

# Analysis of selected top15 group best
imp = head(imp, n=40) #I'm using the first 62 because is the smallest value with the lowest MSE
# Reverse the order X -axis, so the histogram from the draw down
imp
#colnames(imp) <- c("IncMSE",  "IncNodePurity")

# imp species name decomposition
# Add Class level to remain in line
str(TaxASV)
imp <- merge(imp, TaxASV, by.x = 0, by.y = 0)
imp
imp = imp[order(imp[,2],decreasing = T),]
#write.table(imp,file = "importance_classRF1-4.txt",quote = F,sep ='\t', row.names = T, col.names = T)

##PLot the data
imp =read.csv("importance_classRF1-4.csv", na.strings = c("","NA"), header=TRUE)
imp = head(imp, n=40)
imp = imp[order(imp[,2],decreasing = T),]

#count how many ASVs are for each genera
imp$Genus <- as.factor(imp$Genus)
imp %>% tally()
imp %>% count(Genus)

#Draw a time series heat map
# Data filtering 15 Ge feature show
rownames(imp) <- imp$Row.names
str(C_OTU)
TASV <- t(C_OTU)
sampFile <- merge(imp[,-c(1:12)], TASV, by.x =0, by.y = 0)
rownames(sampFile) <- sampFile$Row.names
sampFile <- as.data.frame(t(sampFile[,-c(1)]))
sampMeta <- merge(sampFile[,-c(1:40)], CON, by.x = 0, by.y = 0)
rownames(sampMeta) <- sampMeta$Row.names
sampMeta <- sampMeta[,-c(1)]

sampFile2 = as.data.frame(sampMeta$response,row.names = row.names(sampMeta))
colnames(sampFile2)[1] = "group"
mat_t = sampFile ## OTU Table sample in the rows and bacteria in the columns
mat_t2 = merge(sampFile2, mat_t, by="row.names") #combination between otu table, ID and day
mat_t2 = mat_t2[,-1] #remove the ID
mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
otu_norm_group2 = do.call(rbind, mat_mean)[-1,] #add the mean
colnames(otu_norm_group2) = mat_mean$group #clumn and row names
colnames(otu_norm_group2) <- c("1", "2", "3", "4")
otu_norm_group2 <- as.data.frame(otu_norm_group2)
otu_norm_group2 <- merge(otu_norm_group2, imp[-c(1:9,11:12)], by.x = 0, by.y = 0)
rownames(otu_norm_group2) <- otu_norm_group2$SpeciesASV
otu_norm_group2 <- otu_norm_group2[,-c(1,6)]
otu_norm_group2 <- as.matrix(otu_norm_group2)
pheatmap(otu_norm_group2,scale="row",cluster_cols = F, cluster_rows = F)
#pheatmap(otu_norm_group, scale="column",cluster_cols = T, cluster_rows = F, filename = "heatmap_groups.pdf", width = 5, height = 5)

## Find the maximum value of each group
str(imp)
bak=otu_norm_group2
otu_norm_group2 = otu_norm_group2[as.character(imp$SpeciesASV),] #By initial sort
for (i in 1:length(rownames(otu_norm_group2))) {
  #  i=1
  x=as.data.frame(sort(otu_norm_group2[i,],decreasing = T))
  imp[i,"order"]=rownames(x)[1]
}

taxonomy2 = arrange(imp, desc(order), SpeciesASV)
otu_norm_group3 = otu_norm_group2[match(taxonomy2$SpeciesASV,rownames(otu_norm_group2)),] #By initial sort

newnames <- lapply(
  rownames(otu_norm_group3),
  function(x) bquote(italic(.(x))))

p2=pheatmap(otu_norm_group3,scale="row",cluster_cols = F, cluster_rows = F, angle_col = "0",  labels_row = as.expression(newnames))
p2# Sort by initial, the 40 ASV in the CON samples

#### SCFP animals
str(SCFP)

#Running the RandomForest regression, day should be as numeric
SCFP <- rf.data[rownames(rf.data) %in% rownames(SCFP),]
class(SCFP)
SCFP$Rownames <- rownames(SCFP)
tail(SCFP)
SCFP <- SCFP[,-c(2:7511)]

set.seed(315)
rf2 = randomForest(S_OTU, SCFP$response, importance=TRUE, proximity=TRUE)
print(rf2)
plot(rf2)
which.min(rf2$mse) 
rf3 = randomForest(S_OTU, SCFP$response, importance=TRUE, proximity=TRUE, ntree=238) #use the same number of trees as CON
print(rf3)

## Cross-Validation Selection Features
set.seed(315) #Random data guarantees that the results can be repeated, must
# rfcv is the random forest cross validation function: Random Forest Cross Validation
result2 = rfcv(S_OTU, SCFP$response, cv.fold=10)
#10% of the data is used to test the function
# View the error rate table, 23 the lowest error rate, is the best model
result2$error.cv
# Draw validation, results
with(result2, plot(n.var, error.cv, log="x", type="o", lwd=2))


# Export feature importance
imp2 = as.data.frame(rf3$importance)
imp2 = imp2[order(imp2[,1],decreasing = T),]
head(imp2)

# Simple visualization
varImpPlot(rf2, main = "Top 40-Feature OTU importance",n.var = 40, bg = par("bg"),
           color = par("fg"), gcolor = par("fg"), lcolor = "gray")

# Analysis of selected top15 group best
imp2 = head(imp2, n=40) #I'm using the first 31
# Reverse the order X -axis, so the histogram from the draw down
imp2
#colnames(imp2) <- c("IncMSE",  "IncNodePurity")

# imp species name decomposition
# Add Class level to remain in line
str(TaxASV)
imp2 <- merge(imp2, TaxASV, by.x = 0, by.y = 0)
imp2
imp2 = imp2[order(imp2[,2],decreasing = T),]
#write.table(imp2,file = "importance_classSCFP1-4.txt",quote = F,sep ='\t', row.names = T, col.names = T)

##PLot the data
imp2 = read.csv("importance_classSCFP1-4.csv", na.strings = c("","NA"), header=TRUE)
imp2 = head(imp2, n=40)
imp2 = imp2[order(imp2[,2],decreasing = T),]

#count how many ASVs are for each genera
imp2$Genus <- as.factor(imp2$Genus)
imp2 %>% tally()
imp2 %>% count(Genus)

#. Draw a time series heat map
# Data filtering 15 Ge feature show
rownames(imp2) <- imp2$Row.names
str(imp2)
str(S_OTU)
TASVS <- t(S_OTU)
sampFile3 <- merge(imp2[,-c(1:12)], TASVS, by.x = 0, by.y = 0)
rownames(sampFile3) <- sampFile3$Row.names
sampFile3 <- as.data.frame(t(sampFile3[,-c(1)]))
sampMeta2 <- merge(sampFile3[,-c(1:40)], SCFP, by.x = 0, by.y = 0)
rownames(sampMeta2) <- sampMeta2$Row.names
sampMeta2 <- sampMeta2[,-c(1)]

sampFile4 = as.data.frame(sampMeta2$response,row.names = row.names(sampMeta2))
colnames(sampFile4)[1] = "group"
mat_t = sampFile3 ## OTU Table sample in the rows and bacteria in the columns
mat_t2 = merge(sampFile4, mat_t, by="row.names") #combination between otu table, ID and day
mat_t2 = mat_t2[,-1] #remove the ID
mat_mean2 = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
otu_norm_group4 = do.call(rbind, mat_mean2)[-1,] #add the mean
colnames(otu_norm_group4) = mat_mean$group #clumn and row names
colnames(otu_norm_group4) <- c("1", "2", "3", "4")
otu_norm_group4 <- as.data.frame(otu_norm_group4)
otu_norm_group4 <- merge(otu_norm_group4, imp2[-c(1:9,11:12)], by.x = 0, by.y = 0)
rownames(otu_norm_group4) <- otu_norm_group4$SpeciesASV
otu_norm_group4 <- otu_norm_group4[,-c(1,6)]
otu_norm_group4 <- as.matrix(otu_norm_group4)
pheatmap(otu_norm_group4,scale="row",cluster_cols = F, cluster_rows = F)
#pheatmap(otu_norm_group, scale="column",cluster_cols = T, cluster_rows = F, filename = "heatmap_groups.pdf", width = 5, height = 5)

## Find the maximum value of each group
str(imp2)
bak=otu_norm_group4
otu_norm_group4 = otu_norm_group4[as.character(imp2$SpeciesASV),] #By initial sort
for (i in 1:length(rownames(otu_norm_group4))) {
  #  i=1
  x=as.data.frame(sort(otu_norm_group4[i,],decreasing = T))
  imp2[i,"order"]=rownames(x)[1]
}

taxonomy3 = arrange(imp2, desc(order), SpeciesASV)
otu_norm_group5 = otu_norm_group4[match(taxonomy3$SpeciesASV,rownames(otu_norm_group4)),] #By initial sort
newnames2 <- lapply(
  rownames(otu_norm_group5),
  function(x) bquote(italic(.(x))))

p3=pheatmap(otu_norm_group5,scale="row",cluster_cols = F, cluster_rows = F, legend_labels = c("0",  "28", "56", "84", "112"), angle_col = "0", labels_row = as.expression(newnames2))
p3# Sort by initial, the 40 ASV in the CON samples


#The scale in the heapmap is the Z-score
#Check which ASVs are in the two groups, 
CompC = merge(otu_norm_group3, taxonomy2[,-c(2:9, 11:13)], by.x= 0, by.y="SpeciesASV" ) #40 total
CompC$treatment <- c("CON")
CompS = merge(otu_norm_group5, taxonomy3[,-c(2:9, 11:13)], by.x= 0, by.y="SpeciesASV" ) #40 total
CompS$treatment <- c("SCFP")
Comp = merge(CompC, CompS, by.x = "Row.names.y", by.y = "Row.names.y")
rm(Comp1)
Comp1 <- rbind(CompC, CompS)
str(Comp1)
#write.table(Comp,file = "Comp2.txt",quote = F,sep ='\t', row.names = T, col.names = T)

#Look at the difference between the two groups
##venn diagram
#install.packages("VennDiagram")
library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(2, "Pastel2")

venn.diagram(
  x = list(
    Comp1 %>% filter(treatment=="SCFP") %>% select(Row.names.y) %>% unlist() , 
    Comp1 %>% filter(treatment=="CON") %>% select(Row.names.y) %>% unlist()
  ),
  category.names = c("SCFP" , "CON"),
  filename = '#14_venn_diagramm.png',
  output=TRUE,
  col=c('turquoise2', 'indianred1'),
  fill = c(alpha('turquoise2',0.3), alpha('indianred1',0.3))
  )

list(Comp$Row.names.x)

dev.off()
CompCon <- Comp[,c(1:7)]
CompCon = arrange(CompCon, Row.names.x) #alphabetic order
rownames(CompCon) <- CompCon$Row.names.x
CompCon1 <- CompCon[,-c(1,2,7)] # we use it to make the heatmap
colnames(CompCon1) <- c("1", "2", "3", "4")
newnames3 <- lapply(
  rownames(CompCon1),
  function(x) bquote(italic(.(x))))
p4=pheatmap(CompCon1 ,scale="row",cluster_cols = F, cluster_rows = F,  angle_col = "0", labels_row = as.expression(newnames3))
p4# Sort by initial

rm(CompSC1)
CompSC <- Comp[,c(8:13)]
CompSC = arrange(CompSC, Row.names.y.y) #alphabetic order
rownames(CompSC) <- CompSC$Row.names.y.y
CompSC1 <- CompSC[,-c(1,6)]
colnames(CompSC1) <- c("1", "2", "3", "4")
newnames4 <- lapply(
  rownames(CompSC1),
  function(x) bquote(italic(.(x))))
p5=pheatmap(CompSC1 ,scale="row",cluster_cols = F, cluster_rows = F,  angle_col = "0", labels_row = as.expression(newnames4))
p5# Sort by initial

## Things that are not common
dev.off()
CompCon$Row.names.x
Cnames <- c(CompCon$Row.names.x)  ## the 18 shared ones
otu_norm_1 <- otu_norm_group2[!(row.names(otu_norm_group2) %in% c(Cnames)), ]
p6=pheatmap(otu_norm_1 ,scale="row",cluster_cols = F, cluster_rows = F, , angle_col = "0")
p6# Sort by initial
rownames(CompCon1)
CompCon1$SpeciesASV <- rownames(CompCon1)
imp3 <- anti_join(imp, CompCon1, by=c("SpeciesASV"))  #21 not shared
for (i in 1:length(rownames(otu_norm_1))) {
  #  i=1
  x=as.data.frame(sort(otu_norm_1[i,],decreasing = T))
  imp3[i,"order"]=rownames(x)[1]
}
otu_norm_1 <- as.matrix(otu_norm_1 )
taxonomy4 = arrange(imp3, desc(order), SpeciesASV)
otu_norm_1 = otu_norm_1[match(taxonomy4$SpeciesASV,rownames(otu_norm_1)),] #By initial sort
newnames5 <- lapply(
  rownames(otu_norm_1),
  function(x) bquote(italic(.(x))))
p6=pheatmap(otu_norm_1,scale="row",cluster_cols = F, cluster_rows = F, angle_col = "0", labels_row = as.expression(newnames5))
p6# Sort by initial, the 21 ASV in the CON samples not shared
write.csv(otu_norm_1, "CONotu_norm_1.csv")

### SCFP not shared
dev.off()
CompSC$Row.names.y.y
Snames <- rownames(CompSC) ## the 18 shared ones
otu_norm_2 <- otu_norm_group4[!(row.names(otu_norm_group4) %in% c(Snames)), ]
p7=pheatmap(otu_norm_2 ,scale="row",cluster_cols = F, cluster_rows = F,  angle_col = "0")
p7# Sort by initial
CompSC1$SpeciesASV <- rownames(CompSC1)
imp4 <- anti_join(imp2, CompSC1, by=c("SpeciesASV")) #22 not shared
for (i in 1:length(rownames(otu_norm_2))) {
  #  i=1
  x=as.data.frame(sort(otu_norm_2[i,],decreasing = T))
  imp4[i,"order"]=rownames(x)[1]
}
taxonomy5 = arrange(imp4, desc(order), SpeciesASV)
otu_norm_2 = otu_norm_2[match(taxonomy5$SpeciesASV,rownames(otu_norm_2)),] #By initial sort
newnames6 <- lapply(
  rownames(otu_norm_2),
  function(x) bquote(italic(.(x))))
p7=pheatmap(otu_norm_2,scale="row",cluster_cols = F, cluster_rows = F, angle_col = "0", labels_row = as.expression(newnames6))
p7# Sort by initial, the 21 ASV not shared in the SCP samples
write.csv(otu_norm_2, "SCFPotu_norm_2.csv")

#not common
NotSharesSCFP <- as.data.frame(otu_norm_2)
NotSharedCON <- as.data.frame(otu_norm_1)
rownames(NotSharesSCFP)
rownames(NotSharedCON)
Comp$Row.names.x

##make plots of the RF models and the MSE
library(gridExtra)
library(RColorBrewer)

makeColors <- function(){
  maxColors <- 10
  usedColors <- c()
  possibleColors <- colorRampPalette( brewer.pal( 9 , "Set1" ) )(maxColors)
  
  function(values){
    newKeys <- setdiff(values, names(usedColors))
    newColors <- possibleColors[1:length(newKeys)]
    usedColors.new <-  c(usedColors, newColors)
    names(usedColors.new) <- c(names(usedColors), newKeys)
    usedColors <<- usedColors.new
    
    possibleColors <<- possibleColors[length(newKeys)+1:maxColors]
    usedColors
  }
} 

mkColor <- makeColors()


set.seed(1)
df1 <- data.frame(c=imp$Phylum, x=imp$SpeciesASV,  y=imp$X.IncMSE)
colnames(df1) <- c("Phylum", "SpeciesASV", "IncMSE")
df2 <- data.frame(c=imp2$Phylum, x=imp2$SpeciesASV,  y=imp2$X.IncMSE)
colnames(df2) <- c("Phylum", "SpeciesASV", "IncMSE")

p=ggplot(data = df1, mapping = aes(x=reorder(SpeciesASV, IncMSE), y= IncMSE,fill=Phylum)) +
  geom_bar(stat="identity")+coord_flip() +theme_bw() +xlab("ASV ID") +
  scale_fill_manual(values = mkColor(df1$Phylum)) +
  theme(legend.text = element_text(size=12, face = "italic")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12, face="italic")) 
p
p1=ggplot(data = df2, mapping = aes(x=reorder(SpeciesASV, IncMSE), y= IncMSE,fill=Phylum)) +
  geom_bar(stat="identity")+coord_flip() +theme_bw() +xlab("ASV ID") +
  scale_fill_manual(values = mkColor(df2$Phylum))+
  theme(legend.text = element_text(size=12, face = "italic")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12, face="italic")) 
p1
grid.arrange(p, p1, ncol=2)
