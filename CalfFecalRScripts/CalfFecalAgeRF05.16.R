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
library(afex)
library(emmeans)
library(lme4)

setwd("~/Desktop/eunice/MSc/CalfFecal/CalfFecal01_2021/NewQiime/TrimmedFecal/exported/")
getwd()
rm(list = ls()) # clears R memory

##Data import
#metadata
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
ASV_table <- as.matrix(ASV_table) #7529 we change the sample name and 214 samples

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
with(result, plot(n.var, error.cv, log="x", type="o", lwd=2))

# Export feature importance
imp = as.data.frame(rf1$importance)
imp = imp[order(imp[,1],decreasing = T),]
head(imp)

# Simple visualization
varImpPlot(rf, main = "Top 49-Feature OTU importance",n.var = 40, bg = par("bg"),
           color = par("fg"), gcolor = par("fg"), lcolor = "gray")

# Analysis of selected top49 group best
imp = head(imp, n=40) #I'm using the first 40 because is the smallest value with the lowest MSE

# imp species name decomposition
str(TaxASV)
imp <- merge(imp, TaxASV, by.x = 0, by.y = 0)
class(imp)
imp = imp[order(imp[,2],decreasing = T),] #the same 40 ASV for the taxonomy composition through time

##PLot the data
imp =read.csv("importance_classRF1-4.csv", na.strings = c("","NA"), header=TRUE)
imp = head(imp, n=40)
imp = imp[order(imp[,2],decreasing = T),]
rownames(imp) <- imp$Row.names

#图4.1. Draw a histogram of the importance of species types
# We are going try to predict the microbiome age using the CON as the training and testing in the SCFP samples
# Make a dataframe of training data with OTUs as column and samples as rows
#predictor using only the 40 ASVs
asvT <- as.data.frame(t(ASV_table))
asvT <- as.matrix(t(asvT[rownames(asvT) %in% rownames(imp),]))

#making a new phyloseq object with the most representative data
#New to generate the OTU for the healthy period
asvT
C_OTU <- asvT[rownames(asvT) %in% rownames(CON),]
S_OTU <- asvT[rownames(asvT) %in% rownames(SCFP),]

##making the phylosq object
OTU.physeq = otu_table(as.matrix(C_OTU), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(TaxASV))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(CON)

#We then merge these into an object of class phyloseq.
physeq_deseq2 = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_deseq2 #[ 40 taxa and 105 samples ]
colnames(tax_table(physeq_deseq2))
predictors <- (otu_table(physeq_deseq2))
dim(predictors) # we have 105  40

# Make one column for our outcome/response variable 
response <- sample_data(physeq_deseq2)$response ##classification
class(response)

# Combine them into 1 data frame
rf.data <- data.frame(response, predictors) ## complete dataframe with the sample ID, classification and OTU table
str(rf.data)
tail(rf.data$response)
#write.csv(rf.data,"OtutableRF.csv") 

### other way
##Randon Forest with the 40 ASV
## performing random forest for the 16S data
trainC <- merge(CON, rf.data, by.x = 0, by.y = 0)
rownames(trainC) <- trainC[,1]
trainC = trainC[,-c(1:3)] ## ASV table with the response
str(trainC)
set.seed(316)
rf = randomForest(response.y ~., trainC, importance=TRUE, proximity=TRUE)
print(rf)
plot(rf)
which.min(rf$mse) #128 is better

rf2 = randomForest(response.y ~., trainC, importance=TRUE, proximity=TRUE, ntree=128)
print(rf2)
plot(rf2)
which.min(rf2$mse)

trainPred <- as.data.frame(rf2[["predicted"]])
trainPred$ID <- rownames(CON)
str(trainPred)  ### Data predicted using the CON 
trainPred <- merge(trainPred, CON, by.x = 0, by.y = 0)
trainPred <- trainPred[-c(3,5)]
trainPred$treatment <- c("CON")
#write.csv(trainPred, "CONprd.csv")


#rfcv：Random Forest Cross Validation
set.seed(316) 
result = rfcv(C_OTU, CON$response, cv.fold=10) #just the ASV table of the 40 ASV and the CON samples
result$error.cv
#plot
p1=with(result, plot(n.var, error.cv, log="x", type="o", lwd=2,xlab="n.var.seed=316"))

#testing data, we are going to use the SCFP data
SCFP
S_OTU <- t(S_OTU)
testS <- merge(imp, S_OTU, by.x = 0, by.y = 0) # merging the 40 ASVs with the CON samples
rownames(testS) <- testS[,1]
testS = t(testS[,-c(1:13)]) ##just the ASV counts
testS = merge(SCFP, testS, by.x = 0, by.y = 0)
rownames(testS) <- testS$Row.names
testS <- testS[,-c(1:10,12)]
testS %>% count(d)

predS = predict(rf2, newdata=testS)#predictions using the 40 ASVs calculated from the CON group
predS
pS <- as.data.frame(predS)
pS$ID <- rownames(pS)
scfpPred <- merge(pS, SCFP, by.x = 0, by.y = 0)
scfpPred  <- scfpPred [-c(3:10,12,14)]
write.csv(scfpPred, "pS.csv")

testS %>%
  mutate(lda.pred = (pS$predS)) %>%
  summarise(lda.error = mean(d != lda.pred))


###Correlation
# we need the predictive values of the CON and SCFP to make a correlation between the actual animal age
# we will use the correlation of repeated measures
install.packages("rmcorr")
library(rmcorr)

corre =read.csv("corre.csv", na.strings = c("","NA"), header=TRUE)
corre$id <- as.factor(corre$id)
corre$treatment <- as.factor(corre$treatment)
corre$d <- as.numeric(corre$d)
str(corre)
corre2 <- corre
corre2$d <- as.factor(corre2$d)

my.rmc <- rmcorr(participant = id, measure1 = pred, measure2 = d, CI.level = 0.95, dataset = corre)
my.rmc
plot(my.rmc, overall = TRUE)



#plot
ggplot(data=corre2, aes(x=d, y=pred, color=treatment)) + 
  geom_boxplot() +
  theme_bw()+
  labs(x='Age (months)', y= 'Predicted-Age (months)') +
  labs(color= "Treatment") +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)) 

#difference between treatment
set_sum_contrasts() # important for afex

# we need to add the withing factor in the error term
M1 <- mixed(pred ~ d + d*treatment + (d||id), data = corre, method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
M1
summary(M1)
anova(M1)
plot(M1$full_model)
qqnorm(residuals(M1$full_model))
qqline(residuals(M1$full_model))

#transform data 
summary(corre$pred)
corre <- mutate(corre, pred_log = log10(pred))
corre <- mutate(corre, pred_sqrt = sqrt(pred))

M1 <- mixed(pred_log ~ d + d*treatment+  (d||id), data = corre, method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
M1
MyNorm <- function(x){ (x-mean(x))/sd(x)}
#Add na.rm = TRUE to deal with NAs

corre$d<- MyNorm(corre$d)
M1 <- mixed(pred_log ~ d + d*treatment + (d||id), data = corre, method = "KR", 
           control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
M1
summary(M1)
anova(M1)
plot(M1$full_model)
qqnorm(residuals(M1$full_model))
qqline(residuals(M1$full_model))

##plotting the data with the transformed pred-Log
corre2 <- mutate(corre2, pred_log = log10(pred))
summary(corre2$pred_log)
corre2$d <- as.factor(corre2$d)

ggplot(data=corre2, aes(x=d, y=pred_log, color=treatment)) + 
  geom_boxplot() +
  theme_bw()+
  labs(x='Age (months)', y= 'Predicted-Age Log10 transformed (months)') +
  labs(color= "Treatment") +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)) 

