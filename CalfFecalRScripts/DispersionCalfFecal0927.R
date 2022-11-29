#Betadisper Nasal Samples
library(vegan) 
library(ggplot2)
library(ggpubr)
library(dplyr) #tally function
library(data.table) #transpose function
library(car) ##for type III SS 
library(lme4)
library(afex)
library(emmeans)

setwd("~/Desktop/eunice/MSc/CalfFecal/CalfFecal01_2021/NewQiime/TrimmedFecal/exported/") #sets new working directory for Windows systems (remember to replace … with your filepath)
calffecal <- read.table("CalfFecal_meta.txt", header=TRUE, row.names=1, sep="\t")


DistWU = read.csv("weighted-distance-matrix.tsv", header=TRUE, sep="\t", stringsAsFactors = F)
row.names(DistWU) <- DistWU[,1]
DistWU = DistWU[,-c(1)]
distWU <- as.dist(DistWU)
#check homogeneity of variance

time <- betadisper(distWU, type = c("centroid"), group =calffecal$d)
time
plot(time)
distances <- time[["distances"]]
distances <- as.data.frame(distances)
distances <- merge(distances, calffecal, by.x = 0, by.y=0) # we are combining the data and removing unneccesary data
distances <- distances[,-c(3, 7:12,14)] # here I need to calculate by hand the mean of the groups
unique(calffecal$calf)
write.csv(distances, "distancesHWU.csv")

#now let's try a different approach which is calculated in correction factor
# you sum all the values for d0, d7 and d14 and divided by the total number of observations in each dat
metadata <- read.csv("betadisperCalfFecal.csv", header=TRUE)
x <- metadata$d0
x <- na.exclude(x) 
mean0 <- mean(x)
x <- sum(x)
nx <- 52 #total observations

y <- metadata$d28
y <- na.exclude(y) 
mean28 <- mean(y)
y <- sum(y)
ny <- 55

z <- metadata$d56
z <- na.exclude(z)
mean56 <- mean(z)
z <- sum(z)
nz <- 54

a <- metadata$d84
a <- na.exclude(a)
mean84 <- mean(a)
a <- sum(a)
na <- 52

b <- metadata$d112
b <- na.exclude(b)
mean112 <- mean(b)
b <- sum(b)
nb <- 52

CF <- (x + y + z + a + b) / (nx + ny + nz + na + nb)

#Now we get the SS of the model
# SScondition = the sum of N * (average distance of each day - the gran mean ) ^ 2
str(metadata)
metadata$average <- as.numeric(metadata$average) # the average is from each animal, not for day
GrandMean <- mean(metadata$average)
summary_table <-metadata %>%
  summarize(
    sqrtDev0 = (mean0 - GrandMean) ^ 2, ###SS of the conditions
    sqrtDev28 = (mean28 - GrandMean) ^ 2,##SS of the conditions
    sqrtDev56 = (mean56 - GrandMean) ^ 2,##SS of the conditions
    sqrtDev84 = (mean84 - GrandMean) ^ 2,##SS of the conditions
    sqrtDev112 = (mean112 - GrandMean) ^ 2,##SS of the conditions
  ) 
#write.csv(summary_table,"summary_table.csv")
##We get the SS of the conditions (time points), you multiple by (n) number of subjects under each (ith) condition
metadata
SSconditions <- 60 * sum (summary_table$sqrtDev0 + summary_table$sqrtDev28  + summary_table$sqrtDev56 +
                            summary_table$sqrtDev84 + summary_table$sqrtDev112)
SSconditions
SSconditions
##Now we get the SS of the within groups which is each individual variarion form the group mean
# Sum of the (individual value - mean of the group) ^ 2
summary_table2 <-metadata %>%
  summarize(
    sqrtDevW0 = (d0 - mean0) ^ 2, ###SS of the conditions
    sqrtDevW28 = (d28 - mean28) ^ 2, ###SS of the conditions
    sqrtDevW56 = (d56 - mean56) ^ 2, ###SS of the conditions
    sqrtDevW84 = (d84 - mean84) ^ 2, ###SS of the conditions
    sqrtDevW112 = (d112 - mean112) ^ 2, ###SS of the conditions
  )
##removing NAs
sqrtDevW0 <- as.data.frame(summary_table2$sqrtDevW0)
sqrtDevW0 <- as.data.frame(sqrtDevW0[-c(7,17,19,22,23,27,28,30),])
colnames(sqrtDevW0) <- c("sqrtDevW0")

sqrtDevW28 <- as.data.frame(summary_table2$sqrtDevW28)
sqrtDevW28 <- as.data.frame(sqrtDevW28[-c(1,27,33,53,57),])
colnames(sqrtDevW28) <- c("sqrtDevW28")

sqrtDevW56 <- as.data.frame(summary_table2$sqrtDevW56)
sqrtDevW56 <- as.data.frame(sqrtDevW56[-c(1,22,25,33,39,57),])
colnames(sqrtDevW56) <- c("sqrtDevW56")

sqrtDevW84 <- as.data.frame(summary_table2$sqrtDevW84)
sqrtDevW84 <- as.data.frame(sqrtDevW84[-c(1,22,25,33,37,39,52,57),])
colnames(sqrtDevW84) <- c("sqrtDevW84")

sqrtDevW112 <- as.data.frame(summary_table2$sqrtDevW112)
sqrtDevW112 <- as.data.frame(sqrtDevW112[-c(1,22,25,33,37,39,52,57),])
colnames(sqrtDevW112) <- c("sqrtDevW112")

SSwithin <- sum (sqrtDevW0$sqrtDevW0 + sqrtDevW28$sqrtDevW28  +  sqrtDevW56$sqrtDevW56 +
                   sqrtDevW84$sqrtDevW84 + sqrtDevW112$sqrtDevW112)

##Now we calculate the SS for each of the subjects
## Number of variables (K) * Sum of (individual average - grand mean) ^2
#K = number of variables: 5 time point
summary_table3 <-metadata %>%
  group_by(calf) %>%
  summarize(
    sqrtDevS = (average - GrandMean) ^ 2, ###SS of the conditions
  )


## Now the we have the square value for all the different time point of each animal, we need to sum, this will give us the SS within 
SSsubject <- 5 * sum (summary_table3$sqrtDevS)

##Now we have the SSwithin and SSmodel, we need to calcualte the SSerror
SSerror <- SSwithin - SSsubject

#Now we specify the degrees of freedom of the model and for the error term
#dfmodel = k - 1, k number of variables
#dferror = (k-1) * (n-1), n is the number of subject
dfmodel = 5 - 1 #we have 5 variable or points
dfSubject = 60 -1
dferror = dfSubject * dfmodel

##Now we process the calculate the MS of the model, MSerror and Fvalue
MSconditions = SSconditions / dfmodel ##SSconditions / (k-1)
MSerror = SSerror / dferror # SSerror / (n-1)(k-1)
Fvalue = MSconditions / MSerror
results <- data.frame(MSconditions,MSerror,Fvalue)

#Calculating pvalue
#we need the degrees of freedom and the F value
p.value <- pf(Fvalue, dfmodel, dferror, lower.tail = FALSE)

#creating the data frame
Response <- c("SSconditions", "SSerror")
SS <- c(SSconditions, SSerror)
DF <- c(dfmodel,dferror)
MS <- c(MSconditions, MSerror)
F.value <- c(Fvalue, "NA")
P <- c(p.value, "NA")

df <- data.frame(Response, SS, DF, MS, F.value, P)
print (df)


##now we'll run a mixed model to find a significance
#it works but data is not normal distributed
set_sum_contrasts() # important for afex
str(distances)
levels(distances$calf) #19
# we need to add the withing factor in the error term

m1 <- mixed(distances ~ d + (d||calf), data = distances,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(m1) ## Significant 
plot(m1$full_model)
qqnorm(residuals(m1$full_model))
qqline(residuals(m1$full_model))
emm_options(lmer.df = "asymptotic") # also possible: 'satterthwaite', 'kenward-roger'
emm_1<- emmeans(m1, "d")
emm_1
update(pairs(emm_1), by = NULL, adjust = "holm")
with(distances, shapiro.test(distances))

#now using as factor time and performing permutest with strata function
time <- betadisper(distWU, type = c("centroid"), group = healthy1$d)
time

metadata <- read.csv("betadisperCalfN.csv", header=TRUE)

time <- betadisper(distWU, type = c("centroid"), group = healthy1$d)
time
aveTime <- print.default(format(tapply(time$distances, time$group, mean))) #getting the average distance from the centroids
mean0 <- 0.4131679
mean7 <- 0.4179798
mean14 <- 0.4353913
GrandMean <- (mean0 + mean7 + mean14) / 3 

##this is a test
metadata$average <- as.numeric(metadata$average)
GrandMean2 <- mean(metadata$average) #calculating the grand mean, it has to be numeric

x <- metadata$d0
x <- na.exclude(x) 
x <- mean(x)

y <- metadata$d7
y <- na.exclude(y) 
y <- mean(y)

z <- metadata$d14
z <- na.exclude(z) 
z <- mean(z)

X <- mean(x,y,z)
GrandMean2
#upload the formatted table
metadata <- read.csv("betadisperCalfN.csv", header=TRUE)
str(metadata)
metadata2 <- na.exclude(metadata) 

#Now we get the SS of the model
# SScondition = the sum of N * (average distance of each day - the gran mean ) ^ 2
str(metadata)
metadata$average <- as.numeric(metadata$average) # the average is from each animal, not for day
summary_table <-metadata %>%
  summarize(
    sqrtDev0 = (mean0 - GrandMean) ^ 2, ###SS of the conditions
    sqrtDev7 = (mean7 - GrandMean) ^ 2,##SS of the conditions
    sqrtDev14 = (mean14 - GrandMean) ^ 2,##SS of the conditions
  ) 
#write.csv(summary_table,"summary_table.csv")
##We get the SS of the conditions (time points), you multiple by (n) number of subjects under each (ith) condition
metadata
s1 <- 12 * summary_table$sqrtDev0
s2 <- 19 * summary_table$sqrtDev7
s3 <- 16 * summary_table$sqrtDev14
SSconditions <- sum(s1,s2,s3)
SSconditions
##Now we get the SS of the within groups which is each individual variarion form the group mean
# Sum of the (individual value - mean of the group) ^ 2
summary_table2 <-metadata %>%
  summarize(
    sqrtDevW0 = (d0 - mean0) ^ 2, ###SS of the conditions
    sqrtDevW7 = (d7 - mean7) ^ 2, ###SS of the conditions
    sqrtDevW14 = (d14 - mean14) ^ 2, ###SS of the conditions
  )
##removing NAs
sqrtDevW0 <- as.data.frame(summary_table2$sqrtDevW0)
sqrtDevW0 <- as.data.frame(sqrtDevW0[-c(8,10,13,14,17:20),])
colnames(sqrtDevW0) <- c("sqrtDevW0")

sqrtDevW7 <- as.data.frame(summary_table2$sqrtDevW7)
sqrtDevW7 <- as.data.frame(sqrtDevW7[-c(12),])
colnames(sqrtDevW7) <- c("sqrtDevW7")

sqrtDevW14 <- as.data.frame(summary_table2$sqrtDevW14)
sqrtDevW14 <- as.data.frame(sqrtDevW14[-c(1,2,10,20),])
colnames(sqrtDevW14) <- c("sqrtDevW14")

SSwithin <- sum (sqrtDevW0$sqrtDevW0 + sqrtDevW7$sqrtDevW7  + 
                   sqrtDevW14$sqrtDevW14)

##Now we calculate the SS for each of the subjects
## Number of variables (K) * Sum of (individual average - grand mean) ^2
#K = number of variables: 5 time point
summary_table3 <-metadata %>%
  group_by(calf) %>%
  summarize(
    sqrtDevS = (average - GrandMean) ^ 2, ###SS of the conditions
  )


## Now the we have the square value for all the different time point of each animal, we need to sum, this will give us the SS within 
SSsubject <- 3 * sum (summary_table3$sqrtDevS)

##Now we have the SSwithin and SSmodel, we need to calcualte the SSerror
SSerror <- SSwithin - SSsubject

#Now we specify the degrees of freedom of the model and for the error term
#dfmodel = k - 1, k number of variables
#dferror = (k-1) * (n-1), n is the number of subject
dfmodel = 3 - 1 #we have 3 variable or points
dfSubject = 20 -1
dferror = dfSubject * dfmodel

##Now we process the calculate the MS of the model, MSerror and Fvalue
MSconditions = SSconditions / dfmodel ##SSconditions / (k-1)
MSerror = SSerror / dferror # SSerror / (n-1)(k-1)
Fvalue = MSconditions / MSerror
results <- data.frame(MSconditions,MSerror,Fvalue)

#Calculating pvalue
#we need the degrees of freedom and the F value
p.value <- pf(Fvalue, dfmodel, dferror, lower.tail = FALSE)

#creating the data frame
Response <- c("SSconditions", "SSerror")
SS <- c(SSconditions, SSerror)
DF <- c(dfmodel,dferror)
MS <- c(MSconditions, MSerror)
F.value <- c(Fvalue, "NA")
P <- c(p.value, "NA")

df <- data.frame(Response, SS, DF, MS, F.value, P)
print (df)

### ----------------------------------Performing a Tuker-Kramer pairwise comparisons
### Performing a Tuker-Kramer pairwise comparisons
#we need the mean of each condition (timepoints)
#Sample size in each condition (N)
# Numerator and denominator degrees of freedom
# Q value (using the Standarized range Q table for the alpha used 0.05)
# Critical Value
# MSE of the model

str(calffecal)
calffecal %>% tally()
calffecal %>% count(d)
#0 52
#28 55
#56 54
#84 52
#112 52
#sample size is uneven 

#we already have the distance means from each group
mean0
mean28
mean56
mean84
mean112

#Let's specify the total number of combinations: k(k-1) / 2 
( 5 * (5-1)) / 2 #10 total combination between the different time points
#Make the table
#"d0-d28","d0-d56","d0-d84", "d0-d112", "d28-d56","d28-d84","d28-d112", "d56-d84",
#"d56-d112", "d84-d112"
d0d28= abs(mean0-mean28)
d0d56= abs(mean0-mean56)
d0d84= abs(mean0-mean84)
d0d112= abs(mean0-mean112)
d28d56= abs(mean28-mean56)
d28d84= abs(mean28-mean84)
d28d112= abs(mean28-mean112)
d56d84= abs(mean56-mean84)
d56d112= abs(mean56-mean112)
d84d112= abs(mean84-mean112)

Comp <- c("d0-d28","d0-d56","d0-d84", "d0-d112", "d28-d56","d28-d84","d28-d112", "d56-d84","d56-d112", "d84-d112")
Abs <- c(d0d28, d0d56,d0d84,d0d112,d28d56,d28d84,d28d112,d56d84,d56d112,d84d112)
data <- data.frame(Comp,Abs)
data
##Now we calculate the critical value 
#we need the degrees of freedom
#Numerator: (5 -1)= 4 #variables
#denominator (k-1)(n-1) = (5-1)(60-1) = 236
Q = 3.685
# Q value: there is no 236 df in the table, so we'll use the closet one: 120, with an alpha of 0.05
calffecal %>% count(d)
#0 52
#28 55
#56 54
#84 52
#112 52
#sample size is uneven 
Abs1 <- Q * sqrt((MSerror / 2) * ((1 / 52) + (1/ 55))) #0-28
Abs2 <- Q * sqrt((MSerror / 2) * ((1 / 52) + (1/ 54))) #0-56
Abs3 <- Q * sqrt(MSerror / 52) # they have the similar number of sample size 0-84
Abs4 <- Q * sqrt(MSerror / 52)# they have the similar number of sample size 0-112
Abs5 <- Q * sqrt((MSerror / 2) * ((1 / 55) + (1/ 54))) # 28-56
Abs6 <- Q * sqrt((MSerror / 2) * ((1 / 55) + (1/ 52))) # 28- 84
Abs7 <- Q * sqrt((MSerror / 2) * ((1 / 55) + (1/ 52))) #28 - 112
Abs8 <- Q * sqrt((MSerror / 2) * ((1 / 54) + (1/ 52))) # 56- 84
Abs9 <- Q * sqrt((MSerror / 2) * ((1 / 54) + (1/ 52))) # 56- 112
Abs10 <- Q * sqrt(MSerror / 52) # 84-112

CD <- c(Abs1, Abs2,Abs3,Abs4,Abs5,Abs6,Abs7,Abs8,Abs9,Abs10)
data <- data.frame(Comp,Abs, CD)
str(data)

data$significance <- ifelse(
  data$Abs > data$CD, 'Sig',
  ifelse(data$Abs < data$CD, 'NS',
         'Mid Change'))
data

