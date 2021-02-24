# This is a demo for running the co-occurrence analysis much, much faster

#make sure you have these libraries
#install.packages("zoo")
library(Hmisc)
library(plyr)
library(reshape2)
library(igraph)
library(fdrtool)
library(GGally)
library(intergraph)
library(qiime2R)
library(ggplot2)
library(tidyr)
library(naniar)
library(zoo)

#install.packages("remotes")
#remotes::install_github("jbisanz/qiime2R")

# this is the data
setwd("~/Desktop/eunice/CalfFecal/CalfFecal01_2021/NewQiime/TrimmedFecal/exported/")

####################################################
#I editted here to make it so you don't have to edit the metadata file to remove the second line. 
#Use the same metadata file as you used in QIIME2
####################################################

metadata <- read.delim("CalfFecal_meta.txt", sep = "\t", header = T, quote = "", stringsAsFactors = F)
#metadata <- metadata2[-1,]
str(metadata)
metadata$treatment <- as.factor(metadata$treatment)

ASVs <- read_qza("table.qza")
ASV_s <- as.data.frame(ASVs$data)
ASV_table <- as.data.frame(ASVs$data) #7529 ASVs

#####################################################################
##I added in this chuck. What does it do and why?

ASV_table$ASVnos <- paste0("ASV", 1:nrow(ASV_table))
ASV_table$ASVstring <- rownames(ASV_table)
rownames(ASV_table) <- ASV_table$ASVnos ##We change the ASV name created in Qiime to ASVn
ASVkey <- ASV_table[, (ncol(ASV_table)-1):ncol(ASV_table)] #the key withe the names
ASV_table <- ASV_table[,-(ncol(ASV_table)-1):-ncol(ASV_table)]
######################################################################

##Adding taxonomy
#Taxonomy of each OTU
tax = read.table("taxonomy.tsv", header=TRUE, sep="\t")
tax2 = separate(tax, Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
#This warning means that some cells are empty and that R is replacing empty cells with NA. Now there are other cells that are unclassified that  say, for example `s__` or `g__`. All these need to be removed and replaced with `NA`. 
#All this is OK except that in future use of the taxonomy table, these ASVs will be ignored because they are not classified. Why are ASVs not classified? Its because there is not a close enough match in the database. Just because there is not a good match in the database does not mean they don’t exist, so I wanted to make sure this data was not lost. So in my new code, from lines 300 – 334 I make it so that ASVs that are unclassified at any level are classified as the lowest taxonomic level for which there is a classification.

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
#write.table(tax.final,"taxonomyCalfFecal.txt",sep=",", row.names = FALSE) 

TAX = read.table("taxonomyCalfFecal.txt", header=TRUE, row.names=8, sep="\t")
str(ASVkey)
str(TAX)
TaxASV <- merge(TAX, ASVkey, by.x = 0, by.y = "ASVstring")
#write.table(TaxASV,"TaxASV.txt",sep=",", row.names = FALSE)

dataset <- as.data.frame(t(ASV_table)) #the dataframe is transposed
dataset <- merge(metadata, dataset, by.x = "sample", by.y = 0) #to have all the samples included in the study
datasetn<-dataset
datasetn[datasetn==0]<-NA #change every 0 to NA
# we are going to create a network per treatment
head(dataset[,1:20])
head(datasetn[,1:20])

treatments<-as.vector(unique(dataset$treatment))
treatments
summary(dataset$treatment)
#CON  Pre-CON Pre-SCFP     SCFP 
#106       26       26      107 
#treatments<-treatments[2:5]
#test<-dataset[,-c(1:3)]
#rm(test)
final_results<-data.frame()
i<-4

##This creats the netword of the ASV that are correlated within the same treatment (2 ASVs)
for(i in 1:length(treatments)){
  #subset the data for a particular treatment
  temp<-subset(dataset, treatment==treatments[i])
  tempn<-subset(datasetn, treatment==treatments[i])
  #print(temp$pig)
  # making an object that has all the results in it (both rho and P values)
  results<-rcorr(as.matrix(temp[,-c(1:13)]),type="spearman")
  resultsn<-rcorr(as.matrix(tempn[,-c(1:13)]),type="spearman")
  #make two seperate objects for p-value and correlation coefficients
  rhos<-results$r #correlation value between two ASVs
  ps<-results$P #P-value of the correlation
  ns<-resultsn$n
  # going to melt these objects to 'long form' where the first two columns make up the pairs of OTUs, I am also removing NA's as they are self-comparisons, not enough data, other bad stuff
  ps_melt<-na.omit(melt(ps))
  #creating a qvalue based on FDR
  ps_melt$qval<-p.adjust(ps_melt$value, method="BH")
  #making column names more relevant
  
  names(ps_melt)[3]<-"pval"
  # if you are of the opinion that it is a good idea to subset your network based on adjusted P-values (qval in this case), you can then subset here
  ps_sub<-subset(ps_melt, qval < 0.05)
  # now melting the rhos, note the similarity between ps_melt and rhos_melt
  rhos_melt<-na.omit(melt(rhos))
  names(rhos_melt)[3]<-"rho"
  # now melting the ns
  ns_melt<-(melt(ns))
  names(ns_melt)[3]<-"n"
  #merging together and remove negative rhos
  merged<-merge(ps_sub,rhos_melt,by=c("Var1","Var2"))
  if (treatments[i]=="CON") {
    merged<-merge(merged,subset(ns_melt, n > 24),by=c("Var1","Var2"))
  }   else if (treatments[i]=="SCFP") {
    merged<-merge(merged,subset(ns_melt, n > 24),by=c("Var1","Var2"))
  }   else if (treatments[i]=="Pre-CON") {
    merged<-merge(merged,subset(ns_melt, n > 7),by=c("Var1","Var2"))
  }   else if (treatments[i]=="Pre-SCFP") {
    merged<-merge(merged,subset(ns_melt, n > 7),by=c("Var1","Var2"))
  }   else
    print("Somethings wrong with your treatment designations. Please Check!!")
  merged$treatment<-treatments[i]
  final_results<-rbind(final_results, merged)
  print(paste("finished ",treatments[i],sep=""))
}

##Include taxonomy 
str(TaxASV)
PhylumASV <- (TaxASV[,-c(3,4, 6:8)])
str(PhylumASV)
strong_results<-subset(final_results, abs(rho) >= 0.7)
str(strong_results)
strong_results<-merge(strong_results, PhylumASV, by.x="Var1",by.y="ASVnos")
str(strong_results)
colnames(strong_results) <- c('Var1','Var2', 'pval', 'qval', 'rho', 'n', 'treatment', 'Var1ASV', 'PhylumVar1', 'FamilyVar1')
strong_results<-merge(strong_results, PhylumASV, by.x="Var2",by.y="ASVnos")
colnames(strong_results) <- c('Var2','Var1', 'pval', 'qval', 'rho', 'n', 'treatment', 'Var1ASV', 'PhylumVar1', 'FamilyVar1', 'Var2ASV', 'PhylumVar2', 'FamilyVar2')


#write.csv(subset(strong_results, treatment == "CON"), file="Strong_results_CON_0.7..csv")
#write.csv(subset(strong_results, treatment == "SCFP"), file="Strong_results_SCFP_0.7.csv")
#write.csv(subset(strong_results, treatment == "Pre-CON"), file="Strong_results_PreCON_0.7.csv")
#write.csv(subset(strong_results, treatment == "Pre-SCFP"), file="Strong_results_PreSCFP_0.7.csv")
str(strong_results$treatment)

###ASVs unique betweeen SCFP and CON
#ASVsUnique = read.table("UniqueASVsTRT.txt", header=TRUE, sep="\t")
#ASVsUnique<-merge(ASVsUnique, TaxASV, by.x="ASVs",by.y="ASVnos")
#write.table(ASVsUnique,"ASVsUnique.txt",sep=",", row.names = FALSE)

###ASVs unique betweeen SCFP and CON in each node
ASVsNodes = read.table("UniqueASVnodes.txt", header=TRUE, sep="\t")
ASVsNodes<-merge(ASVsNodes, TaxASV, by.x="ASVs",by.y="ASVnos")
write.table(ASVsNodes,"ASVsNodes.txt",sep=",", row.names = FALSE)

## see which comparison are common in the SCFP and CON
ASVCompCON = read.table("ASVCompCON.txt", header=TRUE, sep="\t")
ASVCompSCFP = read.table("ASVCompSCFP.txt", header=TRUE, sep="\t")
ASVsCommon<-merge(ASVCompCON, ASVCompSCFP, by.x="Comparison",by.y="Comparison")
#write.table(ASVsCommon,"ASVsCommon.txt",sep=",", row.names = FALSE)

##################################################
#Making scatterplots
#
#Make sure your r value is for the treatment group for which you draw these scatterplots

CON_ASVs<-subset(dataset, treatment==treatments[2])
#U_ASVs<-subset(dataset, treatment==treatments[2])

colnames(CON_ASVs[1:10])
head(final_results)

ggplot(CON_ASVs, aes(x = ASV76, y = ASV77)) +
  geom_point()



#I STOPPED EDITTING HERE
#####################################################

head(final_results)
dim(final_results)
hist(final_results$rho)
strong_results<-subset(final_results, rho >= 0.7)
dim(strong_results)
hist(strong_results$rho)






# to make a network from the results
# we need to pass the relevant relationships (pairs across columns Var1 and Var2) to graph.edgelist as a matrix, note I am subsetting for a particular treatment within the arguments
head(final_results)
g<-graph.edgelist(as.matrix(subset(final_results, trt=="carb")[,1:2]),directed=F)

#are relationships different
ggplot(final_results)+geom_density(aes(rho,fill=trt),alpha=0.5)+theme_bw(base_size=17)+theme(aspect.ratio=1)+scale_fill_manual(name="Treatment",values=c("red","black", "green", "gray"))

# now we can calculate stats for the network
final_stats<-data.frame()
for(i in 1:length(unique(final_results$trt))){
  
  temp<-subset(final_results, trt==as.vector(unique(final_results$trt))[i])
  print(head(temp))
  temp.graph<-(graph.edgelist(as.matrix(temp[,c(1,2)]),directed=FALSE))
  E(temp.graph)$weight<-temp$rho
  temp.graph<-simplify(temp.graph)
  stats<-data.frame(row.names((as.matrix(igraph::degree(temp.graph,normalized=TRUE)))),(as.matrix(igraph::degree(temp.graph,normalized=TRUE))),(as.matrix(igraph::betweenness(temp.graph))))
  names(stats)<-c("otus","norm_degree","betweenness")
  stats$trt<-as.vector(unique(final_results$trt))[i]
  stats$clustering_coeff<-igraph::transitivity(temp.graph,type="global")
  stats$clustering_coeff_rand<-igraph::transitivity(igraph::erdos.renyi.game(length(V(temp.graph)),length(E(temp.graph)),type="gnm"))
  stats$cluster_ratio<-stats$clustering_coeff/stats$clustering_coeff_rand
  final_stats<-rbind(final_stats,stats)
  print(paste("finished ",as.vector(unique(final_results$trt))[i],sep=""))
}

# is the structure changing?
ddply(final_stats,.(trt),summarise,mean(cluster_ratio))
for (i in treatments){
  print(head(arrange(subset(final_stats,trt==i), -norm_degree),10))
  print(head(arrange(subset(final_stats,trt==i), -betweenness),10))
}

head(arrange(subset(final_stats,trt=="carb"), -norm_degree),10)
OTU 166 d:Bacteria,p:Firmicutes,c:Negativicutes,o:Selenomonadales,f:Acidaminococcaceae
OTU 2080 d:Bacteria,p:"Bacteroidetes"
OTU 1864 d:Bacteria

head(arrange(subset(final_stats,trt=="non-foaming"), -betweenness),10)
OTU 1864 d:Bacteria
OTU 1601  d:Bacteria,p:"Proteobacteria",c:Betaproteobacteria

# whats the relationship between these statistics?
ggplot(final_stats)+geom_point(aes(x=norm_degree,y=betweenness,color=trt),alpha=0.5)+scale_y_log10()+theme_bw(base_size=17)+labs(x="Normalized Degree",y="Betweenness")+theme(aspect.ratio=1)+scale_colour_manual(name="Treatment",values=c("red","black", "green", "yellow"))




# let's make plots with the most significant relationships
head(final_results)
dim(final_results)
hist(final_results$rho)
strong_results<-subset(final_results, rho >= 0.7)
dim(strong_results)
hist(strong_results$rho)


for (i in treatments){
  print(i)
  temp.graph<-(graph.edgelist(as.matrix(subset(strong_results, trt==i)[,c(1,2)]),directed=FALSE))
  E(temp.graph)$weight<-subset(strong_results, trt==i)$rho
  temp.graph<-simplify(temp.graph)
  temp.graph
  
  gnet<-asNetwork(temp.graph)
  df<-asDF(gnet)
  head(df$vertexes)
  
  
  ggnet(gnet, size=0, method="kamadakawaii")+geom_point()
  
  # lets color by phylum
  vs<-df$vertexes
  phyla<-read.table("~/Desktop/CarbadoxPigs/CAZy/CAZy_co-occurrence/phylum.txt",sep="\t",check.names=F,header=T)
  head(phyla)
  head(vs)
  vs_phyla<-merge(vs, phyla, by.x="vertex.names",by.y="OTUs")
  vs_phyla<-arrange(vs_phyla,intergraph_id)
  print(head(vs_phyla))
  print(ggnet(gnet, label.nodes=T, label.size = 2, size=1, segment.color = "black", method="kamadakawaii")+geom_point(aes(colour=vs_phyla$phylum, fill = i))+scale_colour_manual(values=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","grey50"))+scale_fill_discrete(name=i))
  ggsave(filename = paste("~/Desktop/CarbadoxPigs/CAZy/CAZy_co-occurrence/20160112/",i,"_require5-6_label.pdf", sep=""))
  print(ggnet(gnet, label.nodes=F, label.size = 2, size=1, segment.color = "black", method="kamadakawaii")+geom_point(aes(colour=vs_phyla$phylum, fill = i))+scale_colour_manual(values=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","grey50")))
  ggsave(filename = paste("~/Desktop/CarbadoxPigs/CAZy/CAZy_co-occurrence/20160112/",i,"_require5-6_nolabel.pdf", sep=""))
  print(paste("finished", i))
}

write.csv(strong_results, file="~/Desktop/CarbadoxPigs/CAZy/CAZy_co-occurrence/Strong_results.csv")

library(tidyr)

getwd()
carb_non_unique<-read.table("Strong_results_carb_non_unique.csv",header=T,sep=",",check.names=F)

unique_carb <- subset(carb_non_unique, trt=="carb")
unmelt_carb <- spread(unique_carb, Var1, Var2)

for (i in treatments){
  print(i)
  temp.graph<-(graph.edgelist(as.matrix(subset(carb_non_unique, trt==i)[,c(1,2)]),directed=FALSE))
  E(temp.graph)$weight<-subset(carb_non_unique, trt==i)$rho
  temp.graph<-simplify(temp.graph)
  temp.graph
  
  gnet<-asNetwork(temp.graph)
  df<-asDF(gnet)
  head(df$vertexes)
  
  
  ggnet(gnet, size=0, method="kamadakawaii")+geom_point()
  
  # lets color by phylum
  vs<-df$vertexes
  phyla<-read.table("~/Desktop/CarbadoxPigs/CAZy/CAZy_co-occurrence/phylum.txt",sep="\t",check.names=F,header=T)
  head(phyla)
  head(vs)
  vs_phyla<-merge(vs, phyla, by.x="vertex.names",by.y="OTUs")
  vs_phyla<-arrange(vs_phyla,intergraph_id)
  print(head(vs_phyla))
  print(ggnet(gnet, label.nodes=T, label.size = 2, size=1, segment.color = "black", method="kamadakawaii")+geom_point(aes(colour=vs_phyla$phylum, fill = i))+scale_colour_manual(values=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","grey50"))+scale_fill_discrete(name=i))
  ggsave(filename = paste("~/Desktop/CarbadoxPigs/CAZy/CAZy_co-occurrence/20160112/",i,"_require5-6_carb-non_unique_label.pdf", sep=""))
  print(ggnet(gnet, label.nodes=F, label.size = 2, size=1, segment.color = "black", method="kamadakawaii")+geom_point(aes(colour=vs_phyla$phylum, fill = i))+scale_colour_manual(values=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","grey50")))
  ggsave(filename = paste("~/Desktop/CarbadoxPigs/CAZy/CAZy_co-occurrence/20160112/",i,"_require5-6_carb-non_unique_nolabel.pdf", sep=""))
  print(paste("finished", i))
}

final_strong_stats<-data.frame()
for(i in 1:length(unique(final_results$trt))){
  
  temp<-subset(final_results, trt==as.vector(unique(final_results$trt))[i])
  print(head(temp))
  temp.graph<-(graph.edgelist(as.matrix(temp[,c(1,2)]),directed=FALSE))
  E(temp.graph)$weight<-temp$rho
  temp.graph<-simplify(temp.graph)
  stats<-data.frame(row.names((as.matrix(igraph::degree(temp.graph,normalized=TRUE)))),(as.matrix(igraph::degree(temp.graph,normalized=TRUE))),(as.matrix(igraph::betweenness(temp.graph))))
  names(stats)<-c("otus","norm_degree","betweenness")
  stats$trt<-as.vector(unique(final_results$trt))[i]
  stats$clustering_coeff<-igraph::transitivity(temp.graph,type="global")
  stats$clustering_coeff_rand<-igraph::transitivity(igraph::erdos.renyi.game(length(V(temp.graph)),length(E(temp.graph)),type="gnm"))
  stats$cluster_ratio<-stats$clustering_coeff/stats$clustering_coeff_rand
  final_strong_stats<-rbind(final_strong_stats,stats)
  print(paste("finished ",as.vector(unique(final_results$trt))[i],sep=""))
}

strong_stats_phyla<-merge(final_strong_stats, phyla, by.x="otus", by.y="OTUs")
write.csv(strong_stats_phyla, file="~/Desktop/CarbadoxPigs/Strong_stats_phyla.csv")

for (i in treatments){
  print(head(arrange(subset(final_strong_stats,trt==i), -norm_degree),10))
  print(head(arrange(subset(final_strong_stats,trt==i), -betweenness),10))
}

strong_results_phyla<-merge(strong_results, phyla, by.x="Var1", by.y="OTUs")
strong_results_phyla2<-merge(strong_results_phyla, phyla, by.x="Var2", by.y="OTUs")
write.csv(strong_results_phyla2, file="~/Desktop/CarbadoxPigs/Stronstrong_results_phyla2.csv")


###Old code
temp.graph<-(graph.edgelist(as.matrix(subset(strong_results, trt=="non")[,c(1,2)]),directed=FALSE))
E(temp.graph)$weight<-subset(strong_results, trt=="non")$rho
temp.graph<-simplify(temp.graph)
temp.graph

gnet<-asNetwork(temp.graph)
df<-asDF(gnet)
head(df$vertexes)

ggnet(gnet, size=0, method="kamadakawaii")+geom_point()

# lets color by phylum
vs<-df$vertexes
phyla<-read.table("~/Desktop/otu_phylum.txt",sep="\t",check.names=F,header=T)
head(phyla)
head(vs)
vs_phyla<-merge(vs, phyla, by.x="vertex.names",by.y="otus")
vs_phyla<-arrange(vs_phyla,intergraph_id)
head(vs_phyla)
ggnet(gnet, size=0, method="kamadakawaii")+geom_point(aes(colour=vs_phyla$phylum))+scale_colour_manual(values=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","grey50"))


