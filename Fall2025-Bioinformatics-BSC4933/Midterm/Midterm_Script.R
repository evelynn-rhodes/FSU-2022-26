if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

library(edgeR)
#edge R manual: https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
library(dplyr)
library(ggplot2)
install.packages("tidyverse")
library(tidyr)

CountTable1<-read.csv("GSE272731_Data_CountMatrix.csv",header = T)

head(CountTable1)
dim(CountTable1)

rownames(CountTable1)<-CountTable1$gene_id
head(CountTable1)

str(CountTable1)

groupings<-c("0.1-uM-dist", "0.1-uM-dist", "0.1-uM-dist", "1-uM-dist", "1-uM-dist", "1-uM-dist", "5-uM-dist", "5-uM-dist", "5-uM-dist", "DMSO-dist", "DMSO-dist", "DMSO-dist","DMSO-prox","DMSO-prox", "DMSO-prox")
#groupings<-c("GC1171", "GC1171", "GC1171", "GC1171", "GC1459", "GC1459", "GC1459", "GC1459")

#Make a DGEList object in edgeR so that you can do all further processing!
d1<-DGEList(counts=CountTable1,group=factor(groupings))


#Check out the components of your DGEList object
d1$samples

head(d1$genes)

head(d1$counts)

d1<-calcNormFactors(d1)

d1$samples #Now this has the norm. factors column calculated

dim(d1) #check the dimensions of the d2 without any filtering
keep_filter<-rowSums(cpm(d1)>1)>=4 #decide how you want to filter the data
d1<-d1[keep_filter,] #filter d2 to only include genes that passed filtering
dim(d1) #check the new dimensions

d1<-estimateDisp(d1) #calculate dispersion (uses negative binomial model)

MDS_plot<-plotMDS(d1)

plotBCV(d1)

de.tag<-exactTest(d1,d1$tagwise.dispersion,pair=c("DMSO-dist","DMSO-prox" ))
#de.tag<-exactTest(d3,d3$tagwise.dispersion,pair=c("GC1171","GC1459"))
head(de.tag$table)

#Now, what are the components of your de.tag object?
#LogFC, logPCM, PValue
head(de.tag$table)

head(de.tag$genes)

de.tag$comparison

str(de.tag)

#Top 10 LogFC
topTags(de.tag) #Quick way to see your top differentially expressed genes

#We want to see all the genes in a nice table...
de.tag_sort<-topTags(de.tag,n=nrow(de.tag$table))$table

head(de.tag_sort)

#Now we can filter this to significant differentially expressed genes

de.tag_sort_sig<-de.tag_sort%>%
  filter(FDR<0.05)

dim(de.tag_sort_sig) #over half of the genes!

#Now, let's say we wanted to look at the expression of some of our most significant genes
#First, we need a table of the normalized CPMs (NOT COUNTS!)

cpm_table<-cpm(d1,normalized.lib.sizes = TRUE)
#Normally will make it a matrix, force it as a df
cpm_table<-as.data.frame(cpm_table)
cpm_table$gene_id<-rownames(cpm_table)
head(cpm_table)
dim(cpm_table)

rownames(CountTable1) <- CountTable1$gene_id

# CHANGE THIS STEP

#Next, we need to put it in long format and add the condition groups as a column
cpm_table_long<-cpm_table%>%
  pivot_longer(cols = 1:15,names_to = "condition",values_to = "CPMs")%>%
  mutate(groups=rep(groupings,39781))#The number is how many genes there are after filtering

#this takes the groupings vector and replicates it the number listed times, which is the number of genes in the table

print(cpm_table_long,n=15) #looks good

#Now, we can use this to make a graph of the top tag genes

top10genes<-de.tag_sort[1:10,]
dim(top10genes)
top10geneNames<-top10genes$gene_id
top10geneNames
print(top10geneNames)

cpm_table_long%>%
  filter(gene_id %in% top10geneNames)%>%
  ggplot(aes(x=groups,y=CPMs,color=gene_id))+
  geom_boxplot()+geom_point()+facet_wrap(vars(gene_id),ncol=5,scales = "free_y")
#in all cases except 1, the genes are going up in L2 larvae and essentially not expressed in L1 larvae.
#these are likely developmental genes turning on

#If you want to export a table of results, run the following to generate a text file. Then you can open in excel and save as whatever you want.
write.table(de.tag_sort,"L1vsL2_RNAseq_Results.txt",quote = F,sep = "\t")
#this will save to your working directory (likely your downloads folder, unless you changed it)
#quote=F prevents R from putting quotes around everything. sep="\t" creates a tab separated file.