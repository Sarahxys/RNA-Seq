# Differential expression usign EdgeR
I start learning EdgeR through looking at other's script.

# EdgeR DE analysis from Trinity
Below is the whole script. I broke down the script line by line to learn. 
```
library("limma")
library("edgeR")

setwd("C:/Users/Sarah/Desktop/borealis_sexual_antagonism/Differential expression/EdgeR DE")

data = read.table("borealis_liver.counts.matrix", header=T, row.names=1, com='')
col_ordering = c(5,6,7,8,1,2,3,4)
rnaseqMatrix = data[,col_ordering]
#rnaseqMatrix = round(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix)>=2,]
conditions = factor(c(rep("female_liver", 4), rep("male_liver", 4)))

exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
exp_study = calcNormFactors(exp_study)
exp_study = estimateCommonDisp(exp_study)
exp_study = estimateTagwiseDisp(exp_study)
et = exactTest(exp_study, pair=c("female_liver", "male_liver"))
tTags = topTags(et,n=NULL)

head(tTags$table)
sum(tTags$table$FDR < .1)
plotSmear(et, de.tags = rownames(tTags)[tTags$table$FDR < .01])
abline(h = c(-2, 2), col = "blue")


write.table(tTags, file='liver.counts.matrix.female_liver_vs_male_liver.edgeR.DE_results', sep='	', quote=F, row.names=T)

dev.off()
```
Load in the R package and set up the session directory on my laptop.
```
library("limma")
library("edgeR")

setwd("C:/Users/Sarah/Desktop/borealis_sexual_antagonism/Differential expression/EdgeR DE")
```
Here I read in the data file that is in the working session directory. There is a header so `header = T`. The first row is the transcript ID so I set the first row to be the row.name 
```
data = read.table("borealis_liver.counts.matrix", header=T, sept = "\t",row.names=1, com='')
```



#EdgeR DE analysis from 
