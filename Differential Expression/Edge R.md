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
# Broke down of the script
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
Here a vector is created to store the order of samples. 5-8 is a group and 1-4 is a group, this vector will be useful when we do grouping later
```
col_ordering = c(5,6,7,8,1,2,3,4)
```
Here we re-order the data frame 'data' using the 'col_ordering' vector. The store the re-ordered data frame into another data frame. 
```
rnaseqMatrix = data[,col_ordering]
```
Here we filter the matrix and only retain the row that its row sum is greater than 1. Each row contain expression read count for each transcript, if the read count is smaller than 1, the transcript is probably is not a real transcript or is not expressing. 
```
rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix)>=1,]
```
Here we did two things. We first created a vector that contain "female_liver" 4 times and "male_liver" 4 times, which corresponds to the number number of samples that I have for each condition (female vs male). Then the function factor() is used to convert the vector into factor (categorized vectors). So in this case, "female_liver" will be category 1 and "male_liver" will be category 2 because f comes before m.  In R’s memory, these factors are represented by numbers (1, 2, 3). **I dont fully understand factor but it seems like it is useful in grouping** 
```
conditions = factor(c(rep("female_liver", 4), rep("male_liver", 4)))
```
The function DGEList() creates a DGEList object from a table/matrix of read counts (rows=features, columns=samples), 
group indicator which is a vector or factor giving the experimental group/condition for each sample/library.
```
exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
```

```
exp_study = calcNormFactors(exp_study)
```


# EdgeR DE analysis from 
