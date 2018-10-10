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
The function calcNormFactors in EdgeR calculate normalization factors to scale the raw library sizes. The default method for normalization is TMM, which is the weighted trimmed mean of M-values (to the reference) proposed by Robinson and Oshlack (2010), where the weights are from the delta method on Binomial data. If refColumn is unspecified, the library whose upper quartile is closest to the mean upper quartile is used. If the input object is a DGELIST, the DGELIST, then it is returned as output with the relative normalization factor in object$samples$norm.factors. 
  - note: rows that have zero counts for all comuns are trimmed before normalization factors are computed. Therefore rows with a zero counts do not affect the estimated factos. 
```
exp_study = calcNormFactors(exp_study)
```
the function estimateCommonDisp maximizes the negative binomial conditional common likelihood (CML) to estimate a common dispersion value across all genes. The CML method involves computing a matrix of quantile-quantile normalized counts, called pseudo-counts. The pseudo-counts are adjusted in such a way that the library sizes are equal for all samples, while preserving differences between groups and variability within each group. The pseudo-counts are included in the output of the function, but are intended mainly for internal edgeR use.
```
exp_study = estimateCommonDisp(exp_study)
```


# EdgeR DE analysis from online tutorial

# Compare the two script
### Trinity script
1. read in the raw read count matrix
2. re-order the columns based on the orders of sample grouping (female vs male; column 5, 6, 7, 8 have read count for female and column 1, 2,3 4 for male)
3. filter the matrix: it eleminate rows that have a row sum less than 1, 
    - have a row sum less than 1 means the total number of raw read for both male and female is less than one.
    - this could be reasonable since if the total number of supporting read for a transcript is less than 1, it might not be a real transcript. 
    - meanwhile, it could be un-reasonable since the transcirpt is one of the isoforms of a gene and the gene have a relatively low expression level, the "asigned" or distributed read to that transcript might be less than 1. 
    - as for now, I would said this a reasonable filter
4. creating a condition-grouping vector and use the vector to create a DEGList object from the table of raw read counts (rows = gene/transcript id, columns = samples)
5. calculate the normalizeing factor for between sample normalization
    - edgeR use TMM (weighted trimmed mean of m-values) normalization to account for the difference in sequencing depth between samples
    - the output of 
6. estimate common dispertion
    - 


### Tutorial script
1. read in the raw read count matrix
