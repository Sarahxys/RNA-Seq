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

write.table(tTags, file='liver.counts.matrix.female_liver_vs_male_liver.edgeR.DE_results', sep='	', quote=F, row.names=T)
```
Below is the step involved in this script:
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
7. estimate empirical Bayes Tagwise dispersion values
8. perform comparison between the two group of values = differential expression between the two condition groups
9. extract the useful information from the differential expression analysis
10. output the DE result

## Broke down of the script
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
Here we do some basic filtering before normalization. Trinity v2 edgeR script filters the matrix and only retain the row that its row sum is greater than 1. Each row contain expression read count for each transcript, if the read count is smaller than 1, the transcript is probably is not a real transcript or is not expressing. Trinity v4 edgeR script filters the matrix and only retain the row that its cpm-normalized count greater than 1 and the row sum of cmp-normalized count greater than 2. 
```
#v2
rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix)>=1,]
#v4
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]
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
The function calcNormFactors in EdgeR calculate normalization factors to scale the raw library sizes. This is inter-sample normalization to account for the differnce in library size. The default method for normalization is TMM, which is the weighted trimmed mean of M-values (to the reference) proposed by Robinson and Oshlack (2010), where the weights are from the delta method on Binomial data. If refColumn is unspecified, the library whose upper quartile is closest to the mean upper quartile is used. If the input object is a DGELIST, the DGELIST, then it is returned as output with the relative normalization factor in object$samples$norm.factors. 
  - note: rows that have zero counts for all comuns are trimmed before normalization factors are computed. Therefore rows with a zero counts do not affect the estimated factos. 
```
exp_study = calcNormFactors(exp_study)
```
**I think this is a intra-sample normalization, normalization between genes of a samples (need to be confirmed).** the function estimateCommonDisp maximizes the negative binomial conditional common likelihood (CML) to estimate a common dispersion value across all genes. The CML method involves computing a matrix of quantile-quantile normalized counts, called pseudo-counts. The pseudo-counts are adjusted in such a way that the library sizes are equal for all samples, while preserving differences between groups and variability within each group. The pseudo-counts are included in the output of the function, but are intended mainly for internal edgeR use.
  - if the input object is a DEGLIST, the output add the following components to the input DGEList object:
    - common.dispersion: estimate of the common dispersion
    - pseudo.counts: numeric matrix of pseudo-counts
    - pseudo.lib.size: the common library size to which the pseudo-count have been adjusted
    - AveLogCPM: numeric vector giving log2(AveCPM) for each row of y
```
exp_study = estimateCommonDisp(exp_study)
```
**I dont fully understand the why and how of this step** I think this a continuation of the intra-sample normalization. The function estimateTagwiseDisp estimate tagwise negative binomial dispersion values by an empirical Bays method based on weighted conditional maximum likelihood.
  - The prior values for the dispersions are determined by a global trend. The individual tagwise dispersions are then squeezed towards this trend. The prior degrees of freedom determines the weight given to the prior. The larger the prior degrees of freedom, the more the tagwise dispersions are squeezed towards the global trend. If the number of libraries is large, the prior becomes less important and the tagwise dispersion are determined more by the individual tagwise data.
  - it returned the following component to the input DGElist object:
    - prior.df: prior degrees of freedom.
    - prior.n: estimate of the prior weight.
    - tagwise.dispersion: numeric vector of the tagwise dispersion estimates.
    - span: width of the smoothing window, in terms of proportion of the data set.
```
exp_study = estimateTagwiseDisp(exp_study)

```
This is where the comparison actually happend. The function exactTest tests for differential expression between two groups of count libraries. They implement the exact test proposed by Robinson and Smyth (2008) for a difference in mean between two groups of negative binomial random variables. The functions accept two groups of count libraries, and a test is performed for each row of data. For each row, the test is conditional on the sum of counts for that row. The test can be viewed as a generalization of the well-known exact binomial test (implemented in binomTest) but generalized to overdispersed counts.
  - exactTest produces an object of class DGEExact containing the following components:
    - table: data frame containing columns for the log2-fold-change, logFC, the average log2-counts-per-million, logCPM, and the two-sided p-value PValue
    - comparison: character vector giving the names of the two groups being compared
    - genes: optional data frame containing annotation for each gene; taken from object

```
et = exactTest(exp_study, pair=c("female_liver", "male_liver"))

```
This step is to extract the useful information from the differential expression analysis. The function topTags (Table Of The Top Differentially Expressed Tags) extracts the top DE tags in a data frame for a given pair of groups, ranked by p-value or absolute log-fold change.
  - its default usage is: topTags(object, n=10, adjust.method="BH", sort.by="PValue", p.value=1)
  - n controls the number of tags to display/return. In the trinity script, the limit of number of tags to display is set to null, meaning no limit was set and all tags were returns or display.
  - an object of class TopTags containing the following elements for the top n most differentially expressed tags as determined by `sort.by`:
    - table: a data frame containing the elements 1) logFC, the log-abundance ratio, i.e. fold change, for each tag in the two groups being compared, 2) logCPM, the log-average concentration/abundance for each tag in the two groups being compared, 3) PValue, exact p-value for differential expression using the NB model, 4) FDR, the p-value adjusted for multiple testing as found using p.adjust using the method specified.
    - adjust.method: character string stating the method used to adjust p-values for multiple testing.
    - comparison: a vector giving the names of the two groups being compared.
    - test: character string stating the name of the test.
```
tTags = topTags(et,n=NULL)
```
Output the result of DE into a tsv file
```
write.table(tTags, file='liver.counts.matrix.female_liver_vs_male_liver.edgeR.DE_results', sep='	', quote=F, row.names=T)
```

# EdgeR DE analysis from online tutorial 1
This script is from https://gist.github.com/jdblischak/11384914/a4b57e05fd77a3cd1012977662d7b0b31158dc8f 
```
library("limma")
library("edgeR")

setwd("D:/school_grad school/Project_borealis_sexual_antagonism/documentation/Differential expression/EdgeR DE")

data_raw <- read.table("borealis_liver.counts.matrix", header = TRUE)

cpm_log <- cpm(data_raw, log = TRUE) 
median_log2_cpm <- apply(cpm_log, 1, median)
expr_cutoff <- -3
data_clean <- data_raw[median_log2_cpm > expr_cutoff, ]

group <- substr(colnames(data_clean), 1, 1)
y <- DGEList(counts = data_clean, group = group)
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
results_edgeR <- topTags(et, n = nrow(data_clean), sort.by = "none")

```
## Broke down of the script
Load in the R package and set up the session directory on my laptop.
```
library("limma")
library("edgeR")

setwd("D:/school_grad school/Project_borealis_sexual_antagonism/documentation/Differential expression/EdgeR DE")
```
Here I read in the data file that is in the working session directory. There is a header so `header = T`. The first row is the transcript ID so I set the first row to be the row.name 
```
data_raw <- read.table("borealis_liver.counts.matrix", header = TRUE)
```
We should also remove genes that are unexpressed or very lowly expressed in the samples. One simple method to do this is to choose a cutoff based on the median log2-transformed counts per gene per million mapped reads (cpm). edgeR provides the function, cpm, to compute the counts per million.
```
cpm_log <- cpm(data_raw, log = TRUE) 
median_log2_cpm <- apply(cpm_log, 1, median)
expr_cutoff <- -3
data_clean <- data_raw[median_log2_cpm > expr_cutoff, ]#removing all genes with a median log2 cpm below expr_cutoff
```
Do normalization like what the Trinity edgeR script did. Inter sample normalization with `calNormFactors`, intra-sample normalization with `estimateDisp` , then do DE with `exactTest`, and finally extract result using `topTags`. 
```
group <- substr(colnames(data_clean), 1, 1)
y <- DGEList(counts = data_clean, group = group)
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
results_edgeR <- topTags(et, n = nrow(data_clean), sort.by = "none")
```
# EdgeR DE analysis from online tutorial 2
This script is from http://www.nathalievialaneix.eu/doc/html/solution_edgeR-tomato-withcode.html
```
rawCountTable <- read.table("countData.txt", header=TRUE, sep="\t", row.names=1)
sampleInfo <- read.table("design.csv", header=TRUE, sep=",", row.names=1)

dgeFull <- DGEList(rawCountTable, group=sampleInfo$condition)

#remove genes with zero counts for all samples
dgeFull <- DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) != 0, ],
                   group=dgeFull$samples$group)
#normalization
dgeFull <- calcNormFactors(dgeFull, method="TMM")
dgeFull <- estimateCommonDisp(dgeFull)
dgeFull <- estimateTagwiseDisp(dgeFull)
dgeTest <- exactTest(dgeFull)
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table))
sum(resNoFilt$table$FDR < 0.01)

#extract and sort differentially expressed genes
sigDownReg <- resNoFilt$table[resNoFilt$table$FDR<0.01,]
sigDownReg <- sigDownReg[order(sigDownReg$logFC),]
sigUpReg <- sigDownReg[order(sigDownReg$logFC, decreasing=TRUE),]
write.csv(sigDownReg, file="sigDownReg_tomato.csv")
write.csv(sigUpReg, file="sigUpReg_tomato.csv")

#create a MA plot with 1% differentially expressed genes
plotSmear(dgeTestFilt,
          de.tags = rownames(resFilt$table)[which(resFilt$table$FDR<0.01)])
```

# General pipeline of edgeR DE 
After comparing the example scripts, below is the general pipeline that I conclude: 

1. Read in the raw read count matrix
2. Filter data before normalization
  - filter the matrix: 
    - Trinity script v2: row sum of raw count greater than 1
    - Trinity script v4: cpm normalized count greater than 1; then row sum of cmp normalized count greater than 2. 
    - Tutorial 1 script: calculate cpm_log and find the median cpm_log for each row; the median cpm_log need to be greater than a set cutoff
    - Tutorial 2 script: remove genes with zero counts for all samples
  
3.  Grouping based on experimental condition
    - re-order the columns based on the orders of sample grouping
    - creating a condition-grouping vector and use the vector to create a DEGList object from the table of raw read counts (rows = gene/transcript id, columns = samples)
4. Normalization
    - inter-sample normalization: calculate the normalizeing factor for between sample normalization
      - edgeR use TMM (weighted trimmed mean of m-values) normalization to account for the difference in sequencing depth between samples
    - intra-sample normalization: estimate common dispertion and estimate empirical Bayes Tagwise dispersion values
5. Differential expression
    - cperform comparison between the two group of values = differential expression between the two condition groups
6. Extract the useful information from the differential expression analysis
7. Output the DE result


