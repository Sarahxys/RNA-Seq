

```
library(edgeR)
library(DESeq2)

data = read.table("/home/xue/borealis_transcriptome/borealis_denovo_transcriptome_august2017/sum_expression/borealis_collapsedT_DE/borealis_expression_per")
col_ordering = c(1,2,3,4,5,6,7,8)
rnaseqMatrix = data[,col_ordering]
#rnaseqMatrix = round(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]
conditions = data.frame(conditions=factor(c(rep("female", 4), rep("male", 4))))
rownames(conditions) = colnames(rnaseqMatrix)

#store it in the DESeq DE object
ddsFullCountTable <- DESeqDataSetFromMatrix(
    countData = rnaseqMatrix,
    colData = conditions,
    design = ~ conditions)

dds = DESeq(ddsFullCountTable)
contrast=c("conditions","female","male")
res = results(dds, contrast)
baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$condition == "female"])
baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$condition == "male"])
res = cbind(baseMeanA, baseMeanB, as.data.frame(res))
res = cbind(sampleA="female", sampleB="male", as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1
write.table(as.data.frame(res[order(res$pvalue),]), file='borealis_expression_per_gene.matrix.female_vs_male.DESeq2.DE_results', sep='  ', quote=FALSE)

```

