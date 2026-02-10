rm(list=ls())
setwd("~/Desktop/")
lapply(c("DESeq2", "apeglm"), library, character.only=T)

# NDC sample YOUR CONTROL SAMPLE HERE
sample1 <- read.csv("Manuscript/DEseq2output/FeatureCountsSplit/Processed/TetNDC_T90_BothStrand_counts.csv", header = T, row.names = "Geneid")
# Mutant sample YOUR TREATMENT SAMPLE HERE
sample2 <- read.csv("Manuscript/DEseq2output/FeatureCountsSplit/Processed/Tet1Q_T90_BothStrand_counts.csv", header = T, row.names = "Geneid")

#Create CountData
countdata <- sample1[,c((ncol(sample1)-2):ncol(sample1))]
countdata[,c(4:6)] <- sample2[,c((ncol(sample2)-2):ncol(sample2))]

colnames(countdata) <- c("wtc-a", "wtc-b", "wtc-c", "treatment-a", "treatment-b", "treatment-c")
coldata <- read.csv("Manuscript/DEseq2output/coldata.csv", header = T, stringsAsFactors = F)
rownames(countdata) <- rownames(sample1)
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design= ~ condition)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)
res

res <- lfcShrink(dds, coef="condition_tr_vs_ctr", type="apeglm")

res
write.csv(as.data.frame(res), file="~/Desktop/Manuscript/DEseq2output/Tet1Q_vs_NDC_t90_split.csv")