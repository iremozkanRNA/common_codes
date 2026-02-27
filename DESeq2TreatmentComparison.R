rm(list=ls())
setwd("~/Desktop/")
lapply(c("DESeq2", "apeglm"), library, character.only=T)

# NDC sample your CONTROL SAMPLE HERE
sample1 <- read.table("Mutant/RNAseq/WTK/bam_output/AB/WtKChl_s2.txt", header = T, row.names = "Geneid")
# Mutant sample your TREATMENT SAMPLE HERE
sample2 <- read.table("Mutant/RNAseq/M1/bam_output/AB/MutantChl_s2.txt", header = T, row.names = "Geneid")

#Create CountData
countdata <- sample1[,c((ncol(sample1)-2):ncol(sample1))]
countdata[,c(4:6)] <- sample2[,c((ncol(sample2)-2):ncol(sample2))]
colnames(countdata) <- c("control-a", "control-b", "control-c", "treatment-a", "treatment-b", "treatment-c")

#Read your coldata (Put control first - IMPORTANT) needs to be alphabetical
coldata <- read.csv("Mutant/RNAseq/coldata.csv", header = T, stringsAsFactors = F)
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
write.csv(as.data.frame(res), file="~/Desktop/WtK_vs_M1_Chl_s2.csv")
