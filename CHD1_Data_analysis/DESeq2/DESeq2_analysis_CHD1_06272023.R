library(DESeq2)
library(dplyr)
library(ggplot2)
library(AnnotationDbi)
library(org.Dm.eg.db)

setwd("/nas/longleaf/home/bdmcmi/work/Benjamin/flyHead_analysis/CHD1_Data_analysis/DESeq2")
currentDate <- Sys.Date() %>% format("%Y%m%d")

fullCounts <- read.table("/nas/longleaf/home/bdmcmi/work/Benjamin/flyHead_analysis/CHD1_Data_analysis/featureCounts/sortedByCoord_counts.txt",
                     header = TRUE, row.names = 1)

counts <- fullCounts[,c(6:14)] %>%
            dplyr::rename("SRR11237854_WT_1" = "SRR11237854_Aligned.sortedByCoord.out.bam",
                          "SRR11237855_WT_2" = "SRR11237855_Aligned.sortedByCoord.out.bam",
                          "SRR11237856_WT_3" = "SRR11237856_Aligned.sortedByCoord.out.bam",
                          "SRR11237857_Chd1Def_1" = "SRR11237857_Aligned.sortedByCoord.out.bam",
                          "SRR11237857_Chd1Def_2" = "SRR11237858_Aligned.sortedByCoord.out.bam",
                          "SRR11237857_Chd1Def_3" = "SRR11237859_Aligned.sortedByCoord.out.bam",
                          "SRR11237857_Chd1Rescue_1" = "SRR11237860_Aligned.sortedByCoord.out.bam",
                          "SRR11237857_Chd1Rescue_2" = "SRR11237861_Aligned.sortedByCoord.out.bam",
                          "SRR11237857_Chd1Rescue_3" = "SRR11237862_Aligned.sortedByCoord.out.bam"
                          )

replicate <- factor(rep(c(1, 2, 3), 3))

genotype <- factor(c(rep("WT", 3),
              rep("Chd1Def", 3),
              rep("Chd1Rescue", 3)
                  ))

sampleInfo <- data.frame(replicate, genotype)
rownames(sampleInfo) <- colnames(counts)

rownames(sampleInfo) 
colnames(counts)
all(rownames(sampleInfo) == colnames(counts))

dds <- DESeqDataSetFromMatrix(countData = counts, colData = sampleInfo, design = ~ genotype)
dds$genotype = relevel(dds$genotype, "WT")
dds <- DESeq(dds)

normCounts <- counts(dds, normalized=T)
write.csv(as.data.frame(normCounts), 
          file = paste0("DESeq2_",currentDate,"_normCounts.csv"))

# Plot PCA
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup = "genotype")

pcaData <- plotPCA(vsd, intgroup="genotype", returnData=TRUE)
percentVar <- round(100* attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"%variance"))+
  ylab(paste0("PC2: ",percentVar[2],"%variance"))+
  ggtitle("Chd1 Genotypes") +
  geom_jitter()

Chd1Def_res <- results(dds, contrast = c("genotype", "Chd1Def", "WT"))
summary(Chd1Def_res, 0.05)

# MA plot
plotMA(Chd1Def_res, ylim = c(-5, 5))

## READ THIS BEFORE RUNNING
## Added on 06-29-2023 simply to count number of up and down enriched peaks that correspond to Venn diagrams!!!

Chd1Def_sigCount <- Chd1Def_res %>%
  data.frame() %>%
  dplyr::filter(padj < 0.05) %>% 
  #dplyr::filter(log2FoldChange > 1 | log2FoldChange < -1) %>%
  dplyr::mutate(Enrichment = ifelse(log2FoldChange > 0, 'Up', 'Down'))
table(Chd1Def_sigCount['Enrichment'])

Chd1Def_res %>%
  data.frame() %>%
  dplyr::filter(padj < 0.05) %>% 
  #dplyr::filter(log2FoldChange > 1 | log2FoldChange < -1) %>%
  dplyr::mutate(Enrichment = ifelse(log2FoldChange > 0, 'Up', 'Down')) %>%
  write.csv(., file = paste0("DESeq2_",currentDate,"_Chd1Def_sigCounts.csv")
  )

# This adds extra columns with gene id and descriptions
Chd1Def_res$symbol <- mapIds(org.Dm.eg.db, keys=row.names(Chd1Def_res), column="SYMBOL", keytype="FLYBASE", multiVals="first")
Chd1Def_res$name <- mapIds(org.Dm.eg.db, keys=row.names(Chd1Def_res), column="GENENAME", keytype="FLYBASE", multiVals="first")

Chd1Def_resOrdered <- Chd1Def_res[order(Chd1Def_res$padj),]
write.csv(as.data.frame(Chd1Def_resOrdered), 
          file = paste0("DESeq2_",currentDate,"_Chd1Def_res.csv"))

Chd1Rescue_res <- results(dds, contrast = c("genotype", "Chd1Rescue", "WT"))
summary(Chd1Rescue_res, 0.05)

# MA plot
plotMA(Chd1Rescue_res, ylim = c(-5, 5))

# This adds extra columns with gene id and descriptions
Chd1Rescue_res$symbol <- mapIds(org.Dm.eg.db, keys=row.names(Chd1Rescue_res), column="SYMBOL", keytype="FLYBASE", multiVals="first")
Chd1Rescue_res$name <- mapIds(org.Dm.eg.db, keys=row.names(Chd1Rescue_res), column="GENENAME", keytype="FLYBASE", multiVals="first")

Chd1Rescue_resOrdered <- Chd1Rescue_res[order(Chd1Rescue_res$padj),]
write.csv(as.data.frame(Chd1Rescue_resOrdered), 
          file = paste0("DESeq2_",currentDate,"_Chd1Rescue_res.csv"))


Chd1Rescue_Chd1Def_res <- results(dds, contrast = c("genotype", "Chd1Rescue", "Chd1Def"))
summary(Chd1Rescue_Chd1Def_res, 0.05)

# MA plot
plotMA(Chd1Rescue_Chd1Def_res, ylim = c(-5, 5))

# This adds extra columns with gene id and descriptions
Chd1Rescue_Chd1Def_res$symbol <- mapIds(org.Dm.eg.db, keys=row.names(Chd1Rescue_Chd1Def_res), column="SYMBOL", keytype="FLYBASE", multiVals="first")
Chd1Rescue_Chd1Def_res$name <- mapIds(org.Dm.eg.db, keys=row.names(Chd1Rescue_Chd1Def_res), column="GENENAME", keytype="FLYBASE", multiVals="first")

Chd1Rescue_Chd1Def_resOrdered <- Chd1Rescue_Chd1Def_res[order(Chd1Rescue_Chd1Def_res$padj),]
write.csv(as.data.frame(Chd1Rescue_Chd1Def_resOrdered), 
          file = paste0("DESeq2_",currentDate,"_Chd1Rescue_Chd1Def_res.csv"))
