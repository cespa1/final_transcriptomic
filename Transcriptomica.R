##Para realizar esta parte del código se uso el ambiente conda facilitado en las prácticas 
##conda activate dge_arrays
##Se carga las librerias
library("DESeq2")
library("pheatmap")
library("RColorBrewer")


##Establecer el directorio Apartado2
setwd("/home/vant/transcriptomic-final-exercise/Apartado2/")

##Se carga los inputs
rawcounts <- as.matrix(read.table("input/rawcounts.tsv", header = TRUE, row.names = 1, sep = "\t"))
metadata <- read.table("input/metadata.tsv", header = TRUE, row.names =  1, sep = "\t")

metadata24 <- t(metadata$time == "24h")
metadata48 <-  subset(metadata, metadata$time=="48h")
test <- subset(rawcounts,metadata$time=="48h")

dds <- DESeqDataSetFromMatrix(countData = rawcounts,
                              colData = metadata,
                              design= ~ agent + time)

keep <- rowSums(counts(dds)) >= 10 ## Seleccionar genes con más de 10 cuentas en todos los samples
dds <- dds[keep, ]

vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = c("time","agent"))
plotPCA(vsd, intgroup = "agent")

sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$time, vsd$agent, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

dds2 <- DESeq(dds, test = "Wald")

plotDispEsts(dds2)

plotMA(dds2)

my_results2 <- results(object = dds2,
                      contrast = c("agent", "DPN", "Control"),
                      lfcThreshold = 1,
                      alpha = 0.05,
                      pAdjustMethod = "BH",
                      tidy = TRUE
)

mat2 <- assay(vsd)[head(order(my_results2$padj), 30), ] 
pheatmap(mat2)

select <- order(rowMeans(counts(dds2,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("agent","time")])
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

dds3 <- dds2[metadata24, ]

plotDispEsts(dds3)

plotMA(dds3)

my_results <- results(object = dds3,
                      contrast = c("agent", "Control", "OHT"),
                      lfcThreshold = 1,
                      alpha = 0.05,
                      pAdjustMethod = "BH",
                      tidy = TRUE
)

my_results2 <- results(object = dds2,
                       contrast = c("agent", "Control", "OHT"),
                       lfcThreshold = 1,
                       alpha = 0.05,
                       pAdjustMethod = "BH",
                       tidy = TRUE
)

mat <- assay(vsd)[head(order(my_results$padj), 30), ] 
pheatmap(mat)




dds_gsea <- DESeq(dds)
res <- results(dds_gsea, alpha = 0.05, contrast = c("agent", "DPN", "Control"))


res.ape <- lfcShrink(dds_gsea, coef = "agent_DPN_vs_Control", type = "apeglm",
                     res = res)

rnk <- data.frame(Feature = rownames(res.ape), LFC = res.ape$log2FoldChange)
rnk$Feature <- str_remove(rnk$Feature, "\\..*$")

head(rnk)
##Se guarda el rango con el nombre "DPN_vs_control.rnk"
write.table(rnk, file = "/home/vant/RNA_seq/practicas/transcriptomic-final-exercise/Apartado2/input/DPN_vs_control.rnk", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
