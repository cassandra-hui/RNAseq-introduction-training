## Load count martix

folder  = "my_run/042_d_STAR_map_raw/"

# we skip the 4 first lines, which contains 
# N_unmapped , N_multimapping , N_noFeature , N_ambiguous   

sample_a1_table = read.table(paste0( folder , "sample_a1" , ".ReadsPerGene.out.tab") , 
                             row.names = 1 , skip = 4 )
head( sample_a1_table )


# Loop through for all for "V2"

raw_counts = data.frame( row.names =  row.names(sample_a1_table) )

for( sample in c('a1','a2','a3','a4','b1','b2','b3','b4') ){
  sample_table = read.table(paste0( folder , "sample_" , sample , ".ReadsPerGene.out.tab") , 
                            row.names = 1 , skip = 4 )
  
  raw_counts[sample] = sample_table[ row.names(raw_counts) , "V2" ]
  
}

head( raw_counts )


## Load libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)

## Set up experimental Design

#note: levels let's us define the reference levels
treatment <- factor( c(rep("a",4), rep("b",4)), levels=c("a", "b") )
colData <- data.frame(treatment, row.names = colnames(raw_counts))
colData


## Creating DESeq2 object and some QC

dds <- DESeqDataSetFromMatrix(
  countData = raw_counts, colData = colData, 
  design = ~ treatment)
dim(dds)

# Filter low count
# Keep genes with at least 1 read in at least 4 samples

idx <- rowSums(counts(dds, normalized=FALSE) >= 1) >= 4
dds.f <- dds[idx, ]
dim(dds.f)

# Estimate size factors
dds.f <- DESeq(dds.f)

# Plot PCA

vsd <- varianceStabilizingTransformation(dds.f, blind=TRUE )
pcaData <- plotPCA(vsd, intgroup=c("treatment"))
pcaData + geom_label(aes(x=PC1,y=PC2,label=name))

# Remove outliers

raw_counts_no_outliers = raw_counts[ , !( colnames(raw_counts) %in% c('a4','b3') ) ]

treatment <- factor( c(rep("a",3), rep("b",3)), levels=c("a", "b") )
colData <- data.frame(treatment, row.names = colnames(raw_counts_no_outliers))
colData 


dds <- DESeqDataSetFromMatrix(
  countData = raw_counts_no_outliers, colData = colData, 
  design = ~ treatment)
dim(dds)

# Filter low count with 3 samples now

idx <- rowSums(counts(dds, normalized=FALSE) >= 1) >= 3
dds.f <- dds[idx, ]
dim(dds.f)


dds.f <- DESeq(dds.f)



vsd <- varianceStabilizingTransformation(dds.f, blind=TRUE )
pcaData <- plotPCA(vsd, intgroup=c("treatment"))
pcaData + geom_label(aes(x=PC1,y=PC2,label=name))

plotDispEsts(dds.f)


# extracting results for the treatment versus control contrast
res <- results(dds.f)
summary( res )

head(coef(dds.f)) # the second column corresponds to the difference between the 2 conditions

res.lfc <- lfcShrink(dds.f, coef=2, res=res)
DESeq2::plotMA(res.lfc)

# Volcano Plot

FDRthreshold = 0.01
logFCthreshold = 1.0
# add a column of NAs
res.lfc$diffexpressed <- "NO"
# if log2Foldchange > 1 and pvalue < 0.01, set as "UP" 
res.lfc$diffexpressed[res.lfc$log2FoldChange > logFCthreshold & res.lfc$padj < FDRthreshold] <- "UP"
# if log2Foldchange < 1 and pvalue < 0.01, set as "DOWN"
res.lfc$diffexpressed[res.lfc$log2FoldChange < -logFCthreshold & res.lfc$padj < FDRthreshold] <- "DOWN"

ggplot( data = data.frame( res.lfc ) , aes( x=log2FoldChange , y = -log10(padj) , col =diffexpressed ) ) + 
  geom_point() + 
  geom_vline(xintercept=c(-logFCthreshold, logFCthreshold), col="red") +
  geom_hline(yintercept=-log10(FDRthreshold), col="red") +
  scale_color_manual(values=c("blue", "grey", "red"))

table(res.lfc$diffexpressed)

## Heatmap

vsd.counts <- assay(vsd)

topVarGenes <- head(order(rowVars(vsd.counts), decreasing = TRUE), 20)
mat  <- vsd.counts[ topVarGenes, ] #scaled counts of the top genes
mat  <- mat - rowMeans(mat)  # centering
pheatmap(mat)


# save
write.csv( res ,'051_r_mouseMT.DESeq2.results.csv' )

