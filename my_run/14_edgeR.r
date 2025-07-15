# setup
library(edgeR)
library(ggplot2)

# Experimental Design

treatment <- factor( c(rep("a",4), rep("b",4)), levels=c("a", "b") )
names(treatment) = colnames(raw_counts)

treatment


# edgeR object and QC

dge.all <- DGEList(counts = raw_counts , group = treatment)  

dge.f.design <- model.matrix(~ treatment)

# filtering by expression level. See ?filterByExpr for details
keep <- filterByExpr(dge.all)
dge.f <- dge.all[keep, keep.lib.sizes=FALSE]
table( keep )

# Only keeping 9, DEseq2 was 12


#normalization
dge.f <- calcNormFactors(dge.f)
dge.f$samples


plotMDS( dge.f , col = c('cornflowerblue','forestgreen')[treatment] )

# Remove outliers

raw_counts_no_outliers = raw_counts[ , !( colnames(raw_counts) %in% c('a4','b3') ) ]
head( raw_counts_no_outliers )

treatment <- factor( c(rep("a",3), rep("b",3)), levels=c("a", "b") )
colData <- data.frame(treatment, row.names = colnames(raw_counts_no_outliers))
colData 

dge.all <- DGEList(counts = raw_counts_no_outliers , group = treatment)  

dge.f.design <- model.matrix(~ treatment)

# filtering by expression level. See ?filterByExpr for details
keep <- filterByExpr(dge.all)
dge.f <- dge.all[keep, keep.lib.sizes=FALSE]
table( keep )

#normalization
dge.f <- calcNormFactors(dge.f)
dge.f$samples


plotMDS( dge.f , col = c('cornflowerblue','forestgreen')[treatment] )

# estimate of the dispersion
dge.f <- estimateDisp(dge.f,dge.f.design , robust = T)
plotBCV(dge.f)


# testing for differential expression. 
# This method is recommended when you only have 2 groups to compare
dge.f.et <- exactTest(dge.f)
topTags(dge.f.et) # printing the genes where the p-value of differential expression if the lowest


summary(decideTests(dge.f.et , p.value = 0.01)) # let's use 0.01 as a threshold

## plot all the logFCs versus average count size. Significantly DE genes are  colored
plotMD(dge.f.et)
# lines at a log2FC of 1/-1, corresponding to a shift in expression of x2 
abline(h=c(-1,1), col="blue")
abline(h=c(0), col="grey")

# Volcano Plot

allGenes = topTags(dge.f.et , n = nrow(dge.f.et$table) )$table

FDRthreshold = 0.01
logFCthreshold = 1.0
# add a column of NAs
allGenes$diffexpressed <- "NO"
# if log2Foldchange > 1 and pvalue < 0.01, set as "UP" 
allGenes$diffexpressed[allGenes$logFC > logFCthreshold & allGenes$FDR < FDRthreshold] <- "UP"
# if log2Foldchange < 1 and pvalue < 0.01, set as "DOWN"
allGenes$diffexpressed[allGenes$logFC < -logFCthreshold & allGenes$FDR < FDRthreshold] <- "DOWN"

ggplot( data = allGenes , aes( x=logFC , y = -log10(FDR) , col =diffexpressed ) ) + 
  geom_point() + 
  geom_vline(xintercept=c(-logFCthreshold, logFCthreshold), col="red") +
  geom_hline(yintercept=-log10(FDRthreshold), col="red") +
  scale_color_manual(values=c("blue", "grey", "red"))


# save
write.csv( allGenes , '052_r_mouseMT.edgeR.results.csv')



# Convert Gene IDs

library(clusterProfiler)
library(org.Mm.eg.db)

genes_universe <- bitr(rownames(allGenes), fromType = "ENSEMBL",
                       toType = c("ENTREZID", "SYMBOL"),
                       OrgDb = "org.Mm.eg.db")
genes_universe

