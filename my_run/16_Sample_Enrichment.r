library(AnnotationHub)
library(AnnotationDbi)
library(clusterProfiler)
library(ReactomePA)

library(org.Mm.eg.db)


res = read.csv( 'Ruhland2016.DESeq2.results.csv'  , row.names=1)
#let's define significance as padj <0.01 & abs(lfc) > 1
res$sig = abs(res$log2FoldChange)>1 & res$padj<0.01

table( res$sig )

genes_universe <- bitr(rownames(res), fromType = "ENSEMBL",
                       toType = c("ENTREZID", "SYMBOL"),
                       OrgDb = "org.Mm.eg.db")

head( genes_universe )

#ENSEMBL ENTREZID  SYMBOL
#2 ENSMUSG00000033845    27395  Mrpl15
#4 ENSMUSG00000025903    18777  Lypla1
#5 ENSMUSG00000033813    21399   Tcea1
#7 ENSMUSG00000002459    58175   Rgs20
#8 ENSMUSG00000033793   108664 Atp6v1h
#9 ENSMUSG00000025907    12421  Rb1cc1

dim(genes_universe)
# 15443     3

length(rownames(res))
# 18012
genes_DE <- bitr(rownames(res)[which( res$sig==T )], fromType = "ENSEMBL",
                 toType = c("ENTREZID", "SYMBOL"),
                 OrgDb = "org.Mm.eg.db")
dim(genes_DE)
# 382   3
# GO "biological process (BP)" enrichment
ego_bp <- enrichGO(gene          = as.character(unique(genes_DE$ENTREZID)),
                   universe      = as.character(unique(genes_universe$ENTREZID)),
                   OrgDb         = org.Mm.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)
# couple of minutes to run

head(ego_bp)
dotplot(ego_bp, showCategory = 20)
# sample plot, but with adjusted p-value as x-axis
#dotplot(ego_bp, x = "p.adjust", showCategory = 20)

# Reactome pathways enrichment
reactome.enrich <- enrichPathway(gene=as.character(unique(genes_DE$ENTREZID)),
                                 organism = "mouse",
                                 pAdjustMethod = "BH",
                                 qvalueCutoff = 0.01,
                                 readable=T,
                                 universe = genes_universe$ENTREZID)
# <1 minute to run


dotplot(reactome.enrich, x = "p.adjust")



