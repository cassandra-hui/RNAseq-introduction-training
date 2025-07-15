library(tximportData)

setwd("~/Documents/Training/RNAseq-introduction-training/my_run")


#dir <- system.file("extdata", package = "tximportData")  # Example
dir <- "Salmon_run"
list.files(dir)


samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
samples

#files <- file.path(dir, "033_d_salmon_mouseMT_", samples$sample, "quant.sf")
files <- file.path(dir, paste0("033_d_salmon_mouseMT_", samples$sample), "quant.sf")
names(files) <- samples$sample
all(file.exists(files))


library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")



# library(readr)
# tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))
# head(tx2gene)


library(tximport)
txi <- tximport(files, type = "salmon", 
                tx2gene = tx2gene)
names(txi)

############
# Only mapped to "MT" not all the genes

sf <- read.table(files[4], header = TRUE)
unique(sf$Name)

tx2gene <- data.frame(TXNAME = "MT", GENEID = "MT")
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(txi)
head(txi$counts)

# Transcript level summary
txi.tx <- tximport(files, type = "salmon", txOut = TRUE)
txi.tx

txi.sum <- summarizeToGene(txi.tx, tx2gene)
all.equal(txi$counts, txi.sum$counts)

## Then run DESeq2

library(DESeq2)

sampleTable <- data.frame(condition = factor(rep(c("A", "B"), each = 4)))
rownames(sampleTable) <- colnames(txi$counts)
sampleTable


dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)

dim(dds)

