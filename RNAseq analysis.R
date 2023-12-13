
#install and load the Bioconductor Package Manager and  packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    install.packages("tidyverse")
  BiocManager::install("tximport")
  BiocManager::install('ComplexHeatmap')
  BiocManager::install('EnhancedVolcano')
  BiocManager::install("DESeq2")

library("tidyverse")
library("tximport")
library("DESeq2")
library("ComplexHeatmap")
library("EnhancedVolcano")


#set you working directory as the directory which salmon created the read quantification directores
setwd("{YOUR PATH}/transcript_quant")

#define the condition groups associated with each sample
samples <- data.frame(
  run = c('SRR6671757', 'SRR6671758', 'SRR6671759', 'SRR6671775', 'SRR6671776', 'SRR6671777'), 
  condition = c("control","control","control","meropenem","meropenem","meropenem")
  )

#create an object which points each sample to the location of its .sf read quantification file
files <- file.path(
  paste(
    substring(samples$run, 1),"/quant.sf", 
    sep = ""
    )
  )

names(files) <- samples$run

#create an object which contains the reference file to convert transcript annotations to their gene names
tx2gene <- read.csv("tx2gene.csv")

#concatenate the read quantifications into a single file and convert transcript annotations to gene names with tximport
txi <- tximport(files, type="salmon", tx2gene=tx2gene, importer = read.delim)

head(txi$counts)

#run DESeq2 to calculate differential expression values and write results to a .csv file

dds <- DESeqDataSetFromTximport(txi, colData=samples, design= ~ condition)
dds <- DESeq(dds)
res <- results(dds)
de <- write.table(res,"DESeq2_reads.csv", sep = ",", row.names = TRUE)

#check quality of results with PCA and dispersion model
res_transform <- DESeqTransform(dds)
plotPCA(res_transform)
plotDispEsts(dds)

#make volcano plot to visualise significant DEGs (https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html) 
volcano <- EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue',
    pCutoff = 0.05)

plot(volcano)





