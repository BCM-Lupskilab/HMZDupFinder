# HMZDupFinder

## Prerequisites

R version \>= 4.2 Following R libraries are required:

Shifting level models based segmentation is performed using [SLMSuite](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1734-5).

## install packages
``` r
install.packages("devtools")
devtools::install_github("cluhaowie/HMZDupFinder")
devtools::install_github("cluhaowie/VizCNV/SLMSeg")
library(HMZDupFinder)
library(IRanges)
```

## run on public availiable samples
``` r
library(WES.1KG.WUGSC)
library(BSgenome.Hsapiens.UCSC.hg19)
ref.genome <- BSgenome.Hsapiens.UCSC.hg19
dirPath <- system.file("extdata", package = "WES.1KG.WUGSC")
bamFile <- list.files(dirPath, pattern = '*.bam$',full.names = T)
sampname <- as.character(unlist(read.table(file.path(dirPath, "sampname"))))
bedFile <- file.path(dirPath, "chr22_400_to_500.bed")
bedFile <- "~/Downloads/HMZDupFinder_test/chr22_400_to_500.bed"
# prepare bed file
bedOrdered <- prepareBed(bedFile,ref.genome = ref.genome)

# generate TPM data frame
tpmFile <- list.files(outputDir,pattern = "*tpm.bed$",full.names = T)
tpmDtOrdered <- prepareTPMData(tpmFile,mc.cores = 2)

## compare the TPM profile to identify reference st
perMat <- cor(tpmDtOrdered,method = "pearson")

## write out the Z-TPM and log2 ratio
CalcuZtpm(perMat,tpmDtOrdered,bedOrdered,outputDir,mc.cores = 4)
```
