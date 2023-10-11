---
author: Haowei Du, haoweid\@bcm.edu
date: 2023-04-11
generator: pandoc
title: HMZDupFinder vignette
viewport: width=device-width, initial-scale=1
---

::: {.container-fluid .main-container}
::: {#header}
# HMZDupFinder vignette {#hmzdupfinder-vignette .title .toc-ignore}

#### Haowei Du, <haoweid@bcm.edu> {#haowei-du-haoweidbcm.edu .author}

#### 2023-04-11 {#section .date}
:::

::: {#install-the-package-from-github .section .level1}
# Install the package from github

``` {.r}
# require devtools
install.packages("devtools")
devtools::install_github("cluhaowie/HMZDupFinder")
devtools::install_github("cluhaowie/VizCNV/SLMSeg")
```
:::

::: {#prepare-the-normalized-read-depth .section .level1}
# Prepare the normalized read depth

``` {.r}
library(HMZDupFinder)
library(IRanges)
```

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

``` {.r}
library(WES.1KG.WUGSC)
library(BSgenome.Hsapiens.UCSC.hg19)
```

    ## Loading required package: BSgenome

    ## Loading required package: GenomeInfoDb

    ## Loading required package: GenomicRanges

    ## Loading required package: Biostrings

    ## Loading required package: XVector

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## Loading required package: rtracklayer

``` {.r}
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:Biostrings':
    ## 
    ##     collapse, intersect, setdiff, setequal, union

    ## The following object is masked from 'package:XVector':
    ## 
    ##     slice

    ## The following objects are masked from 'package:GenomicRanges':
    ## 
    ##     intersect, setdiff, union

    ## The following object is masked from 'package:GenomeInfoDb':
    ## 
    ##     intersect

    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` {.r}
library(data.table)
```

    ## 
    ## Attaching package: 'data.table'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     between, first, last

    ## The following object is masked from 'package:GenomicRanges':
    ## 
    ##     shift

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     shift

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, second

``` {.r}
ref.genome <- BSgenome.Hsapiens.UCSC.hg19
dirPath1 <- system.file("extdata", package = "HMZDupFinder")
dirPath2 <- system.file("extdata", package = "WES.1KG.WUGSC")
outputDir <- "./output/"
bamFiles <- list.files(dirPath2, pattern = '*.bam$',full.names = T)
sampname <- as.character(unlist(read.table(file.path(dirPath2, "sampname"))))
bedFile <- file.path(dirPath1, "chr22_400_to_500.bed")
# prepare bed file
bedOrdered <- prepareBed(bedFile,ref.genome = ref.genome)
```

    ## [1] "***Reading bed file***"
    ## [1] "***Annotating GC ratio***"
    ## [1] "Calculating GC-content..."

``` {.r}
## calculate TPM
calcTPMsFromBAMs(bedOrdered,bamFiles,sampname,outputDir,mc.cores = 1)
```

    ## [1] "NA06994"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA06994.mapped.ILLUMINA.bwa.CEU.exome.201205 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA06994.mapped.ILLUMINA.bwa.CEU.exome.20120522.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 29877                                                ||
    ## ||    Successfully assigned alignments : 10881 (36.4%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA10847"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA10847.mapped.ILLUMINA.bwa.CEU.exome.201212 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA10847.mapped.ILLUMINA.bwa.CEU.exome.20121211.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 38109                                                ||
    ## ||    Successfully assigned alignments : 12663 (33.2%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA11840"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA11840.mapped.ILLUMINA.bwa.CEU.exome.201205 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA11840.mapped.ILLUMINA.bwa.CEU.exome.20120522.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 40677                                                ||
    ## ||    Successfully assigned alignments : 12854 (31.6%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA12249"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA12249.mapped.ILLUMINA.bwa.CEU.exome.201205 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA12249.mapped.ILLUMINA.bwa.CEU.exome.20120522.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 34559                                                ||
    ## ||    Successfully assigned alignments : 11677 (33.8%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA12716"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA12716.mapped.ILLUMINA.bwa.CEU.exome.201212 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA12716.mapped.ILLUMINA.bwa.CEU.exome.20121211.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 40153                                                ||
    ## ||    Successfully assigned alignments : 13585 (33.8%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA12750"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA12750.mapped.ILLUMINA.bwa.CEU.exome.201304 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA12750.mapped.ILLUMINA.bwa.CEU.exome.20130415.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 140758                                               ||
    ## ||    Successfully assigned alignments : 27965 (19.9%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA12751"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA12751.mapped.ILLUMINA.bwa.CEU.exome.201212 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA12751.mapped.ILLUMINA.bwa.CEU.exome.20121211.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 39854                                                ||
    ## ||    Successfully assigned alignments : 12370 (31.0%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA12760"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA12760.mapped.ILLUMINA.bwa.CEU.exome.201205 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA12760.mapped.ILLUMINA.bwa.CEU.exome.20120522.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 41954                                                ||
    ## ||    Successfully assigned alignments : 13841 (33.0%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA12761"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA12761.mapped.ILLUMINA.bwa.CEU.exome.201212 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA12761.mapped.ILLUMINA.bwa.CEU.exome.20121211.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 35030                                                ||
    ## ||    Successfully assigned alignments : 10877 (31.1%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA12763"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA12763.mapped.ILLUMINA.bwa.CEU.exome.201212 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA12763.mapped.ILLUMINA.bwa.CEU.exome.20121211.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 35536                                                ||
    ## ||    Successfully assigned alignments : 11423 (32.1%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA18966"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA18966.mapped.ILLUMINA.bwa.JPT.exome.201212 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA18966.mapped.ILLUMINA.bwa.JPT.exome.20121211.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 36024                                                ||
    ## ||    Successfully assigned alignments : 11903 (33.0%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA18967"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA18967.mapped.ILLUMINA.bwa.JPT.exome.201212 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA18967.mapped.ILLUMINA.bwa.JPT.exome.20121211.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 42050                                                ||
    ## ||    Successfully assigned alignments : 12559 (29.9%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA18968"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA18968.mapped.ILLUMINA.bwa.JPT.exome.201205 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA18968.mapped.ILLUMINA.bwa.JPT.exome.20120522.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 36688                                                ||
    ## ||    Successfully assigned alignments : 12649 (34.5%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA18969"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA18969.mapped.ILLUMINA.bwa.JPT.exome.201205 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA18969.mapped.ILLUMINA.bwa.JPT.exome.20120522.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 36995                                                ||
    ## ||    Successfully assigned alignments : 11725 (31.7%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA18970"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA18970.mapped.ILLUMINA.bwa.JPT.exome.201205 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA18970.mapped.ILLUMINA.bwa.JPT.exome.20120522.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 39733                                                ||
    ## ||    Successfully assigned alignments : 13045 (32.8%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA18971"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA18971.mapped.ILLUMINA.bwa.JPT.exome.201205 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA18971.mapped.ILLUMINA.bwa.JPT.exome.20120522.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 38453                                                ||
    ## ||    Successfully assigned alignments : 12492 (32.5%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA18972"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA18972.mapped.ILLUMINA.bwa.JPT.exome.201212 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA18972.mapped.ILLUMINA.bwa.JPT.exome.20121211.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 26012                                                ||
    ## ||    Successfully assigned alignments : 10115 (38.9%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA18973"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA18973.mapped.ILLUMINA.bwa.JPT.exome.201212 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA18973.mapped.ILLUMINA.bwa.JPT.exome.20121211.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 35690                                                ||
    ## ||    Successfully assigned alignments : 11757 (32.9%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA18974"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA18974.mapped.ILLUMINA.bwa.JPT.exome.201205 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA18974.mapped.ILLUMINA.bwa.JPT.exome.20120522.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 36298                                                ||
    ## ||    Successfully assigned alignments : 11824 (32.6%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA18975"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA18975.mapped.ILLUMINA.bwa.JPT.exome.201205 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA18975.mapped.ILLUMINA.bwa.JPT.exome.20120522.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 37632                                                ||
    ## ||    Successfully assigned alignments : 11870 (31.5%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA18976"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA18976.mapped.ILLUMINA.bwa.JPT.exome.201205 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA18976.mapped.ILLUMINA.bwa.JPT.exome.20120522.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 36486                                                ||
    ## ||    Successfully assigned alignments : 11832 (32.4%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA18981"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA18981.mapped.ILLUMINA.bwa.JPT.exome.201205 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA18981.mapped.ILLUMINA.bwa.JPT.exome.20120522.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 33428                                                ||
    ## ||    Successfully assigned alignments : 11218 (33.6%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA18987"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA18987.mapped.ILLUMINA.bwa.JPT.exome.201205 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA18987.mapped.ILLUMINA.bwa.JPT.exome.20120522.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 36112                                                ||
    ## ||    Successfully assigned alignments : 11691 (32.4%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA18990"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA18990.mapped.ILLUMINA.bwa.JPT.exome.201205 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA18990.mapped.ILLUMINA.bwa.JPT.exome.20120522.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 38362                                                ||
    ## ||    Successfully assigned alignments : 12362 (32.2%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA18991"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA18991.mapped.ILLUMINA.bwa.JPT.exome.201205 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA18991.mapped.ILLUMINA.bwa.JPT.exome.20120522.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 31927                                                ||
    ## ||    Successfully assigned alignments : 11084 (34.7%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA19098"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA19098.mapped.ILLUMINA.bwa.YRI.exome.201212 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA19098.mapped.ILLUMINA.bwa.YRI.exome.20121211.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 40178                                                ||
    ## ||    Successfully assigned alignments : 12481 (31.1%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA19119"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA19119.mapped.ILLUMINA.bwa.YRI.exome.201205 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA19119.mapped.ILLUMINA.bwa.YRI.exome.20120522.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 39801                                                ||
    ## ||    Successfully assigned alignments : 13462 (33.8%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA19131"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA19131.mapped.ILLUMINA.bwa.YRI.exome.201205 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA19131.mapped.ILLUMINA.bwa.YRI.exome.20120522.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 42989                                                ||
    ## ||    Successfully assigned alignments : 13127 (30.5%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA19137"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA19137.mapped.ILLUMINA.bwa.YRI.exome.201304 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA19137.mapped.ILLUMINA.bwa.YRI.exome.20130415.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 173640                                               ||
    ## ||    Successfully assigned alignments : 34871 (20.1%)                        ||
    ## ||    Running time : 0.01 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA19138"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA19138.mapped.ILLUMINA.bwa.YRI.exome.201212 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA19138.mapped.ILLUMINA.bwa.YRI.exome.20121211.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 39870                                                ||
    ## ||    Successfully assigned alignments : 12906 (32.4%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA19141"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA19141.mapped.ILLUMINA.bwa.YRI.exome.201205 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA19141.mapped.ILLUMINA.bwa.YRI.exome.20120522.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 36872                                                ||
    ## ||    Successfully assigned alignments : 12175 (33.0%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA19143"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA19143.mapped.ILLUMINA.bwa.YRI.exome.201205 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA19143.mapped.ILLUMINA.bwa.YRI.exome.20120522.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 38251                                                ||
    ## ||    Successfully assigned alignments : 12752 (33.3%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA19144"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA19144.mapped.ILLUMINA.bwa.YRI.exome.201212 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA19144.mapped.ILLUMINA.bwa.YRI.exome.20121211.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 42838                                                ||
    ## ||    Successfully assigned alignments : 14056 (32.8%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA19152"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA19152.mapped.ILLUMINA.bwa.YRI.exome.201205 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA19152.mapped.ILLUMINA.bwa.YRI.exome.20120522.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 39319                                                ||
    ## ||    Successfully assigned alignments : 13722 (34.9%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA19153"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA19153.mapped.ILLUMINA.bwa.YRI.exome.201212 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA19153.mapped.ILLUMINA.bwa.YRI.exome.20121211.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 43819                                                ||
    ## ||    Successfully assigned alignments : 13374 (30.5%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA19159"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA19159.mapped.ILLUMINA.bwa.YRI.exome.201212 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA19159.mapped.ILLUMINA.bwa.YRI.exome.20121211.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 39538                                                ||
    ## ||    Successfully assigned alignments : 12605 (31.9%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA19160"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA19160.mapped.ILLUMINA.bwa.YRI.exome.201205 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA19160.mapped.ILLUMINA.bwa.YRI.exome.20120522.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 38858                                                ||
    ## ||    Successfully assigned alignments : 13302 (34.2%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA19171"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA19171.mapped.ILLUMINA.bwa.YRI.exome.201205 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA19171.mapped.ILLUMINA.bwa.YRI.exome.20120522.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 48717                                                ||
    ## ||    Successfully assigned alignments : 15420 (31.7%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA19200"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA19200.mapped.ILLUMINA.bwa.YRI.exome.201205 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA19200.mapped.ILLUMINA.bwa.YRI.exome.20120522.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 42168                                                ||
    ## ||    Successfully assigned alignments : 14302 (33.9%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA19201"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA19201.mapped.ILLUMINA.bwa.YRI.exome.201212 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA19201.mapped.ILLUMINA.bwa.YRI.exome.20121211.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 44242                                                ||
    ## ||    Successfully assigned alignments : 14673 (33.2%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA19204"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA19204.mapped.ILLUMINA.bwa.YRI.exome.201205 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA19204.mapped.ILLUMINA.bwa.YRI.exome.20120522.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 33022                                                ||
    ## ||    Successfully assigned alignments : 11562 (35.0%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA19206"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA19206.mapped.ILLUMINA.bwa.YRI.exome.201212 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA19206.mapped.ILLUMINA.bwa.YRI.exome.20121211.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 43325                                                ||
    ## ||    Successfully assigned alignments : 13962 (32.2%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA19207"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA19207.mapped.ILLUMINA.bwa.YRI.exome.201212 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA19207.mapped.ILLUMINA.bwa.YRI.exome.20121211.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 40174                                                ||
    ## ||    Successfully assigned alignments : 13357 (33.2%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA19209"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA19209.mapped.ILLUMINA.bwa.YRI.exome.201212 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA19209.mapped.ILLUMINA.bwa.YRI.exome.20121211.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 49423                                                ||
    ## ||    Successfully assigned alignments : 15362 (31.1%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA19210"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA19210.mapped.ILLUMINA.bwa.YRI.exome.201212 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA19210.mapped.ILLUMINA.bwa.YRI.exome.20121211.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 35854                                                ||
    ## ||    Successfully assigned alignments : 13674 (38.1%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## [1] "NA19223"
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.12.3
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           NA19223.mapped.ILLUMINA.bwa.YRI.exome.201205 ... ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : R data.frame                                     ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 3                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : counted                                          ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file .Rsubread_UserProvidedAnnotation_pid17435 ...         ||
    ## ||    Features : 100                                                          ||
    ## ||    Meta-features : 100                                                     ||
    ## ||    Chromosomes/contigs : 1                                                 ||
    ## ||                                                                            ||
    ## || Process BAM file NA19223.mapped.ILLUMINA.bwa.YRI.exome.20120522.bam.ch ... ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 40866                                                ||
    ## ||    Successfully assigned alignments : 13806 (33.8%)                        ||
    ## ||    Running time : 0.00 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//

``` {.r}
# generate TPM data frame
tpmFile <- list.files(outputDir,pattern = "*tpm.bed$",full.names = T)
tpmDtOrdered <- prepareTPMData(tpmFile,mc.cores = 2)
```

    ## [1] "Reading TPM files ..."
    ## [1] "Removing empty elements ..."
    ## [1] "Creating matrix ..."
    ## [1] "Creating data.table (may take a while)..."

``` {.r}
## compare the TPM profile to identify reference st
perMat <- cor(tpmDtOrdered,method = "pearson")

## write out the Z-TPM and log2 ratio
CalcuZtpm(perMat,tpmDtOrdered,bedOrdered,outputDir,mc.cores = 4)
```

    ## [1] "[******Calculate Z-TPM**********]"
    ## [1] "[******Calculate Ratio**********]"
    ## [1] "[******Binding the values**********]"
    ## [1] "[******Write z-tpm and ratio to individual file**********]"
    ## [1] "Writring ./output//NA06994.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA10847.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA11840.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA12249.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA12716.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA12750.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA12751.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA12760.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA12761.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA12763.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA18966.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA18967.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA18968.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA18969.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA18970.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA18971.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA18972.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA18973.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA18974.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA18975.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA18976.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA18981.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA18987.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA18990.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA18991.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA19098.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA19119.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA19131.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA19137.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA19138.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA19141.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA19143.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA19144.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA19152.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA19153.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA19159.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA19160.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA19171.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA19200.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA19201.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA19204.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA19206.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA19207.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA19209.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA19210.ztpm.ratio.bed"
    ## [1] "Writring ./output//NA19223.ztpm.ratio.bed"

    ## [[1]]
    ## NULL
    ## 
    ## [[2]]
    ## NULL
    ## 
    ## [[3]]
    ## NULL
    ## 
    ## [[4]]
    ## NULL
    ## 
    ## [[5]]
    ## NULL
    ## 
    ## [[6]]
    ## NULL
    ## 
    ## [[7]]
    ## NULL
    ## 
    ## [[8]]
    ## NULL
    ## 
    ## [[9]]
    ## NULL
    ## 
    ## [[10]]
    ## NULL
    ## 
    ## [[11]]
    ## NULL
    ## 
    ## [[12]]
    ## NULL
    ## 
    ## [[13]]
    ## NULL
    ## 
    ## [[14]]
    ## NULL
    ## 
    ## [[15]]
    ## NULL
    ## 
    ## [[16]]
    ## NULL
    ## 
    ## [[17]]
    ## NULL
    ## 
    ## [[18]]
    ## NULL
    ## 
    ## [[19]]
    ## NULL
    ## 
    ## [[20]]
    ## NULL
    ## 
    ## [[21]]
    ## NULL
    ## 
    ## [[22]]
    ## NULL
    ## 
    ## [[23]]
    ## NULL
    ## 
    ## [[24]]
    ## NULL
    ## 
    ## [[25]]
    ## NULL
    ## 
    ## [[26]]
    ## NULL
    ## 
    ## [[27]]
    ## NULL
    ## 
    ## [[28]]
    ## NULL
    ## 
    ## [[29]]
    ## NULL
    ## 
    ## [[30]]
    ## NULL
    ## 
    ## [[31]]
    ## NULL
    ## 
    ## [[32]]
    ## NULL
    ## 
    ## [[33]]
    ## NULL
    ## 
    ## [[34]]
    ## NULL
    ## 
    ## [[35]]
    ## NULL
    ## 
    ## [[36]]
    ## NULL
    ## 
    ## [[37]]
    ## NULL
    ## 
    ## [[38]]
    ## NULL
    ## 
    ## [[39]]
    ## NULL
    ## 
    ## [[40]]
    ## NULL
    ## 
    ## [[41]]
    ## NULL
    ## 
    ## [[42]]
    ## NULL
    ## 
    ## [[43]]
    ## NULL
    ## 
    ## [[44]]
    ## NULL
    ## 
    ## [[45]]
    ## NULL
    ## 
    ## [[46]]
    ## NULL
:::

::: {#call-dups-with-custom-defined-parameters .section .level1}
# Call dups with custom defined parameters

``` {.r}
## call homozygous duplication with fine-tune cut-off
hmzdup.seg.mean <- c(0.85,1.15)
ztpm.mean.cutoff <- 1.5 ## with the optimized ztpm cutoff(>4.0), there are no call in the toy data set, here set 1.5 to get the low confidence call
SegNormRDES <- function(df,id,seg.method="slm"){
  ##EDIT: include SLM segmentation
  if (seg.method == "slm") {
    print("segment with SLM")
    df.ls <- base::split(df, df$seqnames)
    res <- lapply(df.ls, function(df) {
      logratio <- log2(df$ratio + 0.001)
      logratio[is.na(logratio)] <- 0
      slm <-
        SLMSeg::HSLM(
          logratio,
          pos_data = (df$start+df$end)/2,
          omega = 0.7,
          FW = 0,
          eta = 1e-5,
          stepeta=1000
        )
      res <- rle(slm[1, ])
      idx <- sapply(seq_along(res$lengths),function(i){
        if(i==1){return(1)}
        start.idx=1+sum(res$lengths[1:(i-1)])
        return(start.idx)
      })
      chr=df$seqnames[idx]
      start=df$start[idx]
      end=c(df$start[c(idx[-1],end(df$start)[1])])
      mean_ztpm <- sapply(seq_along(res$lengths),
                          function(i){
                            if(i==1){
                              return(mean(df$ztpm[1:res$lengths[i]]))
                            }
                            z <- df$ztpm[(1+sum(res$lengths[1:(i-1)])):(sum(res$lengths[1:i]))]
                            return(mean(z))
                          })
      res.dt <- data.table(ID=id,chrom=chr,loc.start=start,loc.end=end,
                           num.mark=res$lengths,seg.mean=res$values,
                           ztpm.mean=mean_ztpm)
    })
    res <- data.table::rbindlist(res)
  }
}

## get all hmz dup call 
ztpmFile <- list.files(outputDir,pattern = "*ztpm.ratio.bed",full.names = T)
tmplist <- lapply(sampname,function(i){
  file <- paste0(outputDir,i,".ztpm.ratio.bed")
  print(i)
  df <- fread(file,stringsAsFactors = F)
  seg <- tryCatch({
    SegNormRDES(df,i,seg.method="slm")
  },error=function(e){
    return(i)
  })
  if(!is.data.frame(seg)){
    return(seg)
  }else{
    hmzdup <- seg%>%
      mutate(cnv=ifelse(dplyr::between(seg.mean,hmzdup.seg.mean[1], hmzdup.seg.mean[2])&ztpm.mean>ztpm.mean.cutoff,"hmzdup","else"))%>%
      filter(cnv=="hmzdup")
    return(hmzdup)
  }
  
})
```

    ## [1] "NA06994"
    ## [1] "segment with SLM"
    ## [1] "NA10847"
    ## [1] "segment with SLM"
    ## [1] "NA11840"
    ## [1] "segment with SLM"
    ## [1] "NA12249"
    ## [1] "segment with SLM"
    ## [1] "NA12716"
    ## [1] "segment with SLM"
    ## [1] "NA12750"
    ## [1] "segment with SLM"
    ## [1] "NA12751"
    ## [1] "segment with SLM"
    ## [1] "NA12760"
    ## [1] "segment with SLM"
    ## [1] "NA12761"
    ## [1] "segment with SLM"
    ## [1] "NA12763"
    ## [1] "segment with SLM"
    ## [1] "NA18966"
    ## [1] "segment with SLM"
    ## [1] "NA18967"
    ## [1] "segment with SLM"
    ## [1] "NA18968"
    ## [1] "segment with SLM"
    ## [1] "NA18969"
    ## [1] "segment with SLM"
    ## [1] "NA18970"
    ## [1] "segment with SLM"
    ## [1] "NA18971"
    ## [1] "segment with SLM"
    ## [1] "NA18972"
    ## [1] "segment with SLM"
    ## [1] "NA18973"
    ## [1] "segment with SLM"
    ## [1] "NA18974"
    ## [1] "segment with SLM"
    ## [1] "NA18975"
    ## [1] "segment with SLM"
    ## [1] "NA18976"
    ## [1] "segment with SLM"
    ## [1] "NA18981"
    ## [1] "segment with SLM"
    ## [1] "NA18987"
    ## [1] "segment with SLM"
    ## [1] "NA18990"
    ## [1] "segment with SLM"
    ## [1] "NA18991"
    ## [1] "segment with SLM"
    ## [1] "NA19098"
    ## [1] "segment with SLM"
    ## [1] "NA19119"
    ## [1] "segment with SLM"
    ## [1] "NA19131"
    ## [1] "segment with SLM"
    ## [1] "NA19137"
    ## [1] "segment with SLM"
    ## [1] "NA19138"
    ## [1] "segment with SLM"
    ## [1] "NA19141"
    ## [1] "segment with SLM"
    ## [1] "NA19143"
    ## [1] "segment with SLM"
    ## [1] "NA19144"
    ## [1] "segment with SLM"
    ## [1] "NA19152"
    ## [1] "segment with SLM"
    ## [1] "NA19153"
    ## [1] "segment with SLM"
    ## [1] "NA19159"
    ## [1] "segment with SLM"
    ## [1] "NA19160"
    ## [1] "segment with SLM"
    ## [1] "NA19171"
    ## [1] "segment with SLM"
    ## [1] "NA19200"
    ## [1] "segment with SLM"
    ## [1] "NA19201"
    ## [1] "segment with SLM"
    ## [1] "NA19204"
    ## [1] "segment with SLM"
    ## [1] "NA19206"
    ## [1] "segment with SLM"
    ## [1] "NA19207"
    ## [1] "segment with SLM"
    ## [1] "NA19209"
    ## [1] "segment with SLM"
    ## [1] "NA19210"
    ## [1] "segment with SLM"
    ## [1] "NA19223"
    ## [1] "segment with SLM"
:::

::: {#get-calls-from-all-samples .section .level1}
# Get calls from all samples

``` {.r}
dupcall <- rbindlist(tmplist)
dupcall
```

    ##         ID chrom loc.start  loc.end num.mark  seg.mean ztpm.mean    cnv
    ## 1: NA10847 chr22  21797070 21799161        1 1.0369728  2.383095 hmzdup
    ## 2: NA12249 chr22  21921983 21947100        1 0.9437214  1.598907 hmzdup
    ## 3: NA18969 chr22  21356071 21369444        1 1.0549529  1.890923 hmzdup
    ## 4: NA18971 chr22  21987420 21988261        1 1.0836136  1.781615 hmzdup
    ## 5: NA18974 chr22  21987420 21988261        1 0.9754562  1.537753 hmzdup
    ## 6: NA19119 chr22  21356071 21369444        1 0.9044307  1.529398 hmzdup
    ## 7: NA19144 chr22  21797070 21799161        1 0.9033155  1.971388 hmzdup
:::

::: {#visualize-the-normalized-read-depth .section .level1}
# Visualize the normalized read depth

``` {.r}
library(ggplot2)
library(dplyr)
library(stringr)
scale_rd <- scale_y_continuous(name="Log2 Ratio",
                                 limits=c(-2.5, 2),
                                 breaks = c(-2,
                                            round(log2(1/2),2),
                                            round(log2(2/2),2),
                                            round(log2(3/2),2),
                                            round(log2(4/2),2),
                                            round(log2(5/2),2),
                                            round(log2(6/2),2)
                                 ))
wg_ratio2plot <- function(df,id){
  temp <- df %>% 
    group_by(seqnames) %>% 
    summarise(max_end = max(end)) %>% 
    mutate(across("seqnames", str_replace, "chr", "")) %>% 
    arrange(as.numeric(seqnames)) %>% 
    mutate(loc_add = lag(cumsum(as.numeric(max_end)), default = 0)) %>% 
    mutate(seqnames = paste0("chr", seqnames))%>%
    mutate(loc_end=loc_add+max_end)

  para.label <- paste0(id)

  df <- df%>%mutate(log2ratio=log2(ratio+0.001))%>%
    mutate(log2ratio=ifelse(log2ratio < -2.5,-2,log2ratio))
  df <- df %>% 
    inner_join(temp, by = "seqnames") %>% 
    mutate(end_cum = loc_add + end) 
  df <- df %>% 
    mutate(start_cum = loc_add+start)
  axis_set <- df %>% 
    group_by(seqnames) %>% 
    summarize(center = mean(end_cum)) %>% 
    arrange((seqnames))

  label_seg_gain <- df %>% 
    filter(dplyr::between(log2ratio,log2(3/2)-0.2,log2(3/2)+0.5))
  label_seg_loss <- df %>% 
    filter(dplyr::between(log2ratio,log2(1/2)-0.5,log2(1/2)+0.5))
  
  wg <-  
    ggplot()+
    geom_point(data=df,aes(x = start_cum, y = log2ratio))+
    geom_point(data = label_seg_gain,aes(x = start_cum, y = log2ratio), color = "red")+
    geom_point(data = label_seg_loss,aes(x = start_cum, y = log2ratio), color = "green")+
    #geom_rect(data = temp[seq(1,24,2),],aes(xmin=loc_add,xmax=loc_end,ymin=-2,ymax=1.58),alpha = 0.3)+
    #geom_point(data = label_seg_gain, shape = 8, color = "red")+
    #geom_point(data = label_seg_gain.sig, shape = 8, color = "blue")+
    #geom_point(data = label_seg_loss, shape= 8, color = "green")+
    theme_minimal() +
    theme( 
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 60, size = 10, vjust = 0.8)
    )+
    scale_rd+
    #scale_size_continuous(range = c(0.5,3))+
    labs(x = NULL)+
    scale_colour_manual(values = c("orange", "purple"))+
    scale_x_continuous(label = axis_set$seqnames, breaks = axis_set$center)+
    coord_cartesian(expand = F)+
    ggtitle(para.label)
  return(wg)
}
df <- fread("output/NA10847.ztpm.ratio.bed",stringsAsFactors = F)
wg_ratio2plot(df,"NA10847")
```

![](data:image/png;base64,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ECBAgAABAgQIECBAgAABAvUFzACtb9Lqc8aNGxdbb711/ueUU04p+plOPfXUfP0LL7yw6PrHHXdcvn4ay29+85ui22iqwh//+MfYa6+9Yu+994733nuvqaL17qXyxx57bHTr1i0Leo4YMSK++MUvZkHFnXbaKR5++OF6dRrLmDp1avyf//N/suDuOuusE8OHD48hQ4ZkbW+++eZx5513Nla1qPzTTjst7zl69OgG6959992x/vrrl/yz0UYbNdiuTAIECBAgQIAAAQIECBAgQIBAaxYwA7Q1v70Cx37RRRdlQbodd9yxwBrLVyzNQL355pvrNHLJJZfET37yk+jUafl/5Z599tk4+eST67Rf6Jennnoqdtttt/joo4/qVUmzSh955JHYY4894tprr4399tuvXpnaGb/73e/ixBNPjM8++6x2dna9YMGCSONMAdp99903rrrqqmyWab2CBWSMHTu2oADyp59+GmnmbqlphRVWKLWqegQIECBAgAABAgQIECBAgACBqhUwA7RqX03zDSwF9g477LBIAbKWSCl4OHfu3KyrtBw8penTp8ftt9+eXS/PHy+88EK25P3jjz8uupm0PP2rX/1qPviZZkumma433nhjnHnmmTFs2LCszXnz5sXIkSPjiiuuaLSPNEv0mGOOyQc/U1vnnHNONuPzz3/+c+y///75urfeemukGbGlpLQNwSGHHBLpHZY7rb766uXuQvsECBAgQIAAAQIECBAgQIAAgRYXEABtcfLKdPjqq69mwb6W6D0FAFNKsz3POuusfJe///3v89elXKRZlF/+8peLXvKe6+sHP/hBvm5a8p5mg/73f/93Njs2LTF//PHH4/vf/35WfOHChZHKp9msS6c04/Oggw7K9hBN99JS/MmTJ8dJJ52UXadg80033ZTN+szVveyyy0oKAKfxvPnmm7lmmvxMs03feuutgn/+8Ic/5NtL+5X+5S9/yX93QYAAAQIECBAgQIAAAQIECBBoKwICoG3lTRbwHL/97W+L2t+ygCbrFUnLvp9++uksP82o/OY3vxldunTJvj/44IPx8ssv16uzrIx0uNBXvvKVOPjgg+ODDz5YVvEG748fPz7uv//+7F46jOmuu+6K1VZbrU7ZmpqaSMvaDzjggCx/zpw5cfnll9cpk75MmDAh0phSWnnllbMyXbt2zb7X/uM73/lOFijN5V1//fW5y4I+r7766sjVKWTrgHTQVTqsqpCfNMv1v/7rv/LjSNskfOlLX8p/d0GAAAECBAgQIECAAAECBAgQaCsCAqBt5U028hxpWfOAAQOyuy2xFD7NdMylPffcM1Jg8Gtf+1q+/3R4UaEpjTfNUkyHCqUT33Pp61//enTs2DH3taDPiy++OF/uyCOPrBf8zN9ccnH66afnv6aA6NLLz5955pn8/RTk7d27d/770he5Z0/5testXW7p7+lwpaOPPjrLHjx4cD4ou3S5Ur7nlvin5fUpffe7343//M//LKUpdQgQIECAAAECBAgQIECAAAECVS8gAFr1r2j5BphmX6Yl6Wl2Y0qvvfZalHIqfCGj+PzzzyPNWsylFABNKc2EzKU0ozK3P2gur7HPtN9nWgKe2+8zBVPTYUppVmSHDsX96j766KP5btJs0qZSCriuu+66WZE0Y/Whhx6qU/zdd9/Nf08nyTeVau+rmZs12lT5dG/RokXZbNd0UFM6mCiZptmdzZXOOOOM/CzdtEfrueee21xNa4cAAQIECBAgQIAAAQIECBAgUHUCxUWRqm74BlSIwPDhw/N7W6byaVbjuHHjCqlaVJnbbrstZs6cmdVZe+21Y6uttsqu06nquUBgul/KXpNpn820tP6II44oakyp8HvvvZcFftN1CgQXstS7dpknn3wyVc2ndOBRLk2cOLHeDNHcvfSZ7ufSeuutl7ts8vN//ud/Ihew/fnPf553bLJSgTf/9a9/1TlRPgU/e/XqVWBtxQgQIECAAAECBAgQIECAAAECrU9AALT1vbOSRpyCarkAXFrSffjhh8cnn3xSUluNVaq9/D3t15mbdZr2rzzwwAPz1WofvpPPbOAi1U9LyNP+nXfeeWdsvPHGDZRadlZuT9JUMu2PmduTtKmaOatU5vnnn69TNB02lJv5mWbUNvY8aYn52Wefna/7jW98I3/d2EUaawp6prTddts1+8FVxx57bKSZuinttttu8e1vfzu79gcBAgQIECBAgAABAgQIECBAoK0KCIC21Te71HN17969zlL4KVOmxMknn7xUqdK/ppPK77vvvnwDhxxySP46XRx66KH570888URMmjQp/72xi7T35a233hpbb711Y0UKyk8zQHMpt7Q9972xz/79++dvpf04a6d0eFI6PT6XfvjDH2ZBy9yJ8ekE+Ycffjg7sX769OlZsU033TRS8LGplA5dSoHi+fPnZwHWdOp9sXudNtV+mqH717/+NSuS2k0HH0kECBAgQIAAAQIECBAgQIAAgbYuIADa1t9wrefbZZdd4gc/+EE+J81cTCezN0e64oorsr0rU1tp5uKGG25Yp9nNNtuszlLu3//+93Xul/NL2kszl9Kp7YWknj175os1NFP2mGOOiUsvvTQ7ACnNqPzFL34Rq6yySqy11lqRgs0777xzfuboV7/61Rg7dmzUbjPfeK2LE088MV566aUsJy1N32CDDWrdXf7LtPdnLqXZqBtttFHuq08CBAgQIECAAAECBAgQIECAQJsVEABts6+24QdLS+Fze1jmlsLnDhlquMayc1M7o0ePzhesPdszn7nkonb+mDFj8ocb1S5TjuvZs2fnmy1k+XsqXPvQoTQzs6F00EEHxZlnnpkdVJS7P2PGjDqHPKWgczr5vk+fPrkiDX7eddddkQsKp4Bpc5/Kfv/998eECROyvtPWAqeddlqD45BJgAABAgQIECBAgAABAgQIEGhrAgKgbe2NLuN50t6Vaa/O3P6cr7/+epx00knLqNX07bTc+9VXX80KpQBjY3tdpv0mV1xxxaxcmlVZ+8T4pntYvru1Z3B27ty5oMZql/vss8/q1XnqqadiwIABcdRRR2VL1lOBtGx+9913j7TcPZ3enlI6bCoFnEeNGpV9b+iPdKp82pM1pd69e2czSxsqtzx555xzTr76/vvvH5tsskn+uwsCBAgQIECAAAECBAgQIECAQFsWEABty2+3kWfbaaed4v/+3/+bv5tmKD7wwAP578Ve1D78aN99943GlpmnvTPTAUK5lJvxmPters/asznT/pqFpNxBQals7sCjXL2///3vMXz48Hj77bezrLRH6bPPPhvTpk3L9ticPHlypAOQ0t6gHTp0iHnz5sVPfvKTRmddppPt33nnnaytP/3pT7HGGmvkumqWzzTzs/b7NfuzWVg1QoAAAQIECBAgQIAAAQIECLQSgU6tZJyG2cwC6RCfu+++Oz9zM81AfO6556JHjx5F9ZSWl9944435Ommp93XXXZf/vvRF7eBeChT+7W9/i+23337pYs36Pe3JmUtz587NXTb5Wbvc0nt3Hn/88fHpp59m9YcNGxaPPvpopJPua6cUND3vvPNi6NCh+dmdyfyb3/xm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CBAgQIAAAQIECBAgQIAAgfILCICW31gPBAgQIECAAAECBAgQIECAAAECBAhUSEAAtELwuiVAgAABAgQIECBAgAABAgQIECBAoPwCAqDlN9YDAQIECBAgQIAAAQIECBAgQIAAAQIVEhAArRC8bgkQIECAAAECBAgQIECAAAECBAgQKL+AAGj5jfVAgAABAgQIECBAgAABAgQIECBAgECFBARAKwSvWwIECBAgQIAAAQIECBAgQIAAAQIEyi8gAFp+Yz0QIECAAAECBAgQIECAAAECBAgQIFAhAQHQCsHrlgABAgQIECBAgAABAgQIECBAgACB8gsIgJbfWA8ECBAgQIAAAQIECBAgQIAAAQIECFRIQAC0QvC6JUCAAAECBAgQIECAAAECBAgQIECg/AICoOU31gMBAgQIECBAgAABAgQIECBAgAABAhUSEACtELxuCRAgQIAAAQIECBAgQIAAAQIECBAov4AAaPmN9UCAAAECBAgQIECAAAECBAgQIECAQIUEBEArBK9bAgQIECBAgAABAgQIECBAgAABAgTKLyAAWn5jPRAgQIAAAQIECBAgQIAAAQIECBAgUCEBAdAKweuWAAECBAgQIECAAAECBAgQIECAAIHyCwiAlt9YDwQIECBAgAABAgQIECBAgAABAgQIVEhAALRC8LolQIAAAQIECBAgQIAAAQIECBAgQKD8AgKg5TfWAwECBAgQIECAAAECBAgQIECAAAECFRIQAK0QvG4JECBAgAABAgQIECBAgAABAgQIECi/gABo+Y31QIAAAQIECBAgQIAAAQIECBAgQIBAhQQEQCsEr1sCBAgQIECAAAECBAgQIECAAAECBMovIABafmM9ECBAgAABAgQIECBAgAABAgQIECBQIQEB0ArB65YAAQIECBAgQIAAAQIECBAgQIAAgfILCICW31gPBAgQIECAAAECBAgQIECAAAECBAhUSEAAtELwuiVAgAABAgQIECBAgAABAgQIECBAoPwCAqDlN9YDAQIECBAgQIAAAQIECBAgQIAAAQIVEhAArRC8bgkQIECAAAECBAgQIECAAAECBAgQKL+AAGj5jfVAgAABAgQIECBAgAABAgQIECBAgECFBARAKwSvWwIECBAgQIAAAQIECBAgQIAAAQIEyi8gAFp+Yz0QIECAAAECBAgQIECAAAECBAgQIFAhAQHQCsHrlgABAgQIECBAgAABAgQIECBAgACB8gsIgJbfWA8ECBAgQIAAAQIECBAgQIAAAQIECFRIQAC0QvC6JUCAAAECBAgQIECAAAECBAgQIECg/AICoOU31gMBAgQIECBAgAABAgQIECBAgAABAhUSEACtELxuCRAgQIAAAQIECBAgQIAAAQIECBAov4AAaPmN9UCAAAECBAgQIECAAAECBAgQIECAQIUEBEArBK9bAgQIECBAgAABAgQIECBAgAABAgTKLyAAWn5jPRAgQIAAAQIECBAgQIAAAQIECBAgUCEBAdAKweuWAAECBAgQIECAAAECBAgQIECAAIHyCwiAlt9YDwQIECBAgAABAgQIECBAgAABAgQIVEhAALRC8LolQIAAAQIECBAgQIAAAQIECBAgQKD8AgKg5TfWAwECBAgQIECAAAECBAgQIECAAAECFRIQAK0QvG4JECBAgAABAgQIECBAgAABAgQIECi/gABo+Y31QIAAAQIECBAgQIAAAQIECBAgQIBAhQQEQCsEr1sCBAgQIECAAAECBAgQIECAAAECBMovIABafmM9ECBAgAABAgQIECBAgAABAgQIECBQIQEB0ArB65YAAQIECBAgQIAAAQIECBAgQIAAgfILCICW31gPBAgQIECAAAECBAgQIECAAAECBAhUSEAAtELwuiVAgAABAgQIECBAgAABAgQIECBAoPwCAqDlN9YDAQIECBAgQIAAAQIECBAgQIAAAQIVEhAArRC8bgkQIECAAAECBAgQIECAAAECBAgQKL+AAGj5jfVAgAABAgQIECBAgAABAgQIECBAgECFBARAKwSvWwIECBAgQIAAAQIECBAgQIAAAQIEyi8gAFp+Yz0QIECAAAECBAgQIECAAAECBAgQIFAhAQHQCsHrlgABAgQIECBAgAABAgQIECBAgACB8gsIgJbfWA8ECBAgQIAAAQIECBAgQIAAAQIECFRIQAC0QvC6JUCAAAECBAgQIECAAAECBAgQIECg/AICoOU31gMBAgQIECBAgAABAgQIECBAgAABAhUSEACtELxuCRAgQIAAAQIECBAgQIAAAQIECBAov4AAaPmN9UCAAAECBAgQIECAAAECBAgQIECAQIUEBEArBK9bAgQIECBAgAABAgQIECBAgAABAgTKLyAAWn5jPRAgQIAAAQIECBAgQIAAAQIECBAgUCEBAdAKweuWAAECBAgQIECAAAECBAgQIECAAIHyCwiAlt9YDwQIECBAgAABAgQIECBAgAABAgQIVEhAALRC8LolQIAAAQIECBAgQIAAAQIECBAgQKD8AgKg5TfWAwECBAgQIECAAAECBAgQIECAAAECFRIQAK0QvG4JECBAgAABAgQIECBAgAABAgQIECi/gABo+Y31QIAAAQIECBAgQIAAAQIECBAgQIBAhQQEQCsEr1sCBAgQIECAAAECBAgQIECAAAECBMovIABafmM9ECBAgAABAgQIECBAgAABAgQIECBQIQEB0ArB65YAAQIECBAgQIAAAQIECBAgQIAAgfILCICW31gPBAgQIECAAAECBAgQIECAAAECBAhUSEAAtELwuiVAgAABAgQIECBAgAABAgQIECBAoPwCAqDlN9YDAQIECBAgQIAAAQIECBAgQIAAAQIVEhAArRC8bgkQIECAAAECBAgQIECAAAECBAgQKL+AAGj5jfVAgAABAgQIECBAgAABAgQIECBAgECFBARAKwSvWwIECBAgQIAAAQIECBAgQIAAAQIEyi8gAFp+Yz0QIECAAAECBAgQIECAAAECBAgQIFAhAQHQCsHrlgABAgQIECBAgAABAgQIECBAgACB8gsIgJbfWA8ECBAgQIAAAQIECBAgQIAAAQIECFRIQAC0QvC6JUCAAAECBAgQIECAAAECBAgQIECg/AICoOU31gMBAgQIECBAgAABAgQIECBAgAABAhUSEACtELxuCRAgQIAAAQIECBAgQIAAAQIECBAov4AAaPmN9UCAAAECBAgQIECAAAECBAgQIECAQIUEOlWo36ro9q677oorr7wyTj/99Bg8eHBRY5ozZ05ccsklTdbZa6+9YuDAgfXKjB8/Pu69996YOnVqdOzYMQYMGBD77rtv0WOo17AMAgQIECBAgAABAgQIECBAgAABAgTqCLTbAOjkyZPj17/+dSxYsCDmzZtXB6WQL6+88krceOONTRbdYost6gVAzz///Ljpppuyep07d47FixfH888/H/fcc08cfvjhccghhzTZppsECBAgQIAAAQIECBAgQIAAAQIECBQu0C4DoBMmTIif//znWfCzcKq6JV9++eUsY+utt47hw4fXvfm/3wYNGlQnf9y4cVnwc8UVV4zjjjsuRowYEYsWLcqCnxdffHFceumlMWTIkNhmm23q1POFAAECBAgQIECAAAECBAgQIECAAIHSBNpVADQtW//d734Xt912W6bVoUOHLABZCl0uAJqCn/vss09BTdx///1ZuV133bVOna9//euRlsU/8cQT2dJ4AdCCOBUiQIAAAQIECBAgQIAAAQIECBAgsEyBdnUI0hFHHJEFP7t27Ro/+9nPYr311lsmUGMF0hL4lDbaaKPGitTLf/vtt7O8bbfdtt69HXbYIct766236t2TQYAAAQIECBAgQIAAAQIECBAgQIBAaQLtKgA6a9as2GOPPeLyyy+P3XbbrTSxJbXSvqGvvfZadOrUKdZff/2snY8++miZe4mmPUFTeuihh7LP2n889thj2dfNNtusdrZrAgQIECBAgAABAgQIECBAgAABAgSWQ6BdLYEfPXp0rLHGGsvB9e+qr7/+esyfPz/WWWeduOaaa+KWW26JmTNnRlpSn/IOOuig2H333ev1k/LSyfOPPvpojBkzJgvGpkOQ0gFIafl7mpm6PIHZeh3KIECAAAECBAgQIECAAAECBAgQINDOBWqWBOAWt1eDQw89NF599dW46KKLIjc7sxCLe++9N0aNGpUv2rt371h33XVj6tSpWSA03dhrr73i1FNPzZfJXUyfPj07gOnFF1/MAqYpPx2EtPHGG8cZZ5wRffv2zRUt+jO189577xVdTwUCBAgQIECAAAECBAgQIECAAAEC1Syw2mqrZauxSxlju5oBWgpQQ3Vy+3/26NEjzjrrrHzwNMWS02zQ888/P5vpOWzYsNh5553rNDFx4sTI7fO55pprZsvp33333Uj7g6Z2lycAmjpKQVCJAAECBAgQIECAAAECBAgQIECAAIF/CwiAlvCbcMghh8SIESNi5ZVXrhOwrKmpif333z+bCZoCoVdffXWdAOiPf/zjbPn7JptsEj/96U9j7bXXznpPgc9f/vKXcdppp2X1jz/++BJG9e8qaQwSAQIECBAgQIAAAQIECBAgQIAAAQL/FhAALeE3Ic38HDx4cKM1d91112wm6JQpU7IZmWlv0AkTJmTBz27dumXL59Oy+VwaOHBgnH322XHwwQfHzTffHHvvvXekvGJT6ifNKpUIECBAgAABAgQIECBAgAABAgQIEPi3QLs6Bb6lXnruoKXPP/88Pv3006zbF154Ifvccssto3bwMzemtdZaK3InwE+aNCmX7ZMAAQIECBAgQIAAAQIECBAgQIAAgeUQMAO0BLwbbrghO+xozz33jAEDBtRr4Z133snyevbsGWm2aEoLFizIPldYYYXss6E/UvmUUuBUIkCAAAECBAgQIECAAAECBAgQIEBg+QXMAC3B8IEHHogxY8bEFVdc0WDtxx57LMsfMmRI/v6gQYOy62effbbBg4rS4UXpZPiUSln+nlX0BwECBAgQIECAAAECBAgQIECAAAECdQQEQOtw1P2STma///7748EHH6xzY5dddsm+jx07Nl599dU6955++ulIM0RTOvzww/P3Nt9880hL42fOnBm/+tWvIp0YXztddtll8eabb2aHKg0dOrT2LdcECBAgQIAAAQIECBAgQIAAAQIECJQoYAl8E3BpL85Ro0ZFx44ds1Pfc0VHjhwZjz/+eHaw0aGHHhpbb7119pNOcx83blxW7Ac/+EFstNFGuSrRtWvX7JT3dML7nXfemc32/MIXvhArrbRSjB8/PiZPnhydOnXKyqQ8iQABAgQIECBAgAABAgQIECBAgACB5RcQAC3BMAVEzznnnGwZ/DXXXBNp1mf6Sal///5x7LHHxrBhw+q1nA5AGj16dJx77rkxceLESAHTXNpqq63iRz/6Uayzzjq5LJ8ECBAgQIAAAQIECBAgQIAAAQIECCynQM2Spdh112IvZ4PtrXo63Gj69Okxa9asWH/99aN79+4FEcyZMyemTZsWNTU1WdCzS5cuBdVTiAABAgQIECBAgAABAgQIECBAgACBwgUEQAu3UpIAAQIECBAgQIAAAQIECBAgQIAAgVYm4BCkVvbCDJcAAQIECBAgQIAAAQIECBAgQIAAgcIFBEALt1KSAAECBAgQIECAAAECBAgQIECAAIFWJiAA2spemOESIECAAAECBAgQIECAAAECBAgQIFC4gABo4VZKEiBAgAABAgQIECBAgAABAgQIECDQygQEQFvZCzNcAgQIECBAgAABAgQIECBAgAABAgQKFxAALdxKSQIECBAgQIAAAQIECBAgQIAAAQIEWpmAAGgre2GGS4AAAQIECBAgQIAAAQIECBAgQIBA4QICoIVbKUmAAAECBAgQIECAAAECBAgQIECAQCsTEABtZS/McAkQIECAAAECBAgQIECAAAECBAgQKFxAALRwKyUJECBAgAABAgQIECBAgAABAgQIEGhlAgKgreyFGS4BAgQIECBAgAABAgQIECBAgAABAoULdCq8qJLVLrB48eJIPxIBAgQIECBAgAABAgQIECBAgACBtiRQU1MT6aeUJABailqV1knBz3feeadKR2dYBAgQIECAAAECBAgQIECAAAECBEoT6NOnT3TqVFoos7RapY1TrTILpCh4r169ytyL5gkQIECAAAECBAgQIECAAAECBAi0rECHDqXv5FmzZNagNdMt+770RoAAAQIECBAgQIAAAQIECBAgQIBACwmUHjptoQHqhgABAgQIECBAgAABAgQIECBAgAABAqUKCICWKqceAQIECBAgQIAAAQIECBAgQIAAAQJVLyAAWvWvyAAJECBAgAABAgQIECBAgAABAgQIEChVQAC0VDn1CBAgQIAAAQIECBAgQIAAAQIECBCoegEB0Kp/RQZIgAABAgQIECBAgAABAgQIECBAgECpAgKgpcqpR4AAAQIECBAgQIAAAQIECBAgQIBA1QsIgFb9KzJAAgQIECBAgAABAgQIECBAgAABAgRKFRAALVVOPQIECBAgQIAAAQIECBAgQIAAAQIEql5AALTqX5EBEiBAgAABAgQIECBAgAABAgQIECBQqoAAaKly6hEgQIAAAQIECBAgQIAAAQIECBAgUPUCAqBV/4oMkAABAgQIECBAgAABAgQIECBAgACBUgUEQEuVU48AAQIECBAgQIAAAQIECBAgQIAAgaoXEACt+ldkgAQIECBAgAABAgQIECBAgAABAgQIlCrQqdSKbaHeXXfdFVdeeWWcfvrpMXjw4KIf6f33349rrrkmXn755fjss89iyJAhseWWW8bOO+/cZFvjxo2LRx99NN58881YtGhRrLPOOrHddtvFbrvt1mQ9NwkQIECAAAECBAgQIECAAAECBAgQKE6gZvGSVFyVtlF68uTJceyxx8aCBQvioosuii222KKoB3vyySfjJz/5ScydOzer16NHj/j444+z66985Stx8sknR8eOHeu0OW/evDjppJNiwoQJWX7Pnj2zz48++ij7TGM455xzokuXLnXq+UKAAAECBAgQIECAAAECBAgQIECAQGkC7XIJfApAnnbaaVnwsxS2d999N84444ws+HnAAQfE9ddfH3fccUf86le/igEDBsTdd98dv//97+s1/dvf/jYLfqYyl156aaQZqOnnkksuif79+8fEiROzYGy9ijIIECBAgAABAgQIECBAgAABAgQIEChJoF0FQOfMmRO//vWvs5mfH374YXToUNrj33vvvdlsz379+sUxxxwTffv2zWZ7Dhs2LA4++ODsRTz00EN1Xkjq+/bbb8/6/MUvfhEbbbRR/v7GG28co0aNyr7feeedkcpKBAgQIECAAAECBAgQIECAAAECBAgsv0BpEcDl77ciLRxxxBFx2223RdeuXeNnP/tZrLfeeiWNIy1333bbbbNg59JB1B122CFr85133okUZM2ltOR+4cKF2UzPhvpNeX369Im0I8Frr72Wq+aTAAECBAgQIECAAAECBAgQIECAAIHlEGhXhyDNmjUr9thjjzj88MNjrbXWijFjxpREt99++0X6aSjlgpdpVugqq6ySL5ICpmkGaG7P0PyN/71Ie5HOnj07+9arV6+lb/tOgAABAgQIECBAgAABAgQIECBAgEAJAu0qADp69OhYY401SmBadpV0wFHaW/T888/PCi8dIK2pqakTEF26xfvuuy8+//zzWHnllSMtrS8lpRPlc0HUUuqrQ4AAAQIECBAgQIAAAQIECBAgQKAaBdJh4ksfOF7oONtVALRcwc+0r+g999yTBTDTi0gHLO25556FvoOYPn16/tCk//zP/4wULC01NTbDtNT21CNAgAABAgQIECBAgAABAgQIECBQaYG0JWWpqV3tAVoq0rLqPfPMM9G7d+/o1KlTts9nOgBpxowZy6qW3Z85c2accMIJkZbnp2Xy++yzT0H1FCJAgAABAgQIECBAgAABAgQIECBAYNkC7WoG6LI5SiuR9hJNszbnz58fV155ZVx++eUxfvz4OO+882Lo0KGNNjpt2rQ48cQT4+23345NNtkk0unwy5PSGFIgViJAgAABAgQIECBAgAABAgQIECDQlgRKXf6eDARAm+E3IbdkfYUVVsgOWEpL2tOenldccUWk5fENpWeffTZOPfXU+Pjjj2ObbbaJM888M7p169ZQ0YLz0jjSGCQCBAgQIECAAAECBAgQIECAAAECBP4tYAl8GX4Tdtlll6zVl156qcHWx44dGz/84Q+z4Gc6lf5Xv/rVcgc/G+xIJgECBAgQIECAAAECBAgQIECAAIF2LmAGaAm/AGlm5xtvvBFHHnlkg6fK52Zhpj1Bl0533HFHnHPOOVn2d7/73TjssMOWLuI7AQIECBAgQIAAAQIECBAgQIAAAQLNJGAGaAmQTz31VPz1r3+Nhx9+uMHaEyZMyPI33HDDOvf//ve/Z7M901L1tPxd8LMOjy8ECBAgQIAAAQIECBAgQIAAAQIEml1AALQJ0nQ40f333x8PPvhgnVIjRozIvqeZoG+99Vadeyn4ef3112d5++23X/7evHnzskORFi9eHEcccUTstdde+XsuCBAgQIAAAQIECBAgQIAAAQIECBAoj0D9Ndrl6adVtjpp0qQYNWpUpFOmckHP9CApsJlmcz7xxBPx7W9/OzvEaMstt4xXXnkl0v6eKcj5jW98I4YNG5Z/7htvvDHS4UgpXXbZZdlP/uZSF+lApC9/+ctL5fpKgAABAgQIECBAgAABAgQIECBAgECxAgKgxYotKZ+WsKcg5bXXXhtXXXVVpCXx6Sel1VdfPY455pjYeeeds++5P1IwNZcWLlyYu2zwc9GiRQ3myyRAgAABAgQIECBAgAABAgQIECBAoDiBmiWzFRcXV0Xp2gILFiyIN998Mz788MNYZ511YrXVVqt92zUBAgQIECBAgAABAgQIECBAgAABAhUUEACtIL6uCRAgQIAAAQIECBAgQIAAAQIECBAor4BDkMrrq3UCBAgQIECAAAECBAgQIECAAAECBCooIABaQXxdEyBAgAABAgQIECBAgAABAgQIECBQXgEB0PL6ap0AAQIECBAgQIAAAQIECBAgQIAAgQoKCIBWEF/XBAgQIECAAAECBAgQIECAAAECBAiUV0AAtLy+WidAgAABAgQIECBAgAABAgQIECBAoIICAqAVxNc1AQIECBAgQIAAAQIECBAgQIAAAQLlFRAALa+v1gkQIECAAAECBAgQIECAAAECBAgQqKCAAGgF8XVNgAABAgQIECBAgAABAgQIECBAgEB5BQRAy+urdQIECBAgQIAAAQIECBAgQIAAAQIEKiggAFpBfF0TIECAAAECBAgQIECAAAECBAgQIFBeAQHQ8vpqnQABAgQIECBAgAABAgQIECBAgACBCgoIgFYQX9cECBAgQIAAAQIECBAgQIAAAQIECJRXQAC0vL5aJ0CAAAECBAgQIECAAAECBAgQIECgggICoBXE1zUBAgQIECBAgAABAgQIECBAgAABAuUVEAAtr6/WCRAgQIAAAQIECBAgQIAAAQIECBCooIAAaAXxdU2AAAECBAgQIECAAAECBAgQIECAQHkFOpW3ea23pMDixYvjk08+acku9UWAAAECBAgQIECAAAECBAgQIECg7ALdunWLDh1Km8spAFr219NyHaQA6KefftpyHeqJAAECBAgQIECAAAECBAgQIECAQAsIdOnSpeQAaM2SoNniFhijLggQIECAAAECBAgQIECAAAECBAgQINDiAqXNG23xYeqQAAECBAgQIECAAAECBAgQIECAAAECxQsIgBZvpgYBAgQIECBAgAABAgQIECBAgAABAq1EQAC0lbwowyRAgAABAgQIECBAgAABAgQIECBAoHgBAdDizdQgQIAAAQIECBAgQIAAAQIECBAgQKCVCAiAtpIXZZgECBAgQIAAAQIECBAgQIAAAQIECBQvIABavJkaBAgQIECAAAECBAgQIECAAAECBAi0EgEB0FbyogyTAAECBAgQIECAAAECBAgQIECAAIHiBQRAizdTgwABAgQIECBAgAABAgQIECBAgACBViIgANpKXpRhEiBAgAABAgQIECBAgAABAgQIECBQvIAAaPFmahAgQIAAAQIECBAgQIAAAQIECBAg0EoEBEBbyYsyTAIECBAgQIAAAQIECBAgQIAAAQIEihcQAC3eTA0CBAgQIECAAAECBAgQIECAAAECBFqJgABoK3lRhkmAAAECBAgQIECAAAECBAgQIECAQPECAqDFm6lBgAABAgSaRWDurFmxaNGiZmlLIwQIECBAgAABAgQIECDQsIAAaMMucgkQIECAQFkEptz/QDy84SYxo0OX6LTKmjGvY5eY2HO1eOTIo2LhggVl6VOjBAgQIECAAAECBAgQaM8CNYuXpPYM4NkJECBAgEBLCTx2zHGx5cV/iM5RE4uX/FNT6zON4ZlVe8XgyROiR9+1WmpI+iFAgAABAgQIECBAgECbFxAAbfOv2AMSIECAQDUITLjgohj8wxOj45LBpMBnY+kfa/SJ7Wa80dht+QQIECBAgAABAgQIECBQpIAl8EWCKU6AAAECBEoRWOmU06JTFvpsPPiZZoV+4Z334qlfnFlKF+oQIECAAAECBAgQIECAQAMCAqANoMgiQIAAAQLNKfDi9X+JQfM+X2aTuZmh8y67fJllFSBAgAABAgQIECBAgACBwgQ6FVasbZa666674sorr4zTTz89Bg8eXPRDzps3L2688cYYP358fPjhhzFo0KDYYostYs8994yOHdMix4bTiy++GDfccEO8/vrr0a1btxg6dGgMHz481l9//YYryCVAgACBVi3w3tiHYmART9BnxrtFlFaUAAECBAgQIECAAAECBJoSaLcB0MmTJ8evf/3rWLDkxN0UyCw2zZo1K37wgx/EG2/8e5+2VVddNe69997s5/HHH4+f//znseKKK9ZrNgVML7jggiy/e/fu8fnnn8czzzwTf/nLX+Lss8+Orbbaql4dGQQIECDQugUWL/l3fTGp46JFxRRXlgABAgQIECBAgAABAgSaEGiXS+AnTJgQp512Whb8bMKmyVu//OUvs+DnF7/4xbjzzjvjtttui+uuuy422GCDeOSRR+LCCy+sVz8FXVN+CoyOGjUq7r777ixgeuyxx8Znn30WJ554YsyYMaNePRkECBAg0LoFVv7C1kU9wPur9CqqvMIECBAgQIAAAQIECBAg0LhAuwqAzpkzJ5v1mQKOacl6hw6lPf7zzz8fTz31VHTp0iXOPPPMWHnllTPhfv36xbnnnpstf7/nnnvi448/riN/xRVXxOLFi+M73/lO7LjjjlFTUxMrrLBCjBw5Mg444ICYP39+3HrrrXXq+EKAAAECrV9gyBGHxzsdapYccbS4oIeZv89XCyqnEAECBAgQIECAAAECBAgsW6C0COCy263KEkcccUQ2U7Nr167xs5/9LNZbb72SxvnQQw9l9XbaaadYaaWV6rSRlsJvu+222dL2FATNpRR8TUHTlPbYY49cdv4zl5dmk6Zl+RIBAgQItB2BTktm/r92xHezM+CXFQT910qd40u/rb+KoO1oeBICBAgQIECAAAECBAi0rEC7CoCmfTtToPHyyy+P3XbbrWTpf/7zn1ndtPy9oZQCoCk9++yz+dsvvPBCNvuzf//+0bdv33x+7mLjjTeOHj16xOzZs2PatGm5bJ8ECBAg0EYEdvjj7+KR7bbNgqBLP1IuKDqtU8fo+cA9kQKmEgECBAgQIECAAAECBAg0j0C7OgRp9OjRscYaayy33FtvvZW10atXw3u05fJzBySlwsuqk8qkemnZfKpXyonwi5YcmpGW9ksECBAgUJ0Cm99xSzwy6r9j9d9fEhvP/f8H8H28JCz69KaDY+AVf44u664bM2fOrM4HMCoCBAgQIECAAAECBAhUSCDFzTp27FhS7+0qANoc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8jM6dIx+8z+JDh38t+NlgilAgAABAgQIECBAoMoE/C2+yl6I4RAgQIAAAQItI/DxQd/OOlrWPqBpmfzLu+4s+Nkyr0UvBAgQIECAAAECBJpdwAzQZifVIAECBAgQINAaBNJBSH9fe0Bs+/a72WL4pfcDTYHRlPfP7l1j47emRueePVvDYxkjAQIECBAgQIAAAQJLCZgBuhSIrwQIECBAgED7EEjL2bd+8bn428D1s0Dn0k+dgp//WKNPrPvcBMHPpXF8J0CAAAECBAgQINCKBMwAbUUvy1AJECBAgACB8gj888qr453zL4zur0+LDosXx8d914peR3w3tvzhceXpUKsECBAgQIAAAQIECLSYQLsLgL744otxww03xOuvvx7dunWLoUOHxvDhw2P99dcvCn3OnDlxySWXNFlnr732ioEDB9YrM378+Lj33ntj6tSp0bFjxxgwYEDsu+++MXjw4HplZRAgQIAAAQIECBAgQIAAAQIECBAgULpAuwqA3njjjXHBBRdkWt27d4/PP/88++nSpUucffbZsdVWWxUs+eyzz8bRRx/dZPkzzzwzdtpppzplzj///LjpppuyvM6dO8fiJbNM0jhqamri8MMPj0MOOaROeV8IECBAgAABAgQIECBAgAABAgQIEChdoFPpVVtXzcmTJ8eFF14YK664Yvz85z+PHXbYIRYsWBC33nprln/iiSfGNddcE2uuuWZBD/byyy9n5bbeeutsBmlDlQYNGlQne9y4cVnwM43huOOOixEjRkQ6gOGee+6Jiy++OC699NIYMmRIbLPNNnXq+UKAAAECBAgQIECAAAECBAgQIECAQGkC7SYAesUVV2SzLb/zne/EjjvumGmtsMIKMXLkyJg+fXqk2aEpGHrUUUcVJJkLgKbl8/vss09Bde6///6s3K677lqnzte//vVIy+KfeOKJbGm8AGhBnAoRIECAAAECBAgQIECAAAECBAgQWKZAuzgFPu3X+dRTT2UYe+yxRz2UXN6dd96ZzQqtV6CBjFdeeSXL3WijjRq423DW22+/nd3Ydttt6xVIM1JTeuutt+rdk0GAAAECBAgQIECAAAECBAgQIECAQGkC7SIA+sILL2SzP/v37x99+/atJ7XxxhtHjx49Yvbs2TFt2rR695fOSEvnX3vttejUqVP+8KSPPvoo5s2bt3TROt+32GKL7PtDDz1UJz99eeyxx7K8zTbbrN49GQQIECBAgAABAgQIECBAgAABAgQIlCbQLpbA52ZV9urVq1GldO/jjz+ON954Ix/UbKxwOkF+/vz5sc4662T7ht5yyy0xc+bM6NChQ5Z30EEHxe67716vesq766674tFHH40xY8ZEmnmaDkFKe4Cm5e9du3aN3XbbrV69YjLSnqISAQIECBAgQIAAAQIECBAgQIAAgbYkkA4QTz+lpHYRAP30008zm6YCoD179szK5Mo2hZnb/zPNFk0HF/Xu3TvSYUhTp07Nfn75y1/GM888E6eeemqdZgYPHhyXX355dgjTH/7wh/jTn/6U3U9ByzQL9YwzzmhwhmqdRpr4ktp55513mijhFgECBAgQIECAAAECBAgQIECAAIHWJ9CnT59sNXYpI28XAdC0B2hKaZl7Y6l79+7Zrblz5zZWJJ+f2/8ztXfWWWdFbml7ms2ZZoOef/752UzPYcOGxc4775yvly4mTpyY3+cznTifltO/++67kfYHTe02tES/TgO+ECBAgAABAgQIECBAgAABAgQIECBQsEC7CIB269YtA/n8888bhcnt39m5c+dGy+RuHHLIITFixIhYeeWV6wQs0zTc/fffP5sFmgKhV199dZ0A6I9//ONs+fsmm2wSP/3pT2PttdfOmkyBzzRr9LTTTsvqH3/88bmuiv5MJ9tLBAgQIECAAAECBAgQIECAAAECBNqSQKnL35NBuwiApiXqKaWDihpLaf/PlHLB0sbKpfw08zMtZ28s7brrrtlM0ClTpkRalp72Bp0wYUIW/Eztjxo1Kls2n6s/cODAOPvss+Pggw+Om2++Ofbee+9IecWm1E/uWYutqzwBAgQIECBAgAABAgQIECBAgACBtijQLk6BzwUFc0HOhl5kLji6yiqrNHS7qLw11lgjK59mnOb2FE0n0ae05ZZbNhikXGuttSJ3AvykSZOysv4gQIAAAQIECBAgQIAAAQIECBAgQGD5BNrFDNDVV189U0onvKfT25deJj579uz44IMPspmagwYNWqboDTfckJ36vueee8aAAQPqlc8dRJQOVsrtO5r2+kxp6b5rV84dxNTUUv3a5V0TIECAAAECBAgQIECAAAECBAgQINC0QLuYAZoOFkqnrH/yySfx5JNP1hMZN25cLFy4MCvTtWvXeveXznjggQdizJgxccUVVyx9K/v+2GOPZZ9DhgzJ388FVp999tlsWXz+xv9epKXyL774YvatlOXvS7fnOwECBAgQIECAAAECBAgQIECAAAECEe0iAJpe9Le+9a3sfY8ePTpqL4VPJ7Bfe+212b2RI0dmn7k/0sns999/fzz44IO5rOxzl112yT7Hjh0br776ap17Tz/9dKQZoikdfvjh+Xubb755pKXxM2fOjF/96leRToyvnS677LJ48803s0OVhg4dWvuWawIECBAgQIAAAQIECBAgQIAAAQIEShSoWRKIqxuJK7Ghaq+WZnh+//vfj7QXZ5oRmoKYaVl6ms2ZgpLbb799nHXWWdky+Nyz3HvvvdmBRR07doyHHnool53NFk0ntaeDjVLaeuuts590mnuaTZpS6isXdM0ylvyRyqd6aSxplucXvvCFWGmllWL8+PExefLk6NSpU1xwwQX5vUBz9XwSIECAAAECBAgQIECAAAECBAgQIFCaQLsJgCaeefPmxXnnnRf33XdfthdoykvBzf322y++973vZcHIlJdLjQVA0/25c+dmy+CvueaaqL1nZ//+/ePYY4+NYcOG5Zqp85lOhj/33HNj4sSJdfK32mqr+NGPfhTrrLNOnXxfCBAgQIAAAQIECBAgQIAAAQIECBAoXaBdBUBzTGnmZ1q6nia/poBlt27dcreK/kxtTZ8+PWbNmhXrr79+dO/evaA25syZE9OmTYuampos6NmlS5eC6ilEgAABAgQIECBAgAABAgQIECBAgEDhAu0yAFo4j5IECBAgQIAAAQIECBAgQIAAAQIECLRmgXZzCFJrfknGToAAAQIECBAgQIAAAQIECBAgQIBAaQKdSqumFgECBAgQIECAQGsWWLhkG5/nLvlzzPr73yNqOsQqX9ouhhx2aHRcciijRIAAAQIECBAgQKAtCVgC35bepmchQIAAAQIECBQg8MTJp0bfcy+MfgsX1Sn9RqeO8c5Jx8ews86sk+8LAQIECBAgQIAAgdYsIADamt+esRMgQIAAAQIEihQYt/d+scOd9zRZ69H99o5dbr6hyTJuEiBAgAABAgQIEGgtAvYAbS1vyjgJECBAgAABAsspMOn3f4wv33l3LF7yT2Mp3dv+lttj8qV/bqyIfAIECBAgQIAAAQKtSsAM0Fb1ugyWAAECBAgQIFC6wFO914qtZn5YUAPj+6wWw959q6CyChEgQIAAAQIECBCoZgEB0Gp+O8ZGgAABAgQIEGgmgc8++CBqVlsrOi5pr2bJP02lNAt0/pICHWe/F5179myqqHsECBAgQIAAAQIEql7AEviqf0UGSIAAAQIECBBYfoH3np0cnbLQZ9PBz9RTCpCuuOTn/X8+v/wda4EAAQIECBAgQIBAhQUEQCv8AnRPgAABAgQIEGgJgS59+hTdTZfVi69TdCcqECBAgAABAgQIECizgCXwZQbWPAECBAgQIECgWgSmrtAt+i1YuGRu57LT6yt0ig0+/2TZBZUgQIAAAQIECBAgUOUCZoBW+QsyPAIECBAgQIBAcwm8NmLngoKfqb/Xdx/RXN1qhwABAgQIECBAgEBFBcwArSi/zgkQIECAAAECLScwd9aseKXfgNh4ztwmO32+W5fYaPrrDkBqUslNAgQIECBAgACB1iJgBmhreVPGSYAAAQIECBBYToGVevWKNcf/PSb1avxk94mrrBxrT3hK8HM5rVUnQIAAAQIECBCoHgEzQKvnXRgJAQIECBAgQKBFBBYtWhRPnvrjWHjtX6L3e+/H4iW9zly9d3Q88Fux3X+PapEx6IQAAQIECBAgQIBASwkIgLaUtH4IECBAgAABAgQIECBAgAABAgQIEGhxAUvgW5xchwQIECBAgAABAgQIECBAgAABAgQItJSAAGhLSeuHAAECBAgQIECAAAECBAgQIECAAIEWFxAAbXFyHRIgQIAAAQIECBAgQIAAAQIECBAg0FICAqAtJa0fAgQIECBAgAABAgQIECBAgAABAgRaXEAAtMXJdUiAAAECBAgQIECAAAECBAgQIECAQEsJdGqpjvRDgAABAgQIECBQPQKfffZZPP744zFjxoxYbbXVYrvttouVV165egZoJAQIECBAgAABAgSaSUAAtJkgNUOAAAECBAgQaA0Cc+fOjdNPPz0uvPDCSEHQXOrUqVMceuihcc4558Qqq6ySy/ZJgAABAgQIECBAoNUL1Cxeklr9U3gAAgQIECBAgACBZQrMmTMndt1113jiiScaLTtw4MB49NFHY80112y0jBsECBAgQIAAAQIEWpOAPUBb09syVgIECBAgQIDAcggcc8wxTQY/U9OvvPJKfPvb316OXlQlQIAAAQIECBAgUF0CZoBW1/swGgIECBAgQIBAWQReffXVGDRoUBS6+OfBBx+M4cOHl2UsGiVAgAABAgQIECDQkgJmgLaktr4IECBAgAABAhUSuOuuuwoOfqYh3n777RUaqW4JECBAgAABAgQINK+AAGjzemqNAAECBAgQIFCVAlOmTClqXFOnTi2qvMIECBAgQIAAAQIEqlVAALRa34xxESBAgAABAgSaUaBbt25Ftda1a9eiyitMgAABAgQIECBAoFoFBECr9c0YFwECBAgQIECgGQW23nrrolrbZpttiiqvMAECBAgQIECAAIFqFXAIUrW+GeMiQIAAAQIECDSjwNy5c2ODDTaI6dOnL7PVLl26RDo0aa211lpm2eYo8NFHH8V1110Xf/vb3+LDDz+Mfv36xe677x777LNPdOzYsTm60AYBAgQIECBAgEA7FhAAbccv36MTIECAAAEC7UsgHWz0ta99bZkPff7558dxxx23zHLNUeCWW26JI488MmbOnFmvuU033TSuv/762GSTTerdk0GAAAECBAgQIECgUAEB0EKllCNAgAABAgQItAGBq666Kgs4zps3r97T1NTUxKhRo+K//uu/6t0rR8aNN94YI0eObLLpXr16xT/+8Y8YOHBgk+XcJECAAAECBAgQINCYgABoYzLyCRAgQIAAAQJtVGDatGlx4YUXxgMPPBAzZsyIVVddNXbcccc45phjYsiQIS3y1B988EGsv/76MXv27GX2l8b28MMPL7OcAgQIECBAgAABAgQaEhAAbUhFHgECBAgQIECAQFkFzjvvvDjhhBMK7uOZZ56JLbfcsuDyChIgQIAAAQIECBDICTgFPifhkwABAgQIECBAoMUEHnnkkaL6KrZ8UY0rTIAAAQIECBAg0KYFBEDb9Ov1cAQIECBAgACB6hRIS+CLScWWL6ZtZQkQIECAAAECBNq2gABo236/no4AAQIECBAgUJUCa665ZlHjWmONNYoqrzABAgQIECBAgACBnIAAaE7CJwECBAgQIECAQIsJ7LrrrkX1VWz5ohpXmAABAgQIECBAoE0LCIC26dfr4QgQIECAAAEC1Slw4IEHRt++fQsa3L777hsbbrhhQWUVIkCAAAECBAgQILC0gADo0iK+EyBAgAABAgQIlF2ga9eucd1110Xnzp2b7GvdddeNP/7xj02WcZMAAQIECBAgQIBAUwI1i5ekpgq413oEFi1aFA4IaD3vy0gJECBAgACBiKeffjqOPfbY+Ne//lWPIy17v+CCC2L11Vevc+/JJ5+Mq666KiZMmBCffvpp9O/fP3bffff47ne/G927d69T1hcCBAgQIECAAIG2IbDKKqtEx44dS3oYAdCS2KqzUgqAzp49uzoHZ1QECBAgQIAAgUYEFi5cGI888kikwGb6u0w6IGn48OExdOjQOjXmz58fp5xySlx++eV18nNf+vTpkwVGt91221yWTwIECBAgQIAAgTYi0LNnTwHQNvIuPQYBAgQIECBAgEAjAkcdddQyl8OnGaBPPfVUDB48uJFWZBMgQIAAAQIECLQ3ATNA29sb97wECBAgQIAAgVYo8Pjjj8f2229f0MjT7NEHH3ywoLIKESBAgAABAgQItH0BhyC1/XfsCQkQIECAAAECrV7g0ksvLfgZxo4dG6+99lrB5RUkQIAAAQIECBBo2wICoG3o/X744YfZvle33XZbzJgxow09mUchQIAAAQIE2rvAP/7xj6IIxo8fX1R5hQkQIECAAAECBNqugABoG3q3aanXwQcfHPvuu2+stdZaccABB8Qbb7zRrE/40ksvxRFHHJGdtppO3lpttdWy/h544IFm7UdjBAgQIECAAIHaAum092JSseWLaVtZAgQIECBAgACB1iUgANq63leTo33//ffr3L/ppptim222iRdffLFOfqlf/vSnP8Wmm24af/7zn+PNN9+MdOr8Bx98EGnG6W677Rbf+973Ip3iKhEgQIAAAQIEmltg3XXXLarJYssX1bjCBAgQIECAAAECrUpAALRVva6mB7vBBhvUK/Duu+9mMzTnz59f714xGbfccksW4FywYEGj1VKA9NRTT230vhsECBAgQIAAgVIF9t5774Kr9urVK7785S8XXF5BAgQIECBAgACBti0gANq232/2dGnZ+jXXXFPyk86bNy+OOeaYguqfe+658cILLxRUViECBAgQIECAQKECaaVJ3759Cyr+k5/8JFZcccWCyipUX+Czzz6LKVOmxMyZM+vflEOAAAECBAgQaIUCAqCt8KWVMuQ77rijlGpZnbS36FtvvVVQ/bQs/uqrry6orEIECBAgQKCtCkyfPj3OOOOMGDFiRGy55ZbxH//xH/Gb3/wmZs+e3VYfuezP1a1bt0grUnr27NlkX9/85jfjhBNOaLKMmw0LPPTQQ7HrrrtGjx49Yv3114/evXvHhhtuGOedd158/vnnDVeSS4AAAQIECBBoBQICoK3gJTXHEKdOnVpyMxMnTiyqbrHli2pcYQIECBAgUOUCo0ePjoEDB8bpp58eY8eOjfS/i/fee2+ceOKJWVApXUulCWy77baRToPfY4896jWwyiqrZEHmtOqlpqam3n0ZTQuceeaZscsuu0T6D9+193R/+eWXs4Dy8OHDY9asWU034i4BAgQIECBAoEoFOlXpuAyrmQXSrIlSU7H/xT8tmZcIECBAgEB7FEirIA477LBGHz0dHpj2skxBph133LHRcm40LpBmJKYg8rRp07Jg6CeffBLpwKMvfelLlr03ztbknTFjxsRPf/rTJsv87W9/iwMPPDDuuuuuJsu5SYAAAQIECBCoRgEB0Gp8K2UYUzoNvtTU0OFKTbWVZr1IBAgQIECgvQmk4ObRRx+9zMdOBwqmIOmLL74YnTr5q9gywRopsM4660T6kZZPIP2H7jQ7uZB09913x3333Re77757IcWVIUCAAAECBAhUjYC/dVfNqyjfQNIysKZmoyyr57RvWefOnaPQmZ377bffspp0nwABAgQItDmBa6+9Nj766KOCnuvVV1+NBx54IPbcc8+CyremQmmf0zRDMx3C2LFjx9h0002zgFmXLl1a02O0m7E+/PDDMWPGjIKfN/2eC4AWzKUgAQIECBAgUCUCAqBV8iKaYxhz585tsJlTTjklhgwZ0uC9QjLTBvgnnXRSpL2hlpXS3lEN7cu1rHruEyBAgACB1i7wxBNPFPUIqXxLBEDTfo7XX3993H777ZH2BO/atWuklSGHHnpobLLJJkWNuanCixcvzvbgTIc/pWXptVOfPn2yewcddFDt7OW+Tocvpu0EUhAvnVi+xhprZAdP7bDDDsvddntp4IUXXijqUYstX1TjChMgQIAAAQIEyiQgAFom2Eo0O2XKlHrdplNQR40aVS+/2Ix0kEP6C+9NN93UaNUUZL3uuusave8GAQIECBBoywKFzv7MGRRbPlevmM90gM3+++8fzz33XJ1q48aNywKSJ598cvb3hA4dlv9czKOOOir+9Kc/1ekn9+W9996Lgw8+ON5666049dRTc9nL9Tlp0qRIAdXJkyfXaScFYLfbbru46qqrothtfOo05AsBAgQIECBAgECbEVj+v+22GYrW/yDDhg2Lfv36xeDBg7Ml70899VT2f26a4//UpCVsN9xw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OHBg1HO66667XKqpqJ1YiUAEAb+px0488cQIe2IxAggggAAC+UMg7gDod999565QT6CV6zNc69y5s02cONHKli3rgqAaCTpnzpxwXVmGAAIIIIAAAggggECBF1ixYoU9+OCDnhxuueUW02wrGgJ+BVSwtnXr1p42O/TQQ+2MM87w1JdOCCCAAAIIpKpA3AHQJUuWuGtq3Lhx1Gtr27atvfvuu1a0aFHbsmWLC5auXLky6jasRAABBBBAAAEEEECgIAq8/vrrbuCAl2tX+ily7XuRok84gZEjR1qVKlXCrcpcplRnr732GgW4MkV4gQACCCCQXwXiDoCWLl3aXfOff/4Z89rbtWtnw4cPNxVJWrNmjXXs2NEFQ2NuSAcEEEAAAQQQQAABBAqQwNy5c31d7bx583z1pzMCAYHatWvbl19+aRqwEq4dffTR9sUXX1isAS/htmUZAggggAACqSZQON4Tql+/vttUVSh3794d86ngv//9b9OUnjvuuMN0o6bp8O+9954bGRrvObAdAggggAACCCCAAALpJLBz505fl+O3v6+d0zntBWrVqmXTp093gdApU6bYb7/9ZpUqVXJBUU2R1wAWGgIIIIAAAukgEHcAtFmzZu4X4tKlS+2JJ56wm2++OabH7bffbsuXL3cJ21UQ6dxzz7VRo0bF3I4OCCCAAAIIIIAAAggUBIHAIAOv1+q3v9f90q9gCbRo0cL0h4YAAggggEC6CsQ9BV43WxdeeKFzUQL2Hj162KeffmobN26MavXMM89kFk1SbtCWLVvanj17om7DSgQQQAABBBBAAAEECoKAioh6bYUKFbIzzzzTa3f6IYAAAggggAACBVYg7gCoxPr27WsHHnigw3v55ZddJcHu3btHxdSN2tixYzNv1hYuXOim0GujvXv3Rt2WlQgggAACCCCAAAIIpLPAUUcdZeedd56nS7z22mtNU5hpCCCAAAIIIIAAAtEFEgqAKvg5Z84cO+OMMzKPUqNGjczXkV4UK1bMBUH79etnJUuWjNSN5QgggAACCCCAAAIIFDiBoUOH2rHHHhv1uk8//XR75JFHovZhJQIIIIAAAggggMA/AgkFQLULJckeN26cLVu2zFV679SpkyfbwoUL26233mqLFy+2rl27WtmyZT1tRycEEEAAAQQQQAABBNJZoEyZMjZ16lS7//773b128LVWr17dBgwY4O6/ixYtGryK1wgggAACCCCAAAIRBDL2TTtPmXnn27dvtxIlSkQ4VRYjgAACCCCAAAIIIFCwBHbv3u0GDKxfv94OOOAAa9iwIZW5C9ZfgRy5WlV7X7Jkifu7dMghh1jlypVz5DjsFAEEEEAAgVQRSKkAaKqgcB4IIIAAAggggAACCCCAQLoJfP3116YCthphHGgZGRl2yimnuJQKjRs3DizmJwIIIIAAAmklQAA0rT7O7Bcza9YsmzFjhm3YsMGqVatmJ554oh166KHZO7IEAQQQQAABBBBAAAEE0lbg9ddftwsvvNB27doV9hqLFy9ub7zxhnXs2DHsehYigAACCCCQnwWSHgD9+OOPbfXq1b5MVBleOUDLlSvnpvbUq1fP1/Z0zi6waNEiu/TSS00B0NCmolVKrq9pVDQEEEAAAQQQQAABBBBIb4EFCxbYUUcdZTt37ox6oSpQO2/ePOP7WFQmViKAAAII5EOBpAdAO3ToYBMmTEiIQtXlu3TpYg8++KCVKlUqoX0VxI3nzp1rxx13nG3ZsiXi5deqVcsFRwmCRiRiBQIIIIAAAggggAACaSHQuXNne/vttz1dS48ePWzEiBGe+tIJAQQQQACB/CKQkgHQAJ4Scn/++edWoUKFwCJ+xhD4+++/3RT3H374IUZPMwWr33///Zj96IAAAggggAACCCCAAAL5U2DHjh1upl2s0Z+BqytTpoxt2rTJ9ttvv8AifiKAAAIIIJDvBZL+W+3dd981/Qn8wmzdurWNHj3avvzyS/vll19ctcHJkydbnz59rHz58g6wQYMGNn78ePdU8qmnnsqccvH9999bz5498z1ybl7AW2+9ZV6CnzonmWuKCw0BBBBAAAEEEEAAAQTSU2DlypUxp74HX7lmka1bty54Ea8RQAABBBDI9wJJD4AuXbrULrroItuzZ48NGzbMZs6cad26dbPmzZu7nJP169d3hXjuvfde++mnn9zyJUuW2Jw5c+yss86y3r17u6DcNddc43DHjh1ry5cvz/fQuXUBkyZN8nUoBaNpCCCAAAIIIIAAAgggkJ4CRYsW9X1h8Wzj+yBsgAACCCCAQC4KJD0A2r9/f9u8ebNdddVVrghPtGvRCNAxY8aYfsHefffd9vvvv7vuSr49aNAga9iwoe3du9dzvppoxyoo69asWePrUv3297VzOiOAAAIIIIAAAggggECeCqi+QsWKFT2fQ82aNUlB5lmLjggggAAC+UUg6QFQ5exUO++88zwZ1K5d25o2bepGjGqafKBlZGTYCSec4N4yAjSgEvun33ypfvvHPgN6IIAAAggggAACCCCAQKoIFCpUyLp37+75dFSMloYAAggggEC6CSQ1AKoCPMrbqVajRg3PVlWqVHF9Fy9enGWbSpUquffr16/Pspw3kQVU/d1P89vfz77piwACCCCAAAIIIIAAAnkvoNl2GgnqpQ0YMMBatmxpEyZM8NKdPggggAACCOQLgaQGQAsXLmyVK1d2F/7VV195Bvj6669dX015D25ffPGFe1u1atXgxbyOInD++ed7nuJyxBFHGAHQKJisQgABBBBAAAEEEEAgDQT0HW3ixImuJoOXy9HMvA4dOrg0ZV760wcBBBBAAIFUF0hqAFQXe9RRR7lrvueee2zTpk0xr//RRx+1tWvXun6BbfVmx44droCSXmuaPM2bQLly5eyFF16I2blEiRL24osv2n77Jf2vQMxj0wEBBBBAAAEEEEAAAQRyV0B5QPUdy0978MEHbeTIkX42oS8CCCCAAAIpKZD06NeNN97oLlTV4Nu3b28fffRR2Atft26d3XrrrXbHHXe49aecckrmE8lFixa5yvHbt2835azR00ead4Gzzz7bxo4dGzF5+UEHHWTTpk2zZs2aed8pPRFAAAEEEEAAAQQQQCDfCvTr188Vq/V7AfrOplRnNAQQQAABBPKzQMa+Kut7k30BvXv3tsGDB2futk6dOla3bl2Xd2bDhg22cuVKlytUAU61+vXr26xZszIDdirMExg9esEFF9ioUaMy98UL7wIy1BPbmTNnmtyVSuCkk05yweVixYp53xE9EUAAAQQQQAABBBBAIF8LKAfo6tWr47qG6dOnW9u2bePalo0QQAABBBBIBYEcCYAqpnr//febnjLGmmZx6qmnumBpvXr1nIdGhgaKImkk48svv2ylS5dOBSvOAQEEEEAAAQQQQAABBBDIdwJ//fWXFS9ePO7zfv75561nz55xb8+GCCCAAAII5LVAjgRAAxe1YsUK++CDD0xFjr755htbsGCBFSlSxI34VMGjXr16WevWrQPd3c81a9bYuHHjrHHjxnbsscdaRkZGlvW8QQABBBBAAAEEEEAAAQQQ8C6we/du9z0s3sl/w4YN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zEoEEEAAAQQQQAABBBBAAAEEEEAAAQRSQIARoCnwIeTEKSig+Vb5mVZ4wb6p8P9tmu7uprz/933hhSXtrbIzCH4GgPiJAAIIIIAAAggggAACCCCAAAIIIJB2AowATbuPNOsF/b3nb7vts+tsUrF3bGu1TW5l6bXl7eS/zrKH/vWkFS1UNOsGvEMAAQQQQAABBBBAAAEEEEAAAQQQQCCNBAiAptGHyaUggAACCCCAAAIIIIAAAggggAACCCCAQFYBpsBn9eAdAggggAACCCCAAAIIIIAAAggggAACCKSRAAHQNPowuRQEEEAAAQQQQAABBBBAAAEEEEAAAQQQyCpAADSrB+8QQAABBBBAAAEEEEAAAQQQQAABBBBAII0ECICm0YfJpSCAAAIIIIAAAggggAACCCCAAAIIIIBAVgECoFk9eIcAAggggAACCCCAAAIIIIAAAggggAACaSRAADSNPkwuBQEEEEAAAQQQQAABBBBAAAEEEEAAAQSyChAAzerBOwQQQAABBBBAAAEEEEAAAQQQQAABBBBIIwECoGn0YXIpCCCAAAIIIIAAAggggAACCCCAAAIIIJBVgABoVg/eIYAAAggggAACCCCAAAIIIIAAAggggEAaCRAATaMPk0tBAAEEEEAAAQQQQAABBBBAAAEEEEAAgawCBECzevAOAQQQQAABBBBAAAEEEEAAAQQQQAABBNJIgABoGn2YXAoCCCCAAAIIIIAAAggggAACCCCAAAIIZBUgAJrVg3cIIIAAAggggAACCCCAAAIIIIAAAgggkEYCBEDT6MPkUhBAAAEEEEAAAQQQQAABBBBAAAEEEEAgq0DhrG8L1rvx48fbyy+/bH369LFGjRr5vvi//vrL3nzzTZs9e7Zt3LjR6tevb02bNrVTTz3VChUqFHF/3333nb3xxhv2888/W6lSpaxx48bWvn17q1u3bsRtWIEAAggggAACCCCAAAIIIIAAAggggAAC/gUy9u5r/jfL/1vMnz/frrvuOvv7779t0KBBLnDp56o2bdpkvXr1spUrV7rNKlasaBs2bHCv27RpY/fee68VLVo02y4VMB04cKBbXrp0adu5c6f7U6JECevXr58deeSR2bZhAQIIIIAAAggggAACCCCAAAIIIIAAAgjEJ1Agp8B/++23duedd7rgZ3xsZg888IALfrZs2dLef/99e/fd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OzF6nxMqUqUyCBAgQIAAAQIECJQvIAAt39wZCRAgQIAAAQIECBAgQIAAAQIECBAoScAj8CVBOw0BAgQIECBAgAABAgQIECBAgAABAuULCEDLN3dGAgQIECBAgAABAgQIECBAgAABAgRKEhCAlgTtNAQIECBAgAABAgQIECBAgAABAgQIlC8gAC3f3BkJECBAgAABAgQIECBAgAABAgQIEChJQABaErTTECBAgAABAgQIECBAgAABAgQIECBQvoAAtHxzZyRAgAABAgQIECBAgAABAgQIECBAoCQBAWhJ0E5DgAABAgQIECBAgAABAgQIECBAgED5AgLQ8s2dkQABAgQIECBAgAABAgQIECBAgACBkgQEoCVBOw0BAgQIECBAgAABAgQIECBAgAABAuULCEDLN3dGAgQIECBAgAABAgQIECBAgAABAgRKEhCAlgTtNAQIECBAgAABAgQIECBAgAABAgQIlC8gAC3f3BkJECBAgAABAgQIECBAgAABAgQIEChJQABaErTTECBAgAABAgQIECBAgAABAgQIECBQvoAAtHxzZyRAgAABAgQIECBAgAABAgQIECBAoCQBAWhJ0E5DgAABAgQIECBAgAABAgQIECBAgED5AgLQ8s2dkQABAgQIECBAgAABAgQIECBAgACBkgQEoCVBOw0BAgQIECBAgAABAgQIECBAgAABAuULCEDLN3dGAgQIECBAgAABAgQIECBAgAABAgRKEhCAlgTtNAQIECBAgAABAgQIECBAgAABAgQIlC8gAC3f3BkJECBAgAABAgQIECBAgAABAgQIEChJQABaErTTECBAgAABAgQIECBAgAABAgQIECBQvoAAtHxzZyRAgAABAgQIECBAgAABAgQIECBAoCQBAWhJ0E5DgAABAgQIECBAgAABAgQIECBAgED5AgLQ8s2dkQABAgQIECBAgAABAgQIECBAgACBkgQEoCVBOw0BAgQIECBAgAABAgQIECBAgAABAuULCEDLN3dGAgQIECBAgAABAgQIECBAgAABAgRKEhCAlgTtNAQIECBAgAABAgQIECBAgAABAgQIlC8gAC3f3BkJECBAgAABAgQIECBAgAABAgQIEChJQABaErTTECBAgAABAgQIECBAgAABAgQIECBQvoAAtHxzZyRAgAABAgQIECBAgAABAgQIECBAoCQBAWhJ0E5DgAABAgQIECBAgAABAgQIECBAgED5AgLQ8s2dkQABAgQIECBAgAABAgQIECBAgACBkgQEoCVBOw0BAgQIECBAgAABAgQIECBAgAABAuULCEDLN3dGAgQIECBAgAABAgQIECBAgAABAgRKEhCAlgTtNAQIECBAgAABAgQIECBAgAABAgQIlC8gAC3f3BkJECBAgAABAgQIECBAgAABAgQIEChJQABaErTTECBAgAABAgQIECBAgAABAgQIECBQvoAAtHxzZyRAgAABAgQIECBAgAABAgQIECBAoCQBAWhJ0E5DgAABAgQIECBAgAABAgQIECBAgED5AgLQ8s2dkQABAgQIECBAgAABAgQIECBAgACBkgQEoCVBOw0BAgQIECBAgAABAgQIECBAgAABAuULCEDLN3dGAgQIECBAgAABAgQIECBAgAABAgRKEhCAlgTtNAQIECBAgAABAgQIECBAgAABAgQIlC8gAC3f3BkJECBAgAABAgQIECBAgAABAgQIEChJQABaErTTECBAgAABAgQIECBAgAABAgQIECBQvkB9+ad0xq0lsG7duli6dOnWGt64BAgQIECAAAECBAgQIECAAAECBLaJQK9evaJLl827l1MAuk2mbOucNAWgS5Ys2TqDG5UAAQIECBAgQIAAAQIECBAgQIDANhLo0aPHZgegdetDs3XbqG6nJUCAAAECBAgQIECAAAECBAgQIECAwFYV2Lz7RrdqSQYnQIAAAQIECBAgQIAAAQIECBAgQIBAMQIC0GIcjUKAAAECBAgQIECAAAECBAgQIECAQAUFBKAVnBQlESBAgAABAgQIECBAgAABAgQIECBQjIAAtBhHoxAgQIAAAQIECBAgQIAAAQIECBAgUEEBAWgFJ0VJBAgQIECAAAECBAgQIECAAAECBAgUIyAALcbRKAQIECBAgAABAgQIECBAgAABAgQIVFBAAFrBSVESAQIECBAgQIAAAQIECBAgQIAAAQLFCAhAi3E0CgECBAgQIECAAAECBAgQIECAAAECFRQQgFZwUpREgAABAgQIECBAgAABAgQIECBAgEAxAgLQYhyNQoAAAQIECBAgQIAAAQIECBAgQIBABQUEoBWcFCURIECAAAECBAgQIECAAAECBAgQIFCMgAC0GEejECBAgAABAgQIECBAgAABAgQIECBQQQEBaAUnRUkECBAgQIAAAQIECBAgQIAAAQIECBQjIAAtxtEoBAgQIECAAAECBAgQIECAAAECBAhUUEAAWsFJURIBAgQIECBAgAABAgQIECBAgAABAsUICECLcTQKAQIECBAgQIAAAQIECBAgQIAAAQIVFBCAVnBSlESAAAECBAgQIECAAAECBAgQIECAQDECAtBiHI1CgAABAgQIECBAgAABAgQIECBAgEAFBQSgFZwUJREgQIAAAQIECBAgQIAAAQIECBAgUIyAALQYR6MQIECAAAECBAgQIECAAAECBAgQIFBBAQFoBSdFSQQIECBAgAABAgQIECBAgAABAgQIFCMgAC3G0SgECBAgQIAAAQIECBAgQIAAAQIECFRQQABawUlREgECBAgQIECAAAECBAgQIECAAAECxQgIQItxNAoBAgQIECBAgAABAgQIECBAgAABAhUUEIBWcFKURIAAAQIECBAgQIAAAQIECBAgQIBAMQIC0GIcjUKAAAECBAgQIECAAAECBAgQIECAQAUFBKAVnBQlESBAgAABAgQIECBAgAABAgQIECBQjIAAtBhHoxAgQIAAAQIECBAgQIAAAQIECBAgUEEBAWgFJ0VJBAgQIECAAAECBAgQIECAAAECBAgUIyAALcbRKAQIECBAgAABAgQIECBAgAABAgQIVFBAAFrBSVESAQIECBAgQIAAAQIECBAgQIAAAQLFCAhAi3E0CgECBAgQIECAAAECBAgQIECAAAECFRQQgFZwUpREgAABAgQIECBAgAABAgQIECBAgEAxAgLQYhyNQoAAAQIECBAgQIAAAQIECBAgQIBABQUEoBWcFCURIECAAAECBAgQIECAAAECBAgQIFCMgAC0GEejECBAgAABAgQIECBAgAABAgQIECBQQQEBaAUnRUkECBAgQIAAAQIECBAgQIAAAQIECBQjIAAtxtEoBAgQIECAAAECBAgQIECAAAECBAhUUEAAWsFJURIBAgQIECBAgAABAgQIECBAgAABAsUICECLcTQKAQIECBAgQIAAAQIECBAgQIAAAQIVFBCAVnBSlESAAAECBAgQIECAAAECBAgQIECAQDECAtBiHI1CgAABAgQIECBAgAABAgQIECBAgEAFBQSgFZwUJREgQIAAAQIECBAgQIAAAQIECBAgUIyAALQYR6MQIECAAAECBAgQIECAAAECBAgQIFBBAQFoBSdFSQQIECBAgAABAgQIECBAgAABAgQIFCMgAC3G0SgECBAgQIAAAQIECBAgQIAAAQIECFRQQABawUlREgECBAgQIECAAAECBAgQIECAAAECxQgIQItxNAoBAgQIECBAgAABAgQIECBAgAABAhUUEIBWcFKURIAAAQIECBAgQIAAAQIECBAgQIBAMQIC0GIcjUKAAAECBAgQIECAAAECBAgQIECAQAUFBKAVnBQlESBAgAABAgQIECBAgAABAgQIECBQjIAAtBhHoxAgQIAAAQIECBAgQIAAAQIECBAgUEEBAWgFJ0VJBAgQIECAAAECBAgQIECAAAECBAgUIyAALcbRKAQIECBAgAABAgQIECBAgAABAgQIVFBAAFrBSVESAQIECBAgQIAAAQIECBAgQIAAAQLFCAhAi3E0CgECBAgQIECAAAECBAgQIECAAAECFRQQgFZwUpREgAABAgQIECBAgAABAgQIECBAgEAxAgLQYhyNQoAAAQIECBAgQIAAAQIECBAgQIBABQUEoBWcFCURIECAAAECBAgQIECAAAECBAgQIFCMgAC0GEejECBAgAABAgQIECBAgAABAgQIECBQQQEBaAUnRUkECBAgQIAAAQIECBAgQIAAAQIECBQjIAAtxtEoBAgQIECAAAECBAgQIECAAAECBAhUUEAAWsFJURIBAgQIECBAgAABAgQIECBAgAABAsUICECLcTQKAQIECBAgQIAAAQIECBAgQIAAAQIVFBCAVnBSlESAAAECBAgQIECAAAECBAgQIECAQDECAtBiHI1CgAABAgQIECBAgAABAgQIECBAgEAFBQSgFZwUJREgQIAAAQIECBAgQIAAAQIECBAgUIyAALQYR6MQIECAAAECBAgQIECAAAECBAgQIFBBAQFoBSdFSQQIECBAgAABAgQIECBAgAABAgQIFCMgAC3G0SgECBAgQIAAAQIECBAgQIAAAQIECFRQQABawUlREgECBAgQIECAAAECBAgQIECAAAECxQgIQItxNAoBAgQIECBAgAABAgQIECBAgAABAhUUEIBWcFKURIAAAQIECBAgQIAAAQIECBAgQIBAMQIC0GIcjUKAAAECBAgQIECAAAECBAgQIECAQAUFBKAVnBQlESBAgAABAgQIECBAgAABAgQIECBQjIAAtBhHoxAgQIAAAQIECBAgQIAAAQIE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kydPjlNPPTXaet3UT37ykzjqqKPaHMtOAgQIECBAoGMKdIg7QPNMzZNPPpnd7Zn+5jc9/r654Wc612c/+9nslC2FrmmRpHSXadpSMNuwDRo0KIYPHx4ffPBBPPXUUw3NjT8ffvjhWLNmTdand+/eje0+ECBAgAABAgQIECDQtsCXvvSleOyxx+LII49s1vFjH/tY/O53v4vzzjuv2T4NBAgQIECAQOcQ6BQBaFpc6Kqrrop169bF+PHj4/jjj881u2mBpAceeCAefPDBJv3TXZ977bVXzJ49O/tlasOdt99+e6TV5ffcc8847LDDNtwVX/3qV7Pvv/jFL7JH8Bt2zps3L+64447s6ymnnNLQ7CcBAgQIECBAgAABAjkFRowYEY8++mjMnTs30pNfv/3tb2PWrFkxY8aMGD16dM5RdCNAgAABAgQ6okDd+lBwXUe8sA2vKYWSN9xwQ9bUtWvXDXc1+5wWRGr4m+P77rsvrrjiikjHTJ06tUnfRx55JC6++OLsrs1jjjkmW00y/XL1xBNPZI+3p0ds/uqv/qrJMekOz7POOiteeOGFSHeEpuNWr16dPb6fQtMjjjgivv/977e4an2TgXwhQIAAAQIECBAgQIAAAQIECBAgQCCXQId4B2h7Vzp9+vTGLimEbGtLj7Dn2dIKkumu0hRYpsfX05+0pTtDL7jggmbhZ9qXgtRrrrkmO+7++++PFMw2tH/lK1+JCRMmCD8zEf8gQIAAAQIECBAgQIAAAQIECBAgUIxAp7gDtBiq1kdJd2+mxY3SQke77rprrhAz3fn5pz/9KXssf4899ohNXYm+9WrsIUCAAAECBAgQIECAAAECBAgQIECgQUAA2iDhJwECBAgQIECAAAECBAgQIECAAAECHU6gUyyC1OFmzQURIECAAAECBAgQIECAAAECBAgQIJBLQACai0knAgQIECBAgAABAgQIECBAgAABAgRqUUAAWouzpmYCBAgQIECAAAECBAgQIECAAAECBHIJCEBzMelEgAABAgQIECBAgAABAgQIECBAgEAtCghAa3HW1EyAAAECBAgQIECAAAECBAgQIECAQC4BAWguJp0IECBAgAABAgQIECBAgAABAgQIEKhFAQFoLc6amgkQIECAAAECBAgQIECAAAECBAgQyCUgAM3FpBMBAgQIECBAgAABAgQIECBAgAABArUoIACtxVlTMwECBAgQIECAAAECBAgQIECAAAECuQQEoLmYdCJAgAABAgQIECBAgAABAgQIECBAoBYFBKC1OGtqJkCAAAECBAgQIECAAAECBAgQIEAgl4AANBeTTgQIECBAgAABAgQIECBAgAABAgQI1KKAALQWZ03NBAgQIECAAAECBAgQIECAAAECBAjkEhCA5mLSiQABAgQIECBAgAABAgQIECBAgACBWhQQgNbirKmZAAECBAgQIECAAAECBAgQIECAAIFcAgLQXEw6ESBAgAABAgQIECBAgAABAgQIECBQiwIC0FqcNTUTIECAAAECBAgQIECAAAECBAgQIJBLQACai0knAgQIECBAgAABAgQIECBAgAABAgRqUUAAWouzpmYCBAgQIECAAAECBAgQIECAAAECBHIJCEBzMelEgAABAgQIECBAgAABAgQIECBAgEAtCghAa3HW1EyAAAECBAgQIECAAAECBAgQIECAQC4BAWguJp0IECBAgAABAgQIECBAgAABAgQIEKhFAQFoLc6amgkQIECAAAECBAgQIECAAAECBAgQyCUgAM3FpBMBAgQIECBAgAABAgQIECBAgAABArUoIACtxVlTMwECBAgQIECAAAECBAgQIECAAAECuQQEoLmYdCJAgAABAgQIECBAgAABAgQIECBAoBYFBKC1OGtqJkCAAAECBAgQIECAAAECBAgQIEAgl4AANBeTTgQIECBAgAABAgQIECBAgAABAgQI1KKAALQWZ03NBAgQIECAAAECBAgQIECAAAECBAjkEhCA5mLSiQABAgQIECBAgAABAgQIECBAgACBWhQQgNbirKmZAAECBAgQIECAAAECBAgQIECAAIFcAgLQXEw6ESBAgAABAgQIECBAgAABAgQIECBQiwIC0FqcNTUTIECAAAECBAgQIECAAAECBAgQIJBLQACai0knAgQIECBAgAABAgQIECBAgAABAgRqUUAAWouzpmYCBAgQIECAAAECBAgQIECAAAECBHIJCEBzMelEgAABAgQIECBAgAABAgQIECBAgEAtCghAa3HW1EyAAAECBAgQIECAAAECBAgQIECAQC4BAWguJp0IECBAgAABAgQIECBAgAABAgQIEKhFAQFoLc6amgkQIECAAAECBAgQIECAAAECBAgQyCUgAM3FpBMBAgQIECBAgAABAgQIECBAgAABArUoIACtxVlTMwECBAgQIECAAAECBAgQIECAAAECuQQEoLmYdCJAgAABAgQIECBAgAABAgQIECBAoBYFBKC1OGtqJkCAAAECBAgQIECAAAECBAgQIEAgl4AANBeTTgQIECBAgAABAgQIECBAgAABAgQI1KKAALQWZ03NBAgQIECAAAECBAgQIECAAAECBAjkEhCA5mLSiQABAgQIECBAgAABAgQIECBAgACBWhQQgNbirKmZAAECBAgQIECAAAECBAgQIECAAIFcAgLQXEw6ESBAgAABAgQIECBAgAABAgQIECBQiwIC0FqcNTUTIECAAAECBAgQIECAAAECBAgQIJBLQACai0knAgQIECBAgAABAgQIECBAgAABAgRqUUAAWouzpmYCBAgQIECAAAECBAgQIECAAAECBHIJCEBzMelEgAABAgQIECBAgAABAgQIECBAgEAtCghAa3HW1EyAAAECBAgQIECAAAECBAgQIECAQC4BAWguJp0IECBAgAABAgQIECBAgAABAgQIEKhFAQFoLc6amgkQIECAAAECBAgQIECAAAECBAgQyCUgAM3FpBMBAgQIECBAgAABAgQIECBAgAABArUoIACtxVlTMwECBAgQIECAAAECBAgQIECAAAECuQQEoLmYdCJAgAABAgQIECBAgAABAgQIECBAoBYFBKC1OGtqJkCAAAECBAgQIECAAAECBAgQIEAgl4AANBeTTgQIECBAgAABAgQIECBAgAABAgQI1KKAALQWZ03NBAgQIECAAAECBAgQIECAAAECBAjkEhCA5mLSiQABAgQIECBAgAABAgQIECBAgACBWhQQgNbirKmZAAECBAgQIECAAAECBAgQIECAAIFcAgLQXEw6ESBAgAABAgQIECBAgAABAgQIECBQiwIC0FqcNTUTIECAAAECBAgQIECAAAECBAgQIJBLQACai0knAgQIECBAgAABAgQIECBAgAABAgRqUUAAWouzpmYCBAgQIECAAAECBAgQIECAAAECBHIJCEBzMelEgAABAgQIECBAgAABAgQIECBAgEAtCghAa3HW1EyAAAECBAgQIECAAAECBAgQIECAQC4BAWguJp0IECBAgAABAgQIECBAgAABAgQIEKhFAQFoLc6amgkQIECAAAECBAgQIECAAAECBAgQyCUgAM3FpBMBAgQIECBAgAABAgQIECBAgAABArUoIACtxVlTMwECBAgQIECAAAECBAgQIECAAAECuQQEoLmYdCJAgAABAgQIECBAgAABAgQIECBAoBYFBKC1OGtqJkCAAAECBAgQIECAAAECBAgQIEAgl4AANBeTTgQIECBAgAABAgQIECBAgAABAgQI1KKAALQWZ03NBAgQIECAAAECBAgQIECAAAECBAjkEqjP1UunmhBYt25dLFmypCZqVSQBAgQIEOjsAh988EEsXLgwevbsGQMGDOjsHK6fAAECBAgQIECAQJsCffv2ja5du7bZp7WddetDs3Wt7dReWwJpKteuXVtbRauWAAECBAh0MoHXXnstbrjhhnjiiScar3z//fePcePGxSGHHNLY5gMBAgQIECBAgAABAn8R6NKlS9TV1f2lYRM+CUA3AUtXAgQIECBAgMCWCMycOTPOOeecWL16dbNh0t9m/+M//mMcd9xxzfZpIECAAAECBAgQIEBg8wUEoJtv50gCBAgQIECAQG6B+fPnx/jx42PBggVx+OGHx5gxY7JH36dPnx433nhj1p7+RvvWW2+NIUOG5B5XRwIECBAgQIAAAQIE2hawCFLbPvYSIECAAAECBAoRSMFmCj+PPvro+OEPfxgHHnhgDBo0KLvjMwWg/fv3j/Q6m9tuu62Q8xmEAAECBAgQIECAAIH/FRCA+jeBAAECBAgQILCVBZYtWxb3339/pPcWTZw4sdm7i3bZZZc488wzsyoefPDBePvtt7dyRYYnQIAAAQIECBAg0HkEBKCdZ65dKQECBAgQILCNBJ5//vlYunRpDB8+PHbeeecWqxg9enT07t07uwt0wwWSWuyskQABAgQIECBAgACB3AIC0NxUOhIgQIAAAQIENk9g1apV2YFr1qxpdYBu3brFpz71qWz/s88+22K/NM71118f7733Xov7NRIgQIAAAQIECBAg0FxAANrcRAsBAgQIECBAoFCBvn37ZuO9/PLLMXfu3FbHPvTQQ7N9L7zwQot97rzzzpg8eXL87d/+bXZHaYudNBIgQIAAAQIECBAg0ERAANqEwxcCBAgQIECAQPECacGjwYMHZ4+3X3LJJfHWW2+1eJJhw4Zl7fPmzYuFCxc26ZPu+kwLKaXt2GOPzR6Xb9LBFwIECBAgQIAAAQIEWhQQgLbIopEAAQIECBAgUJxAXV1dTJgwIVsEKd0B+uijj7Y4+J577hk9evTI9r300ktN+qSV4tNiSv369YuxY8c22ecLAQIECBAgQIAAAQKtCwhAW7exhwABAgQIECBQmMCoUaPi0ksvjZEjR8bf/M3ftDhu165dY5999sn2vfLKK419Uhg6ZcqU7Pv48eNju+22a9znAwECBAgQIECAAAECbQvUt73bXgIECBAgQIAAgaIEjjnmmDj66KMj3RHa2rbvvvvGzJkzI70vtGGbNGlS9nHo0KFxwgknNDT7SYAAAQIECBAgQIBADgF3gOZA0oUAAQIECBAgsKkCq1evjnfffTdmz54d6XPD1lb4mfo0vAd0xowZ2SEPPPBANHyeOHFipLtEbQQIECBAgAABAgQI5BdwB2h+Kz0JECBAgAABAu0KpLDz7rvvjltuuSUWLVqU9U+PrH/uc5+Ls88+u/Edn60NdPDBB2e75s+fH08//XTccMMN2fdPf/rT8YlPfKK1w7QTIECAAAECBAgQINCKQN269Vsr+zQTIECAAAECBAhsgkD6teqiiy6Kxx57rMWjDjjggPjhD3+YLWTUYof/azzjjDNizpw50adPn/jwww+je/fu8ctf/jIGDRrU1mH2ESBAgAABAgQIECDQgoBH4FtA0USAAAECBAgQ2ByBX/ziF1n4uf3220d6XP1Xv/pV3HTTTTF69OhsuPRuz6uuuqrdoRtWeU/hZ9pOPfVU4We7ajoQIECAAAECBAgQaFnAHaAtu2glQIAAAQIECGySQHpk/eSTT476+vq49tprY7/99mtyfHqU/fbbb88WQEp3c+65555N9m/8ZcKECTFr1qzYaaedYvLkydG7d++Nu/hOgAABAgQIECBAgEAOAXeA5kDShQABAgQIECDQnsC9994ba9eujeOPP75Z+JmOHTduXHYXZ3pMPgWa7W3f/va3Y+DAgfGNb3xD+Nkelv0ECBAgQIAAAQIE2hAQgLaBYxcBAgQIECBAIK/Af/3Xf2VdjzjiiBYPSe/xPPHEE7N9TzzxRLT3GvYhQ4bEzTffHMcee2yL42kkQIAAAQIECBAgQCCfgAA0n5NeBAgQIECAAIE2BVatWpXtX7NmTav9jjrqqGzfwoULY+7cuS32mzp1aqQ/aevXr1/2yHz2xT8IECBAgAABAgQIENgsAQHoZrE5iAABAgQIECDQVKBv375Zw/333990xwbfBg8eHLvuumvW8sILL2yw538/pkWPfvzjH8c///M/x2233dZsvwYCBAgQIECAAAECBDZdQAC66WaOIECAAAECBAg0EzjhhBOytgcffDBb7KhZh/9rGDZsWPbppZdeatbl1ltvjXR3aM+ePeO4445rtl8DAQIECBAgQIAAAQKbLiAA3XQzRxAgQIAAAQIEmgkcc8wxsf/++2ft//Zv/xZLly5t1ic1tBaAvvHGG3HXXXdlx5x++ukxYMCAFo/XSIAAAQIECBAgQIDApgkIQDfNS28CBAgQIECAQIsC9fX1cdVVV8VBBx2UPcLeu3fvFvs1BKB/+tOfslXjGzpde+21sXr16uwR+dNOO62h2U8CBAgQIECAAAECBLZQoH4Lj3c4AQIECBAgQIDA/wmk0PPqq69uc+GifffdN+u9YsWK+POf/xx77bVXPPPMM/H4449n7WeffXb06NGDKQECBAgQIECAAAECBQm4A7QgSMMQIECAAAECBNJK8L/5zW9i3bp1rWL0798/Bg4cmO2fMWNGdtdnCk3Tlu4eTY/S2wgQIECAAAECBAgQKE7AHaDFWRqJAAECBAgQ6MQCy5Yti4suuiimTZuWLWR05plntqpx8MEHx3333RcPPfRQrFy5Ml599dXsrtGJEye2eowdBAgQIECAAAECBAhsnoA7QDfPzVEECBAgQIAAgSYC99xzTxZ+psYUhra1HXnkkdnuFJbedNNN2ee0ivw+++zT1mH2ESBAgAABAgQIECCwGQIC0M1AcwgBAgQIECBAYGOBOXPmZE2HHXZYtHcn51FHHRVDhw7N+qfV4vv06RPjx4/feEjfCRAgQIAAAQIECBAoQEAAWgCiIQgQIECAAIHOLfD+++/Hxz72sQzhgAMOaBejS5cucdZZZzX2Gzt2bOywww6N330gQIAAAQIECBAgQKA4AQFocZZGIkCAAAECBDqhQAo/v/a1r2WPv/fs2TOefPLJXArpTtHTTjst9thjjzj55JNzHaMTAQIECBAgQIAAAQKbLlC3fpXS1pcp3fTxHEGAAAECBAgQ6FQCaTGjK664osk1f/WrX40JEyZE165dm7S39GXRokWRVoa3ESBAgAABAgQIECCwdQQEoFvH1agECBAgQIBAJxKYMWNGTJo0KV566aXGq04rvV966aUebW8U8YEAAQIECBAgQIDAthEQgG4bd2clQIAAAQIEOphAeqhmypQp2aru7733XnZ1ffv2jfR+zy9/+ctRX1/fwa7Y5RAgQIAAAQIECBCoDQEBaG3MkyoJECBAgACBGhFIq7rfeuutceedd8bq1auzqgcPHhznnntujBw5skauQpkECBAgQIAAAQIEOo6AALTjzKUrIUCAAAECBCok8Oabb8a1114bjz32WGNVhx9+eJxzzjmRAlEbAQIECBAgQIAAAQLlCAhAy3F2FgIECBAgQKCDCUybNi1b+X327NkxbNiwOOGEE2K33XZrdpWpX3o/6KuvvprtSwsj/eAHP3A3aDMpDQQIECBAgAABAgS2joAAdOu4GpUAAQIECBDowAJ33XVXXH311U2uMAWbX/jCF2L8+PGx4447NtmXHoW/++674+c//3l079497rjjjujVq1eTPr4QIECAAAECBAgQILB1BASgW8fVqAQIECBAgEAHFXjggQfiu9/9bnZ1H/nIR6Jfv37x8ssvx6pVq7K23r17xxlnnBGnnHJKdOvWrYnCkiVL4vXXX48DDjigSbsvBAgQIECAAAECBAhsPQEB6NazNTIBAgQIECBQwwJr166NLl26NLuC0047LdL7PdPq7ulPXV1dLFy4MNJdob/+9a9j5cqV2TG77757tvDREUcc0WwMDQQIECBAgAABAgQIlCcgAC3P2pkIECBAgACBGhGYOnVqtpL7xIkT46CDDmqsOr3vc9y4cTF06NC45ZZbGtsbPrz99tvxL//yL5GOb9hGjBiRBaFDhgxpaPKTAAECBAgQIECAAIESBZrf1lDiyZ2KAAECBAgQIFA1gfQoewoxX3nllSy4vPjii+Odd97Jypw/f37285BDDmmx7LQI0uWXX569HzSFpGl75plnsjtF04rwNgIECBAgQIAAAQIEyhcQgJZv7owECBAgQIBAhQXSezvPPffcGDRoUFblww8/HGPGjImf/exnjYsbbbfddm1ewSc+8Ym4+eab48ILL4z+/fvHmjVrIj1SbyNAgAABAgQIECBAoHwBj8CXb+6MBAgQIECAQA0IpDtB77zzzuxR+GXLlmUVpzBz0aJFse+++2Yruue5jPfffz9uv/32LERtLzjNM54+BAgQIECAAAECBAhsmoAAdNO89CZAgAABAgQ6mcCCBQvipptuiilTpjS58nPOOSdOPfXUJm2+ECBAgAABAgQIECBQPQEBaPXmREUECBAgQIBABQVefPHFmDRpUvzxj39srO7oo4+O888/P3baaafGNh8IECBAgAABAgQIEKiWgAC0WvOhGgIECBAgQKDiAg888EBcf/318T//8z9Zpb169YozzjgjTjnllOjevXvFq1ceAQIECBAgQIAAgc4nIADtfHPuigkQIECAAIGcAqtXr465c+fGG2+8EWmF9yFDhkSPHj1i+fLl2Xs9J0+eHCtXrsxGS4sm/f3f/32MGjUq5+i6ESBAgAABAgQIECBQhoAAtAxl5yBAgAABAgRqTuDZZ5+NK6+8Mt58883G2n/84x/HiBEjGr+/++67cd1110VaKb5h++QnPxnnnXde7L333g1NfhIgQIAAAQIECBAgsA0FBKDbEN+pCRAgQIAAgWoKpPDzggsuiLVr1zYWuMsuu2Srwnfp0qWxreHD9OnTs/eDzp49O2tKfb7//e/HEUcc0dDFTwIECBAgQIAAAQIEtpGAAHQbwTstAQIECBAgUE2BxYsXx+mnnx6LFi2Kww8/PL71rW9FXV1d1NfXR//+/ZsU/dxzz8X2228fQ4cOzcLSe++9N1sxvmvXrnHHHXdEej+ojQABAgQIECBAgACBbSsgAN22/s5OgAABAgQIVEzgV7/6Vdx4443xqU99KruLM4WZLW3z58+PsWPHRtqfFkVK7whN24cffhivv/56DB8+vKXDtBEgQIAAAQIECBAgULJA82e4Si7A6QgQIECAAAECVRKYOnVqVs5JJ52UhZut1ZbuEF2zZk0sWLAgu+uzoV+fPn2Enw0YfhIgQIAAAQIECBCogIAAtAKToAQCBAgQIECgOgLp7s20DRw4sM2i9tlnn7j44ouzPr///e9jyZIlbfa3kwABAgQIECBAgACBbSMgAN027s5KgAABAgQIVFSg4b2dr776arsVjhw5Mnv0fdWqVTFv3rx2++tAgAABAgQIECBAgED5AgLQ8s2dkQABAgQIEKiwwEc/+tGsuvQe0OXLl7dbaUNg2m5HHQgQIECAAAECBAgQ2CYCAtBtwu6kBAgQIECAQFUF0grwaXvrrbfin/7pn2LhwoWtlvrOO+/E3Llzs9Xh99prr1b72UGAAAECBAgQIECAwLYTEIBuO3tnJkCAAAECBCoocOCBB8bXvva1rLJp06bF3/3d38WDDz6YLXi0YbkrVqyIyy67LNatWxejR4+O+vr6DXf7TIAAAQIECBAgQIBARQTq1v/Svq4itSiDAAECBAgQIFAZgeuvvz4mT57cWE9aFCkFnennBx98EPfcc0+8+eabkR6Zv+6668Kj8I1UPhAgQIAAAQIECBColIAAtFLToRgCBAgQIECgSgLpDtBJkyZFawsi7bbbbln4OWDAgCqVrRYCBAgQIECAAAECBDYQEIBugOEjAQIECBAgQKAlgenTp8fdd98dr7zySvbIe+/evePEE0+Mz3/+89GtW7eWDtFGgAABAgQIECBAgEBFBASgFZkIZRAgQIAAAQIECBAgQIAAAQIECBAgULyARZCKNzUiAQIECBAgQIAAAQIECBAgQIAAAQIVERCAVmQilEGAAAECBAgQIECAAAECBAgQIECAQPECAtDiTY1IgAABAgQIECBAgAABAgQIECBAgEBFBASgFZkIZRAgQIAAAQIECBAgQIAAAQIECBAgULyAALR4UyMSIECAAAECBAgQIECAAAECBAgQIFARAQFoRSZCGQQIECBAgAABAgQIECBAgAABAgQIFC8gAC3e1IgECBAgQIAAAQIECBAgQIAAAQIECFREQABakYlQBgECBAgQIECAAAECBAgQIECAAAECxQsIQIs3NSIBAgQIECBAgAABAgQIECBAgAABAhUREIBWZCKUQYAAAQIECBAgQIAAAQIECBAgQIBA8QIC0OJNjUiAAAECBAgQIECAAAECBAgQIECAQEUEBKAVmQhlECBAgAABAgQIECBAgAABAgQIECBQvIAAtHhTIxIgQIAAAQIECBAgQIAAAQIECBAgUBEBAWhFJkIZBAgQIECAAAECBAgQIECAAAECBAgULyAALd7UiAQIECBAgAABAgQIECBAgAABAgQIVERAAFqRiVAGAQIECBAgQIAAAQIECBAgQIAAAQLFCwhAizc1IgECBAgQIECAAAECBAgQIECAAAECFREQgFZkIpRBgAABAgQIECBAgAABAgQIECBAgEDxAgLQ4k2NSIAAAQIECBAgQIAAAQIECBAgQIBARQQEoBWZCGUQIECAAAECBAgQIECAAAECBAgQIFC8gAC0eFMjEiBAgAABAgQIECBAgAABAgQIECBQEQEBaEUmQhkECBAgQIAAAQIECBAgQIAAAQIECBQvIAAt3tSIBAgQIECAAAECBAgQIECAAAECBAhUREAAWpGJUAYBAgQIECBAgAABAgQIECBAgAABAsULCECLNzUiAQIECBAgQIAAAQIECBAgQIAAAQIVERCAVmQilEGAAAECBAgQIECAAAECBAgQIECAQPECAtDiTY1IgAABAgQIECBAgAABAgQIECBAgEBFBASgFZkIZRAgQIAAAQIECBAgQIAAAQIECBAgULyAALR4UyMSIECAAAECBAgQIECAAAECBAgQIFARAQFoRSZCGQQIECBAgAABAgQIECBAgAABAgQIFC8gAC3e1IgECBAgQIAAAQIECBAgQIAAAQIECFREQABakYlQBgECBAgQIECAAAECBAgQIECAAAECxQsIQIs3NSIBAgQIECBAgAABAgQIECBAgAABAhUREIBWZCKUQYAAAQIECBAgQIAAAQIECBAgQIBA8QIC0OJNjUiAAAECBAgQIECAAAECBAgQIECAQEUEBKAVmQhlECBAgAABAgQIECBAgAABAgQIECBQvMD/B0qNyHs1z8slAAAAAElFTkSuQmCC){width="672"}
:::
:::
