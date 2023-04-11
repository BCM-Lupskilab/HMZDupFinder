##' Calculate TPMs from BAM files
##'
##' @description
##' `calcTPMsFromBAMs()` is used to calculate TPM value from BAM file.
##' It takes in a vector of paths to bed files (bedFile), a vector of paths 
##' to BAM files (bamFiles), a vector of sample identifiers (sampleNames), 
##' an output directory (outputDir), and a numeric value of cores (mc.cores).
##'  NOTE: Index files (.bai) should be available for each BAM file
##'
##'
##' @param bedOrdered   a ordered bed data table
##' tab separeated bed file containing capture design probes with the following 4 columns:
##'                    - chromosme
##'                    - start
##'                    - stop
##'                    - gene name *
##' @param bamFiles	  	  a vector of paths to bam files
##' @param sampleNames	  a vector of sample identifiers
##' @param outputDir	  a character value of output directory
##' @param mc.cores	      a numeric value of cores
##' @export

calcTPMsFromBAMs <- function (bedOrdered, bamFiles, sampleNames, outputDir, mc.cores=4){
    bed <- as.data.table(bedOrdered)
    df <- data.frame(cbind(1:nrow(bed), bed[,1:3], 1:nrow(bed)))
    colnames(df) <- c("GeneID", "Chr", "Start", "End", "Strand")
    if (!dir.exists(outputDir)){dir.create(outputDir)}
    out <- pbmclapply(1:length(bamFiles), function(i){
        file <- bamFiles[i]
        prefix <- sampleNames[i]
        print(prefix)
        outputFile <- paste0(outputDir, prefix, ".tpm.bed")
        res <-  featureCounts(files=file, annot.ext=df, allowMultiOverlap=TRUE, nthreads=3,isPairedEnd=TRUE)
        counts <- data.frame(bed[,1:3],V4=res$counts)
        colnames(counts)[4]<-"counts"
        RPK<-as.numeric(counts$count)*1e3/as.numeric(df$End - df$Start)
        scale<-sum(RPK)/1e6
        counts$TPM <-  RPK/scale
        setnames(counts,c("chr","start","end","counts","TPM"))
        fwrite(counts[,c("chr","start","end","counts","TPM")], file=outputFile, sep="\t", col.names=T, row.names=F, quote=F)
    },mc.cores=mc.cores)
}
