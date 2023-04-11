##------------------------------------------------------------------------------
##' Reorder BED and TPM matrix
##'
##' @
##' @param bedFile    tab separeated bed file containing capture design probes with the following 4 columns:
##'                    - chromosme
##'                    - start
##'                    - stop
##'                    - gene name
##' @param tpmDt	 object returned by prepareTPMData
##' @return a list of reordered bed and tmpDt
##' @export
##------------------------------------------------------------------------------
reordertpmDt <- function(bedFile, tpmDt){
    bed <- read.csv(bedFile, sep="\t", stringsAsFactors=F, header=F)
    ord <- order(bed$V1, bed$V2)
    bedOrdered <- bed[ord,]
    tpmDtOrdered <- tpmDt[ord,]
    toRemIdx <- which(duplicated(paste(bedOrdered$V1,"_", bedOrdered$V2, "_",bedOrdered$V3)))
    if (length(toRemIdx)>0){
        bedOrdered <- bedOrdered[-(toRemIdx),]
        tpmDtOrdered<- tpmDtOrdered[-(toRemIdx),]
    }
    colnames(tpmDtOrdered) <- colnames(tpmDt)
    list(bedOrdered=bedOrdered, tpmDtOrdered=tpmDtOrdered)
}
