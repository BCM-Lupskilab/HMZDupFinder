##------------------------------------------------------------------------------
##' Calculate Z-TPM with nearest references and write bed file for visualizetion
##'
##'
##' @param perMat cor matrix of samples
##' @param tpmDtOrdered	 object returned by reorderBedAndRpkmDt
##' @export
##------------------------------------------------------------------------------

CalcuZtpm <- function(perMat,tpmDtOrdered, mc.cores = 4){
    print("[******Calculate Z-TPM**********]")
    candidateExon <- pbmclapply(perMat,.calZscore,tpmDtOrdered,mc.cores=mc.cores)
    candidateZscore<-bind_cols(candidateExon)
    return(candidateZscore)
}

