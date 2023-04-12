
#' prepareTPMData
#' read individual TPM bed file and prepare matrix for z-TPM calculation.
#' @param tpmFile a vector of paths of tpm.bed file
#' @param chrom a vector of chromosome to read
#' @param mc.cores number of core for parallel running
#' @importFrom data.table setnames
#' @importFrom data.table setDT
#'
#' @return a data.table object with sample id named columns; each column stores
#' the TPM values
#' @export
#'
#'
prepareTPMData <- function(tpmFile, chrom=NULL, mc.cores=2){
    print ("Reading TPM files ...")
    tpmList <- pbmclapply(tpmFile, function(file){
        if (file.info(file)$size ==0) {return(NULL)}
        if(is.null(chrom)){
            return(fread(file, header = TRUE)$TPM)
        } else {
            ## need the tpmfile has a header name chr
           return(fread(file, header = TRUE)[chr %in% chrom]$TPM)
        }

    }, mc.cores=mc.cores)
    fids <- gsub(".tpm.bed*","",basename(tpmFile))
    names(tpmList) <- fids
    print ("Removing empty elements ...")
    tpmList2 <- Filter(function(x){!is.null(x)}, tpmList)
    print ("Creating matrix ...")
    tpmDf <- as.data.table(do.call(cbind,tpmList2))
    setnames(tpmDf, names(tpmList2))
    rm(tpmList);rm(tpmList2);gc();
    print ("Creating data.table (may take a while)...")
    setDT(tpmDf)
    tpmDf
}
