
#' Z-TPM internal function
#' 
#' @keywords internal
#' @param group a vector of numeric cor value
#' @param tpmDtOrdered an output from reordertpmDt (tpmDtOrdered)
#'
#' @return ratio
#' @noRd
.calRatio <- function(group,tpmDtOrdered){
    min.cor <- 0.9
    min.nref <- 65
    group<-as.numeric(unlist(group))
    gsize <- length(which(group>min.cor))
    Idx <- sort(group, decreasing = TRUE, index.return = TRUE)$ix[1:min(gsize, min.nref)]
    ll<-tpmDtOrdered[,..Idx]
    l <- as.numeric(unlist(ll[,1]))
    Vmean <- rowMeans(as.matrix(ll),na.rm = T)
    ratio <- l/Vmean
    names(ratio) <- colnames(tpmDtOrdered)[Idx[1]]
    return(ratio)
}
