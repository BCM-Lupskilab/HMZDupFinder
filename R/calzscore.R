#' Z-TPM internal function
#' @description A function that calculate Z-TPM based on nearest references.
#'
#' @keywords internal
#' @param  group a vector of numeric cor value
#' @param tpmDtOrdered an output from reordertpmDt (tpmDtOrdered)
#'
#' @return Z-TPM
#' @noRd
#'

.calZscore <- function(group,tpmDtOrdered){
    min.cor <- 0.9
    min.nref <- 65
    group<-as.numeric(unlist(group))
    gsize <- length(which(group>min.cor))
    Idx <- sort(group, decreasing = TRUE, index.return = TRUE)$ix[1:min(gsize, min.nref)]
    ll <- tpmDtOrdered[,..Idx]
    l <- as.numeric(unlist(ll[,1]))
    Vmean <- rowMeans(as.matrix(ll),na.rm = T)
    Vsd <- rowSds(as.matrix(ll),na.rm = T)
    Zscore <- as.numeric((l-Vmean)/Vsd)
    names(Zscore) <- colnames(tpmDtOrdered)[Idx[1]]
    return(Zscore)
}


