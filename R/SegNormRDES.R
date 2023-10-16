#' SegNormRDES: Performs segmentation on a given data frame using the specified segmentation method.
#'
#' @param df A data frame to be segmented.
#' @param id An identifier for the data.
#' @param seg.method A character string specifying the segmentation method to be used. Default is "slm" (Segmented Linear Model).
#'
#' @return A data table with the calculated values for each sequence.
#'
#' @examples
#' SegNormRDES(df, id, seg.method="slm")
#'
#' @export

SegNormRDES <- function(df,id,seg.method="slm"){
    ##EDIT: include SLM segmentation
    if (seg.method == "slm") {
        print("segment with SLM")
        df.ls <- split(df, df$seqnames)
        res <- lapply(df.ls, function(df) {
            logratio <- log2(df$ratio + 0.001)
            logratio[is.na(logratio)] <- 0
            slm <-
                HSLM(
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
        res <- rbindlist(res)
    }
}
