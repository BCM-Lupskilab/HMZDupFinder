


#' mapTxtoLine(bed.gr, refGene.gr)
#'
#' @description The function first finds the overlaps between the probes and
#' the reference transcripts using the findOverlaps function.
#' It then extracts the transcript IDs from the reference transcripts that
#' overlap with the probes, and calculates the number of lines in the probe file
#' that correspond to each transcript.
#'
#' @param bed.gr a GRanges object that represents the genomic coordinates of a
#' set of probes.
#' @param refGene.gr a GRanges object that represents the genomic coordinates
#' of a set of reference transcripts, metatable includes a transcript_id column
#' @importFrom GenomicRanges findOverlaps
#' @return a list where each element corresponds to a transcript, and contains a
#'  vector of line numbers in the probe file that correspond to that transcript.
#' @export
#'
mapTxtoLine <- function(bed.gr, refGene.gr) {
  ov <- findOverlaps(bed.gr, refGene.gr)
  transcript_ids <- refGene.gr$transcript_id[subjectHits(ov)]
  #exon_id <- refGene.gr$exon_number[subjectHits(ov)]
  tx.idx <- rle(transcript_ids)$lengths
  line_nums <- lapply(seq_along(tx.idx),function(i){
    c((sum(tx.idx[1:i])-tx.idx[i]+1),sum(tx.idx[1:i]))
  })
  line_num_bed <- lapply(line_nums,function(i){
    unique(c(queryHits(ov)[i[1]:i[2]]))
  })
  names(line_num_bed) <- rle(transcript_ids)$values
  return(line_num_bed)
}








