#' prepareBed
#' @description takes in a bed file and a reference genome, and returns a
#' GRanges object with a metadata column 'gcbias'.
#' The function first checks if the reference genome is a valid path to a fasta
#' file or a BSgenome object. If it is not, the function throws an error.
#' Next, the function reads in the bed file using the fread function from the
#' data.table package. It then orders the bed file by chromosome and start
#' position, and removes any duplicated capture designs.
#' Finally, the function annotates the GC ratio for each interval in the bed
#' file. If the reference genome is a BSgenome object, the function uses the
#' getSeq function to extract the sequence for each interval, calculates the GC
#' content, and calculates the GC bias. If the reference genome is a fasta file,
#'  the function uses the scanFa function from the Rsamtools package to extract
#'  the sequence for each interval, calculates the GC content, and calculates
#'  the GC bias.
#' @param bedFile a path to capture design or exon file
#' @param ref.genome a path to indexed fa file or a BSgenome object
#' @importFrom GenomicRanges GRanges
#' @importFrom Biostrings getSeq
#' @importFrom Biostrings letterFrequency
#' @importFrom Rsamtools scanFa
#' @importFrom IRanges IRanges
#'
#' @return The function returns the GRanges object with the metadata column
#' gcbias', excluding any intervals where the GC bias is not available.
#' @export

prepareBed <- function(bedFile,ref.genome){
    if(class(ref.genome)!="BSgenome"&!is.character(ref.genome)){
        stop("Please provide a valid path to reference.fa or BSgenome object")}
    print("***Reading bed file***")
    bed <- fread(bedFile,header = F,stringsAsFactors = F)
    ord <- order(bed$V1, bed$V2)
    bedOrdered <- bed[ord,]

    toRemIdx <- which(duplicated(paste(bedOrdered$V1,"_", bedOrdered$V2, "_",bedOrdered$V3)))
    if(length(toRemIdx)!=0){
        print("***Removing duplicated capture design***")
        bedOrdered <- bedOrdered[-(toRemIdx),]
    }
    print("***Annotating GC ratio***")
    interval.gr <- GRanges(bedOrdered$V1,ranges = IRanges(bedOrdered$V2,bedOrdered$V3))
    if(class(ref.genome)=="BSgenome"){
        print("Calculating GC-content...")
        x <- getSeq(ref.genome,interval.gr)
        GC.count <- letterFrequency(x,"GC")
        all.count <- letterFrequency(x,"ATGC")
        interval.gr$gc_bias <- as.vector(ifelse(all.count==0,NA,GC.count/all.count))
        # exclude unavailable regions

    } else if (file.exists(ref.genome)){
        print("Calculating GC-content...")
        x <- scanFa(ref.genome, interval.gr)
        GC.count <- letterFrequency(x,"GC")
        all.count <- letterFrequency(x,"ATGC")
        interval.gr$gc_bias <- as.vector(ifelse(all.count==0,NA,GC.count/all.count))
        # exclude unavailable regions

    }
    interval.gr$mid <- (start(interval.gr)+end(interval.gr))/2
    return(interval.gr[which(!is.na(interval.gr$gc_bias))])

}
