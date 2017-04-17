#' Trim 5' ends of reads matching the marker sequence.
#' 
#' \code{trim_overreading} removes 5' nucleotide sequences matching the marker.
#' 
#' @description 
#' 
#' @usage 
#' 
#' @param marker
#' 
#' @param reads
#' 
#' @param percentID
#' 
#' @example 
#' 
#' 
#' @author Christopher Nobles, Ph.D.

trim_leading <- function(marker, reads, percentID, filter, buffer = 1){
  require(IRanges)
  require(GenomicRanges)
  require(Biostrings)
  
  submat <- banmat()

  if(nchar(marker) > 0){
    # Find all reads with marker sequence
    alignments <- pairwiseAlignment(
      DNAStringSet( reads, start = 1, end = (nchar(marker) + buffer) ),
      as.character(marker), 
      type = "overlap", 
      substitutionMatrix = submat, 
      gapOpening = 3, 
      gapExtension = 1, 
      scoreOnly = FALSE)
  
    alignments <- alignments[alignments@score/nchar(marker) >= percentID]
  
    # Cut reads where marker aligned
    reads[names(alignments@pattern@unaligned)] <- DNAStringSet(
      reads[names(alignments@pattern@unaligned)],
      start = end(alignments@pattern@range)+1)
    if(filter) reads <- reads[names(alignments@pattern@unaligned)]
    return(reads)
  }else{
    return(reads)
  }
}