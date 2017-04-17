#' Trim 3' ends of reads when the sequence has overread the desired sequence.
#' 
#' \code{trim_overreading} removes 3' ends of nucleotide sequences matching 
#' the marker sequence.
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

trim_overreading <- function(marker, reads, percentID){
  require(IRanges)
  require(GenomicRanges)
  require(Biostrings)
  
  submat <- banmat()
  max_marker_length <- nchar(marker)
  
  # Construct all partial marker sequences used for end alignments
  markers <- sapply(0:(nchar(marker)-1), function(i){
    substr(marker, 1, nchar(marker) - i)
  })
  
  # Find all alignments of marker sequences in given reads
  alignments <- do.call(
    rbind,
    lapply(markers, function(marker, reads, percentID){
      pwa <- pairwiseAlignment(
        reads, 
        as.character(marker), 
        type = "overlap",
        substitutionMatrix = submat,
        gapOpening = 3, 
        gapExtension = 1, 
        scoreOnly = FALSE)
      pwa <- pwa[pwa@score/nchar(marker) >= percentID]
      data.frame(pwa@pattern@range, names = names(pwa@pattern@unaligned))}, 
    reads = reads, 
    percentID = percentID))
  
  # Filter partial alignments to only contain end of reads
  alignments$read_length <- width(reads)[
    match(alignments$names, names(reads))]
  alignments <- alignments[
    alignments$width == nchar(marker) |
      alignments$end == alignments$read_length,]
  
  # Determine first alignment from those remaining for each read
  aln_ranges <- reduce(GRanges(
    seqnames = alignments$names,
    ranges = IRanges(
      start = alignments$start,
      end = alignments$read_length),
    strand = rep("*", nrow(alignments))),
    ignore.strand = TRUE)
  
  # Cut reads where marker aligned
  reads[seqnames(aln_ranges)] <- DNAStringSet(
    reads[seqnames(aln_ranges)], start = 1, end = start(aln_ranges)-1)
  reads
}