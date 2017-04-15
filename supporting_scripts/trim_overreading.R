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
  
  # Construct all partial marker sequences used for end alignments
  markers <- sapply(0:(nchar(marker)-1), function(i){
    substr(marker, 1, nchar(marker) - i)
  })
  
  # Find all alignments of marker sequences in given reads
  alignments <- do.call(
    rbind,
    lapply(markers, function(marker, reads, percentID){
      mismatch <- round( nchar(marker) - percentID*nchar(marker) )
      vmp <- vmatchPattern(marker, reads, max.mismatch = mismatch)
      as.data.frame(unlist(vmp))
    }, reads = reads, percentID = percentID))
  
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
  trim_reads <- reads[seqnames(aln_ranges)]
  trim_reads <- DNAStringSet(trim_reads, start = 1, end = start(aln_ranges)-1)
  
  # Gather all trimmed and untrimmed reads and reorder to match input
  trimmed_reads <- c(trim_reads, reads[!names(reads) %in% names(trim_reads)])
  trimmed_reads[names(reads)]
}