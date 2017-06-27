#' Trim 3' ends of sequences when the sequence has overread into
#' synthetic sequence
#' \code{trim_overreading} removes 3' ends of nucleotide sequences completely or
#' partially matching the trimSequence.
#' @param seqs DNAStringSet of reads or unique sequences
#' @param trimSequence character string of lenth 1, such as "GAAAATC". This
#' string will be used to match the end of sequences, upon which the matching
#' portion will be trimmed from the start of the match to the end of the
#' sequence. Ambiguous nucleotides within the sequence can be used for alignment,
#' but matching sequences are not recored.
#' @param percentID numeric between 0 and 1.0 denoting the minimum percent
#' identity acceptable for an matching alignment.
#' @param maxSeqLength numeric/integer the maximum length to consider of the
#' trimSequence to use for alignments. Using the full length of sequence
#' avaliable for matching can many times be computationally intensive and
#' unnessesarily time consuming. Further, identical results can be obtained
#' using only a portion of the sequence. For example, setting the maxSeqLength
#' to 15L will only use the first 15 nucleotides of the trimSequence.
#' @author Christopher Nobles, Ph.D.

trim_overreading <- function(seqs, trimSequence,
                             percentID, maxSeqLength = NULL){
  require(IRanges)
  require(GenomicRanges)
  require(Biostrings)
  
  if(is.null(names(seqs))){
    noNames <- TRUE
    names(seqs) <- as.character(1:length(seqs))
  }else{
    noNames <- FALSE
  }
  
  if(!is.null(maxSeqLength)){
    trimSequence <- Biostrings::DNAStringSet(
      trimSequence, start = 1L, end = min(nchar(trimSequence), maxSeqLength))
  }
  
  trimSeqs <- sapply(0:(nchar(trimSequence)-1), function(i){
    substr(trimSequence, 1, nchar(trimSequence) - i)
  })
  
  alignments <- do.call(
    rbind,
    lapply(trimSeqs, function(trimSeq, seqs, percentID){
      mismatch <- round( nchar(trimSeq) - percentID*nchar(trimSeq) )
      vmp <- Biostrings::vmatchPattern(
        trimSeq, seqs, max.mismatch = mismatch, fixed = FALSE)
      IRanges::as.data.frame(unlist(vmp))
    }, seqs = seqs, percentID = percentID))
  
  alignments$seqLength <- width(seqs)[
    match(alignments$names, names(seqs))]
  alignments <- alignments[
    alignments$width == nchar(trimSequence) |
      alignments$end == alignments$seqLength,]
  
  alnRanges <- GenomicRanges::reduce(GenomicRanges::GRanges(
    seqnames = alignments$names,
    ranges = IRanges::IRanges(
      start = alignments$start,
      end = alignments$seqLength),
    strand = rep("*", nrow(alignments))),
    ignore.strand = TRUE)
  
  trimmedSeqs <- Biostrings::DNAStringSet(
    seqs[GenomicRanges::seqnames(alnRanges)],
    start = 1L,
    end = GenomicRanges::start(alnRanges)-1)
  
  allSeqs <- c(trimmedSeqs, seqs[!names(seqs) %in% names(trimmedSeqs)])
  allSeqs <- allSeqs[names(seqs)]
  if(noNames) names(allSeqs) <- NULL
  return(allSeqs)
}
