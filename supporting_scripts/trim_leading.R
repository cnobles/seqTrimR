#' Trim beginning or leading ends of nucleotide sequences
#'
#' @param seqs DNAStringSet of reads or unique sequences
#' @param trimSequence character string of lenth 1, such as "GAAAATC". This
#' string will be used to match to the beginning of sequences, upon which
#' non-matching sequences will be discarded and the matching portion will be
#' trimmed from the leading side of the sequence. Ambiguous nucleotides within
#' the sequence will be used to determine random sequences and can be used for
#' alignment or collecting random nucleotide sequences embedded within the
#' trimSequence structure.
#' @param phasing integer/numeric value, denoting the number of nucleotides used
#' for phasing while sequencing. This number of nucleotides will be removed from
#' the beginning of the sequence before any alignment.
#' @param maxMisMatch integer/numeric value or vector. Values indicate the
#' number of mismatches allowed within the DNA segment. Vectors must be equal in
#' length to the number of non-ambiguous segments within the trimSequence. If no
#' ambiguous segments are present in trimSequence or ignoreAmbiguousNts is TRUE,
#' vectors for maxMisMatch will be summed to determine the maximum allowable
#' mismatch over the entire string.
#' @param collectRandomID logical should random / ambiguous protions of the
#' trimSequence be collected? They will be returned as a listed DNAStringSet
#' under 'randomSequences' in order from left to right of appearance within
#' trimSequence.
#' @param ignoreAmbiguousNts logical To ignore ambiguous nucleotides within the
#' trimSequence string. If theses nucleotides are ignored, then randomIDs cannot
#' currently be collected. Rather, adjust maxMisMatch to obtain a suitable
#' alignment for non-ambiguous segments.
#' @return DNAStringSet of sequences with trimSequence removed or a listed
#' object with the first postion being the DNAStringSet of trimmed sequences and
#' the second being the random sequences collected during trimming.
#' @author Christopher Nobles, Ph.D.

trim_leading <- function(seqs, trimSequence, phasing = 0L, maxMisMatch = 1L,
                         collectRandomID = FALSE, ignoreAmbiguousNts = FALSE,
                         noFiltering = FALSE){
  require(BiocGenerics)
  require(Biostrings)
  stopifnot(class(seqs) %in% "DNAStringSet")
  stopifnot(!is.null(names(seqs)))
  if(ignoreAmbiguousNts & collectRandomID){
    message("\nCurrently this function cannot collect random IDs
            if it ignores ambiguous nucleotides.
            Switching collectRandomID to FALSE.")
    collectRandomID <- FALSE
  }
  if(length(maxMisMatch) > 1 & ignoreAmbiguousNts){
    message("\nSum of maxMisMatch is being used for maxMisMatch since
            ignoreAmbiguousNts is being chosen.")
    maxMisMatch <- sum(maxMisMatch)
  }
  
  # Phasing will ignore the first number of nucleotides of the sequence
  seqs <- Biostrings::DNAStringSet(
    seqs,
    start = 1L + phasing)
  
  # Determine the structure of random sequences within the trimSequence
  if(!ignoreAmbiguousNts){
    trimSegments <- unlist(strsplit(trimSequence, "[N]+"))
    tSegRanges <- sapply(trimSegments, vmatchPattern, subject = trimSequence)
    tSegRanges <- IRanges(
      start = sapply(tSegRanges, function(x) as.integer(IRanges::start(x))),
      end = sapply(tSegRanges, function(x) as.integer(IRanges::end(x))),
      names = names(tSegRanges))
  }else{
    trimSegments <- trimSequence
    tSegRanges <- unlist(vmatchPattern(trimSegments, trimSequence))
    names(tSegRanges) <- trimSegments
  }
  # Set allowable mismatches for each range
  if(length(maxMisMatch) == 1 | ignoreAmbiguousNts){
    tSegRanges@metadata$misMatch <- rep(maxMisMatch, length(tSegRanges))
  }else if(length(maxMisMatch) == length(tSegRanges) & is.numeric(maxMisMatch)){
    tSegRanges@metadata$misMatch <- maxMisMatch
  }else{
    stop("\nThe variable maxMisMatch needs to be either a
         single integer or integer vector of length equal
         to fixed fragments within the trimSequence.")
  }
  
  # Serially align the segment(s) from trimSequence to seqs
  aln <- do.call(c, lapply(1:length(tSegRanges), function(i, tSegRanges, seqs){
    tSeq <- names(tSegRanges[i])
    misMatch <- tSegRanges@metadata$misMatch[i]
    alnSeqs <- DNAStringSet(
      seqs,
      start = ifelse(start(tSegRanges[i]) == 1L, 1, start(tSegRanges[i]) - 1),
      end = end(tSegRanges[i]) + 1)
    aln <- unlist(vmatchPattern(
      tSeq, alnSeqs, max.mismatch = misMatch, fixed = FALSE))
    shift(
      aln,
      shift = ifelse(
        start(tSegRanges[i]) == 1L, 0L, start(tSegRanges[i]) - 2))
  },
  tSegRanges = tSegRanges,
  seqs = seqs))
  
  # Identify and select only sequences that have all required alignments
  seqCount <- table(names(aln))
  if(any(seqCount > length(trimSegments))){
    stop("\nAlignment too permissive. Ambiguous mapping of sequences.
         Please adjust maxMisMatch criteria.")
  }
  matchedSeqs <- seqs[names(seqCount[seqCount == length(tSegRanges)])]
  aln <- aln[names(aln) %in% names(seqCount[seqCount == length(tSegRanges)])]
  matchedRanges <- split(aln, names(aln))
  trimRanges <- unlist(reduce(matchedRanges, min.gapwidth = nchar(trimSequence)))
  
  # Trim sequences with the trimSequence alignment position(s)
  trimmedSeqs <- DNAStringSet(
    matchedSeqs[names(trimRanges)],
    start = end(trimRanges) + 1)
  
  if(noFiltering){
    trimmedSeqs <- c(trimmedSeqs, seqs[!names(seqs) %in% names(trimmedSeqs)])
    trimmedSeqs <- trimmedSeqs[names(seqs)]
  }
  
  if(!collectRandomID | !grepl("N", trimSequence)){
    return(trimmedSeqs)
  }else{
    gapRanges <- gaps(matchedRanges)
    randomSeqs <- lapply(
      1:(length(trimSegments) - 1), function(i, gapRanges, matchedSeqs){
        starts <- sapply(1:length(gapRanges), function(j) start(gapRanges[[j]][i]))
        ends <- sapply(1:length(gapRanges), function(j) end(gapRanges[[j]][i]))
        DNAStringSet(
          matchedSeqs[names(gapRanges)],
          start = starts,
          end = ends)
      },
      gapRanges = gapRanges,
      matchedSeqs = matchedSeqs)
    return(list(
      "trimmedSequences" = trimmedSeqs,
      "randomSequences" = randomSeqs))
  }
}
