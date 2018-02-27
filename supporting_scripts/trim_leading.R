#' Trim beginning or leading ends of nucleotide sequences
#'
#' @param seqs ShortReadQ object of reads or unique sequences
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
  # Checks and requirements
  suppressMessages(require(BiocGenerics))
  suppressMessages(require(Biostrings))
  stopifnot(class(seqs) %in% c("ShortReadQ", "ShortRead"))
  stopifnot(!is.null(ShortRead::id(seqs)))
  
  # Change scientific notation switch to inhibit indexing errors
  ori.scipen <- getOption("scipen")
  options(scipen = 99)
  
  if(ignoreAmbiguousNts & all(collectRandomID != FALSE)){
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
  seqs <- narrow(seqs, start = 1L + phasing)
  
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
  
  # Remove seqs that do not have enough sequence for analysis
  # Cutoff = length(trimSequence)
  seqs <- seqs[width(seqs) >= nchar(trimSequence)+2]
  
  # Serially align the segment(s) from trimSequence to seqs
  aln <- lapply(1:length(tSegRanges), function(i, tSegRanges, seqs){
      tSeq <- names(tSegRanges[i])
      misMatch <- tSegRanges@metadata$misMatch[i]
      alnSeqs <- narrow(
        seqs,
        start = ifelse(start(tSegRanges[i]) == 1L, 1, start(tSegRanges[i]) - 1),
        end = end(tSegRanges[i]) + 1)
      aln <- vmatchPattern(
        tSeq, ShortRead::sread(alnSeqs), max.mismatch = misMatch, fixed = FALSE)
    
      if(any(lengths(aln) > 1)){
        stop("\nAlignment too permissive. Ambiguous mapping of sequences.
           Please adjust maxMisMatch criteria.")}
      idx <- lengths(aln) == 1 
      return(list("match" = aln, "idx" = idx))
    },
    tSegRanges = tSegRanges,
    seqs = seqs)

  # Identify and trim only sequences that have all required alignments
  matchedIndex <- table(unlist(lapply(lapply(aln, "[[", "idx"), which)))
  matchedIndex <- as.numeric(names(matchedIndex)[
    matchedIndex == length(tSegRanges)])
  matchedSeqs <- seqs[matchedIndex]
  tShift <- ifelse(length(tSegRanges) > 1, 2, 1)
  matchedStarts <- start(tail(tSegRanges, n = 1)) - tShift + 
    unlist(endIndex(aln[[length(aln)]]$match)) + 1
  trimmedSeqs <- narrow(matchedSeqs, start = matchedStarts)

  if(noFiltering){
    unmatchedIndex <- which(!1:length(seqs) %in% matchedIndex)
    untrimmedSeqs <- seqs[unmatchedIndex]
    trimmedSeqs <- append(trimmedSeqs, untrimmedSeqs)
    trimmedSeqs <- trimmedSeqs[order(c(matchedIndex, unmatchedIndex))]
  }
  
  if(!all(collectRandomID != FALSE) | !grepl("N", trimSequence)){
    # Return scipen option to original value
    options(scipen = ori.scipen)
    
    return(trimmedSeqs)
  }else{
    randomSets <- lapply(1:(length(aln)-1), function(k, aln, matchedIndex){
      gapRanges <- do.call(c, lapply(k:(k+1), function(i){
        matched <- aln[[i]]$match
        idx <- which(lengths(matched) == 1)
        tRanges <- unlist(matched)
        names(tRanges) <- idx
        tShift <- ifelse(i > 1, start(tSegRanges[i]) - 1, start(tSegRanges[i]))
        shift(tRanges, shift = tShift - 1)
      }))
      gapRanges <- gapRanges[names(gapRanges) %in% matchedIndex]
      gapRanges <- split(gapRanges, names(gapRanges))
      gapRanges <- unlist(gaps(gapRanges))
      gapRanges[as.character(matchedIndex)]
    }, aln = aln, matchedIndex = matchedIndex)
    
    randomSeqs <- lapply(randomSets, function(ir, matchedIndex, seqs){
      narrow(seqs[matchedIndex], start = start(ir), end = end(ir))
    }, matchedIndex = matchedIndex, seqs = seqs)
    
    # Return scipen option to original value
    options(scipen = ori.scipen)
    
    return(list(
      "trimmedSequences" = trimmedSeqs,
      "randomSequences" = randomSeqs))
  }
}
