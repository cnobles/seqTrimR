#' Trim beginning or leading ends of nucleotide sequences
#'
#' @param seqs ShortReadQ object of reads or unique sequences
#' @param trim.sequence character string of lenth 1, such as "GAAAATC". This
#' string will be used to match to the beginning of sequences, upon which
#' non-matching sequences will be discarded and the matching portion will be
#' trimmed from the leading side of the sequence. Ambiguous nucleotides within
#' the sequence will be used to determine random sequences and can be used for
#' alignment or collecting random nucleotide sequences embedded within the
#' trim.sequence structure.
#' @param phasing integer/numeric value, denoting the number of nucleotides used
#' for phasing while sequencing. This number of nucleotides will be removed from
#' the beginning of the sequence before any alignment.
#' @param max.mismatch integer/numeric value or vector. Values indicate the
#' number of mismatches allowed within the DNA segment. Integer / numeric 
#' vectors can be used to indicate the number of allowable mismatches within 
#' segments of non-ambiguous nucleotices. The length of input vecors must be 
#' equal to the number of non-ambiguous segments within the trim.sequence. If 
#' no ambiguous segments are present in trim.sequence vectors for max.mismatch 
#' will be summed to determine the maximum allowable mismatch over the entire 
#' string.
#' @param collect.random logical should random / ambiguous protions of the
#' trim.sequence be collected? They will be returned as a listed DNAStringSet
#' under 'randomSequences' in order from left to right of appearance within
#' trim.sequence.
#' @param filter logical (default: TRUE). If TRUE, sequences not matching the 
#' trim.sequence will be dropped or filtered from the output. If FALSE, all
#' input sequences are returned, but matching sequences are trimmed.
#' 
#' @return DNAStringSet of sequences with trim.sequence removed or a listed
#' object with the first postion being the DNAStringSet of trimmed sequences and
#' the second being the random sequences collected during trimming.
#' 
#' @author Christopher Nobles, Ph.D.
#'  

trim_leading <- function(seqs, trim.sequence, phasing = 0L, max.mismatch = 1L,
                         collect.random = FALSE, filter = TRUE){
  # Checks and requirements
  stopifnot(class(seqs) %in% c("ShortReadQ", "ShortRead"))
  stopifnot(!is.null(ShortRead::id(seqs)))
  
  # Change scientific notation switch to inhibit indexing errors
  ori.scipen <- getOption("scipen")
  options(scipen = 99)
  
  if(length(max.mismatch) > 1){
    message("\nSum of max.mismatch is being used for max.mismatch since
            ignoreAmbiguousNts is being chosen.")
    max.mismatch <- sum(max.mismatch)
  }
  
  # Phasing will ignore the first number of nucleotides of the sequence
  seqs <- narrow(seqs, start = 1L + phasing)
  
  # Determine the structure of ambiguous sequences within the trim.sequence
  ambi_present <- stringr::str_detect(trim.sequence, pattern = "[^A^T^G^C]")
  if(ambi_present & collect.random){
    trim_segments <- unlist(strsplit(trim.sequence, "[^A^T^G^C]+"))
    trim_seg_ir <- sapply(trim_segments, vmatchPattern, subject = trim.sequence)
    trim_seg_ir <- IRanges(
      start = sapply(trim_seg_ir, function(x) as.integer(IRanges::start(x))),
      end = sapply(trim_seg_ir, function(x) as.integer(IRanges::end(x))),
      names = names(trim_seg_ir))
  }else{
    trim_segments <- trim.sequence
    trim_seg_ir <- unlist(vmatchPattern(trim_segments, trim.sequence))
    names(trim_seg_ir) <- trim_segments
  }
  
  # Set allowable mismatches for each range
  if(length(max.mismatch) == 1){
    trim_seg_ir@metadata$misMatch <- rep(max.mismatch, length(trim_seg_ir))
  }else if(
    length(max.mismatch) == length(trim_seg_ir) & is.numeric(max.mismatch)){
      trim_seg_ir@metadata$misMatch <- max.mismatch
  }else{
    stop("\nThe variable max.mismatch needs to be either a
         single integer or integer vector of length equal
         to fixed fragments within the trim.sequence.")
  }
  
  # Remove seqs that do not have enough sequence for analysis
  # Cutoff = length(trim.sequence)
  seqs <- seqs[width(seqs) >= nchar(trim.sequence)]
  lead_seqs <- narrow(
    seqs, 
    start = 1, 
    end = ifelse(
      width(seqs) >= nchar(trim.sequence) + 1, 
      rep(nchar(trim.sequence) + 1, length(seqs)), width(seqs)))
  
  # Align whole sequence to 5' end of sequence
  aln <- vmatchPattern(
    trim.sequence, ShortRead::sread(lead_seqs), 
    max.mismatch = sum(max.mismatch), fixed = FALSE)
  
  matched_idx <- which(lengths(aln) == 1)
  
  # Serially align segment(s) from trim.sequence to seqs
  if(length(trim_seg_ir) > 1){
    aln_segs <- lapply(seq_along(trim_seg_ir), function(i, trim_seg_ir, seqs){
        tSeq <- names(trim_seg_ir[i])
        misMatch <- trim_seg_ir@metadata$misMatch[i]
        alnSeqs <- narrow(
          seqs,
          start = ifelse(
            start(trim_seg_ir[i]) == 1L, 1, start(trim_seg_ir[i]) - 1),
          end = end(trim_seg_ir[i]) + 1)
        aln <- vmatchPattern(
          tSeq, ShortRead::sread(alnSeqs), 
          max.mismatch = misMatch, fixed = FALSE)
      
        if(any(lengths(aln) > 1)){
          stop("\nAlignment too permissive. Ambiguous mapping of sequences.
             Please adjust max.mismatch criteria.")}
        idx <- lengths(aln) == 1 
        return(list("match" = aln, "idx" = idx))
      },
      trim_seg_ir = trim_seg_ir,
      seqs = lead_seqs)
    
    # Identify and trim only sequences that have all required alignments
    seg_idx <- table(unlist(lapply(lapply(aln, "[[", "idx"), which)))
    seg_idx <- as.numeric(names(seg_idx)[seg_idx == length(trim_seg_ir)])
    matched_idx <- matched_idx[matched_idx %in% seg_idx]
  }
  
  # Isolate sequences matching the input criteria
  matched_seqs <- seqs[matched_idx]
  trim_shift <- ifelse(length(trim_seg_ir) > 1, 2, 1)
  matched_starts <- unlist(aln@ends[matched_idx]) + 1
  trimmed_seqs <- narrow(matched_seqs, start = matched_starts)

  if(!filter){
    unmatched_idx <- which(!seq_along(seqs) %in% matched_idx)
    untrimmed_seqs <- seqs[unmatched_idx]
    trimmed_seqs <- append(trimmed_seqs, untrimmed_seqs)
    trimmed_seqs <- trimmed_seqs[order(c(matched_idx, unmatched_idx))]
  }
  
  if(!collect.random | !ambi_present){
    # Return scipen option to original value
    options(scipen = ori.scipen)
    return(trimmed_seqs)
  }else{
    if(length(trim_seg_ir) > 1){
      random_sets <- lapply(
        seq_along(aln_seg), function(k, aln_seg, matched_idx){
          gap_ranges <- do.call(c, lapply(k:(k+1), function(i){
            matched <- aln_seg[[i]]$match
            idx <- which(lengths(matched) == 1)
            t_ranges <- unlist(matched)
            names(t_ranges) <- idx
            trim_shift <- ifelse(
              i > 1, start(trim_seg_ir[i]) - 1, start(trim_seg_ir[i]))
            shift(t_ranges, shift = trim_shift - 1)
          }))
          gap_ranges <- gap_ranges[names(gap_ranges) %in% matched_idx]
          gap_ranges <- split(gap_ranges, names(gap_ranges))
          gap_ranges <- unlist(gaps(gap_ranges))
          gap_ranges[as.character(matched_idx)]
        }, aln_seg = aln_seg, matched_idx = matched_idx)
    }else{
      matched_region <- GRanges(
        seqnames = as.character(matched_idx), ranges = unlist(aln[matched_idx]))
      non_ambi_region <- GRanges(
        seqnames = as.character(matched_idx), 
        ranges = rep(trim_seg_ir, length(matched_idx)))
      random_region <- GenomicRanges::setdiff(matched_region, non_ambi_region)
      random_sets <- ranges(random_region)
      names(random_sets) <- seqnames(random_region)
      random_sets <- list(random_sets[as.character(matched_idx)])
    }
    
    random_seqs <- lapply(random_sets, function(ir, matchedIndex, seqs){
      narrow(seqs[matchedIndex], start = start(ir), end = end(ir))
    }, matchedIndex = matchedIndex, seqs = seqs)
    
    # Return scipen option to original value
    options(scipen = ori.scipen)
    
    return(list(
      "trimmedSequences" = trimmed_seqs,
      "randomSequences" = random_seqs))
  }
}
