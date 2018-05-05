#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
suppressMessages(library("argparse"))
suppressMessages(library("pander"))
panderOptions("table.style", "simple")

code_dir <- dirname(
  sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))

#' Set up and gather command line arguments
parser <- ArgumentParser(
  description = "R-based nucleotide sequence trimmer. Trim both leading and overreading ends of sequences.")
parser$add_argument(
  "seqFile", nargs = 1, type = "character",
  help = "Sequence file to trim, either fasta or fastq format.")
parser$add_argument(
  "-o", "--output", nargs = 1, type = "character", help = "Output file name.")
parser$add_argument(
  "-l", "--leadTrimSeq", nargs = 1, type = "character", default = ".",
  help = "Sequence to trim from 5' end of reads, or the leading sequence. See README for sequence flexibility.")
parser$add_argument(
  "-r", "--overTrimSeq", nargs = 1, type = "character", default = ".",
  help = "Sequence to trim from 3' end of reads, or the overreading sequence. See README for sequence flexibility.")
parser$add_argument(
  "--phasing", nargs = 1, type = "integer", default = 0, 
  help = "Number of nucleotides to remove from 5' end of sequence before trimming. Default = 0.")
parser$add_argument(
  "--maxMisMatch", nargs = 1, type = "integer",
  help = "Maximum allowable mismatches in leading or overreading trim sequences.")
parser$add_argument(
  "--leadMisMatch", nargs = "+", type = "integer", default = 0,
  help = "Maximum allowable mismatches in leading trim sequence. Default = 0. If multiple matching regions for leading seq, separate mismatches allowed using a space.")
parser$add_argument(
  "--overMisMatch", nargs = 1, type = "integer", default = 0,
  help = "Maximum allowable mismatches in overreading trim sequence. Default = 0.")
parser$add_argument(
  "--overMaxLength", nargs = 1, type = "integer", default = 20,
  help = "Maximum length to consider of the overTrimSeq to use for alignments. See README for in depth explanation of this feature. Default 20 nts.")
parser$add_argument(
  "--overMinLength", nargs = 1, type = "integer", default = 3,
  help = "Minimum length to consider of the overTrimSeq to use for alignments. See README for in depth explanation of this feature. Default 3 nts.")
parser$add_argument(
  "--minSeqLength", nargs = 1, type = "integer", default = 30,
  help = "Minimum length of trimmed sequence. Any trimmed sequence with a length below this value will be filtered out. Default = 30")
parser$add_argument(
  "--collectRandomIDs", nargs = "+", type = "character", default = FALSE,
  help = "Option to collect random nucleotide sequences from trimmed portions. If used, provide an output file name.")
parser$add_argument(
  "--noFiltering", action = "store_true",
  help = "Will not filter reads based on leadTrimSeq, the default behavior.")
parser$add_argument(
  "--noQualTrimming", action = "store_true",
  help = "Will not quality trim reads, the default behavior.")
parser$add_argument(
  "--badQualBases", nargs = 1, type = "integer", default = 5,
  help = "Number of bases below threshold in sliding window before read is trimmed. Default = 5.")
parser$add_argument(
  "--qualSlidingWindow", nargs = 1, type = "integer", default = 10,
  help = "Slinding window size for which to assess quality scores below threshold. Default = 10.")
parser$add_argument(
  "--qualThreshold", nargs = 1, type = "character", default = '?',
  help = "Quality threshold for trimming, minimum allowable score. Default = '?', Q30.")
parser$add_argument(
  "--compress", action = "store_true", help = "Output fastq files are gzipped.")
parser$add_argument(
  "-c", "--cores", nargs = 1, default = 1, type = "integer", 
  help = "Max cores to be used. If 0 or 1 (default), program will not utilize parallel processing.")

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

if(is.null(args$seqFile)){
  stop("Please choose a sequence file (fasta or fastq).")
}

if(!is.null(args$maxMisMatch)){
  args$leadMisMatch <- args$maxMisMatch
  args$overMisMatch <- args$maxMisMatch
}

if(args$overMaxLength == 0){
  args$overMaxLength <- nchar(args$overTrimSeq)
}

if(all(args$collectRandomIDs != FALSE)){
  if(!grepl("N", args$leadTrimSeq)){
    message("No random nucleotides (Ns) found in leadTrimSeq. Turning off collection of randomIDs.")
    args$collectRandomIDs <- FALSE
  }
}

if(args$leadTrimSeq == "."){
  args$leadTrimSeq <- ""
}

if(args$overTrimSeq == "."){
  args$overTrimSeq <- ""
}

if(args$cores == 0){
  args$cores <- 1
}

input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply(1:length(args), function(i){
    paste(args[[i]], collapse = ", ")}))
input_table <- input_table[
  match(c("seqFile :", "output :", "leadTrimSeq :", "overTrimSeq :", 
          "phasing :", "maxMisMatch :", "leadMisMatch :", "overMisMatch :", 
          "overMaxLength :", "overMinLength :", "minSeqLength :", 
          "collectRandomIDs :", "ignoreAmbiguousNts :", "noFiltering :", 
          "noQualTrimming :", "badQualBases :", "qualSlidingWindow :", 
          "qualThreshold :", "compress :", "cores :"),
        input_table$Variables),]
pandoc.title("seqTrimR Inputs")
pandoc.table(data.frame(input_table, row.names = NULL), 
             justify = c("left", "left"), 
             split.tables = Inf)

# Load additional R-packages
if(args$cores > 1){
  addPacks <- c("stringr", "ShortRead", "parallel")
}else{
  addPacks <- c("stringr", "ShortRead")
}

addPacksLoaded <- suppressMessages(
  sapply(addPacks, require, character.only = TRUE))
if(!all(addPacksLoaded)){
  pandoc.table(data.frame(
    "R-Packages" = names(addPacksLoaded), 
    "Loaded" = addPacksLoaded, 
    row.names = NULL))
  stop("Check dependancies.")
}

if(args$cores > 1){
  if(args$cores > parallel::detectCores()){
    message("Requested cores is greater than availible for system. Changing cores to max allowed.")
    args$cores <- detectCores()
  }
}
  
# Load supporting scripts
source(file.path(
  code_dir, "supporting_scripts", "trim_leading.R"))
source(file.path(
  code_dir, "supporting_scripts", "trim_overreading.R"))
source(file.path(
  code_dir, "supporting_scripts", "write_seq_files.R"))
source(file.path(
  code_dir, "supporting_scripts", "utility_funcs.R"))
if(!all(
  c("trim_leading", "trim_overreading", "write_seq_files", "log_seq_data", 
    "serial_append_S4") %in% ls())){
    stop("Cannot load supporting scripts. You may need to clone from github again.")
}

# Determine sequence file type
seqType <- str_extract(args$seqFile, "fa[\\w]*")
if(!seqType %in% c("fa", "fasta", "fastq")){
  stop("Unrecognized sequence file type, please convert to '*.fasta' or '*.fastq'. Gzip compression is acceptable as well.")
}
seqType <- ifelse(seqType %in% c("fa", "fasta"), "fasta", "fastq")

# Determine sequence output file type
outType <- str_extract(args$output, "fa[\\w]*")
if(!outType %in% c("fa", "fasta", "fastq")){
  stop("Unrecognized output file type, please choose '*.fasta' or '*.fastq'.")
}
outType <- ifelse(outType %in% c("fa", "fasta"), "fasta", "fastq")

# Determine random output file type
if(all(args$collectRandomIDs != FALSE)){
  randomType <- str_extract(args$collectRandomIDs, "fa[\\w]*")
  if(all(!randomType %in% c("fa", "fasta", "fastq"))){
    stop("Unrecognized randomID output file type, please choose '*.fasta' or '*.fastq'.")
  }
  randomType <- ifelse(randomType %in% c("fa", "fasta"), "fasta", "fastq")
}

# Read sequence file
if(seqType == "fasta"){
  seqs <- ShortRead::readFasta(args$seqFile)
}else{
  seqs <- ShortRead::readFastq(args$seqFile)
}

# Log info
input_tbl <- log_seq_data(seqs)
pandoc.table(input_tbl, caption = "Input sequence information.")

# If no reads remaining, terminate and write output
if(length(seqs) == 0){
  message(
    "No reads remaining to trim. Terminating script after writing output.")
  write_seq_files(
    seqs = seqs, seqType = outType, file = args$output,
    compress = args$compress)
  q()
}

# Quality trimming, trim from left to remove consecutive bad quality bases.
## Below block sets the OpenMP threads to the cores specified in args.
if(!args$noQualTrimming & seqType == "fastq"){
  nthreads <- .Call(ShortRead:::.set_omp_threads, as.integer(args$cores))
  on.exit(.Call(ShortRead:::.set_omp_threads, nthreads))

  seqs <- ShortRead::trimTailw(
    object = seqs, 
    k = args$badQualBases, 
    a = args$qualThreshold, 
    halfwidth = round(args$qualSlidingWindow/2))
  
  # Log info
  qual_trimmed_tbl <- log_seq_data(seqs)
  pandoc.table(
    qual_trimmed_tbl, 
    caption = "Sequence information remaining after quality trimming.")
}

# If no reads remaining, terminate and write output
if(length(seqs) == 0){
  message(
    "No reads remaining to trim. Terminating script after writing output.")
  write_seq_files(
    seqs = seqs, seqType = outType, file = args$output,
    compress = args$compress)
  q()
}

# Trim sequences, either on a single core or multiple cores
if(args$cores <= 1){
  # Trim 5' end or leading end. Conditionals present for added features.
  if(nchar(args$leadTrimSeq) > 0){
    trimmedSeqs <- trim_leading(
      seqs,
      trim.sequence = args$leadTrimSeq,
      phasing = args$phasing,
      max.mismatch = args$leadMisMatch,
      collect.random = all(args$collectRandomIDs != FALSE),
      filter = !args$noFiltering)
  }else{
    trimmedSeqs <- seqs
  }
  
  # Collect random sequences if desired.
  if(all(args$collectRandomIDs != FALSE)){
    randomSeqs <- trimmedSeqs$randomSequences
    trimmedSeqs <- trimmedSeqs$trimmedSequences
  }
  
  # Log info
  lead_trimmed_tbl <- log_seq_data(trimmedSeqs)
  pandoc.table(
    lead_trimmed_tbl, 
    caption = "Sequence information remaining after lead trimming.")
  
  if(nchar(args$overTrimSeq) > 0){
    # Determine percent identity from allowable mismatch.
    percentID <- (nchar(args$overTrimSeq) - args$overMisMatch) / 
      nchar(args$overTrimSeq)
  
    # Trim 3' end or overreading protion of sequences.
    trimmedSeqs <- trim_overreading(
      trimmedSeqs, 
      trimSequence = args$overTrimSeq, 
      percentID = percentID, 
      maxSeqLength = args$overMaxLength,
      minSeqLength = args$overMinLength)
    
    # Log info
    over_trimmed_tbl <- log_seq_data(trimmedSeqs)
    pandoc.table(
      over_trimmed_tbl, 
      caption = "Sequence information remaining after overreading trimming.")
  }
}else{
  # Split sequences up evenly across cores for trimming
  split.seqs <- split(
    seqs, ceiling(seq_along(seqs)/(length(seqs)/args$cores)))
  # Set up buster the cluster
  buster <- parallel::makeCluster(args$cores)
  
  # Trim 5' end or leading section of sequence while capturing random sequences,
  # if desired. Added features required workflow changes.
  if(nchar(args$leadTrimSeq) > 0){
    trimmedSeqs <- parLapply(
      buster,
      split.seqs,
      trim_leading,
      trimSequence = args$leadTrimSeq,
      phasing = args$phasing,
      maxMisMatch = args$leadMisMatch,
      collectRandomID = all(args$collectRandomIDs != FALSE),
      ignoreAmbiguousNts = args$ignoreAmbiguousNts,
      noFiltering = args$noFiltering
    )
  
    if(all(args$collectRandomIDs != FALSE)){
      randomSeqs <- lapply(trimmedSeqs, "[[", "randomSequences")
      randomSeqs <- lapply(1:length(randomSeqs[[1]]), function(i){
        serial_append_S4(
          lapply(1:length(randomSeqs), function(j){
            randomSeqs[[j]][[i]]
        }))
      })

      trimmedSeqs <- lapply(trimmedSeqs, "[[", "trimmedSequences")
    }
  }else{
    trimmedSeqs <- split.seqs
  }

  trimmedSeqs <- serial_append_S4(trimmedSeqs)
  
  # Log info
  lead_trimmed_tbl <- log_seq_data(trimmedSeqs)
  pandoc.table(
    lead_trimmed_tbl,
    caption = "Sequence information remaining after lead trimming.")
  
  # The method for overread trimming sequentially aligns shorter fragments of 
  # the overTrimSeq, and solely requiring mismatches could lead to some issues.
  # Therefore the same percent identity is requried across all alignments, 
  # however long.
  if(nchar(args$overTrimSeq) > 0){  
    trimmedSeqs <- split(
      trimmedSeqs, 
      ceiling(seq_along(trimmedSeqs)/(length(trimmedSeqs)/args$cores)))
    
    percentID <- (nchar(args$overTrimSeq) - args$overMisMatch) / 
      nchar(args$overTrimSeq)
  
    # Trim 3' end or overreading protion of the sequence.
    trimmedSeqs <- parLapply(
      buster,
      trimmedSeqs,
      trim_overreading,
      trimSequence = args$overTrimSeq, 
      percentID = percentID, 
      maxSeqLength = args$overMaxLength,
      minSeqLength = args$overMinLength)

    trimmedSeqs <- serial_append_S4(trimmedSeqs)
    
    # Log info
    over_trimmed_tbl <- log_seq_data(trimmedSeqs)
    pandoc.table(
      over_trimmed_tbl, 
      caption = "Sequence information remaining after overreading trimming.")
  }
  # Stop buster before he gets out of control.
  stopCluster(buster)
}

# If no reads remaining, terminate and write output
if(length(seqs) == 0){
  message(
    "No reads remaining to trim. Terminating script after writing output.")
  write_seq_files(
    seqs = trimmedSeqs, seqType = outType, file = args$output,
    compress = args$compress)
  q()
}

# Filter sequences by minimum length.
trimmedSeqs <- trimmedSeqs[width(trimmedSeqs) >= args$minSeqLength]

# Log info
len_trimmed_tbl <- log_seq_data(trimmedSeqs)
pandoc.table(
  len_trimmed_tbl, 
  caption = "Sequence information remaining after minimum length trimming.")

# Collect RandomIDs if requested
if(all(args$collectRandomIDs != FALSE)){
  randomSeqs <- lapply(1:length(randomSeqs), function(i, ids){
    randomSeqs[[i]][which(as.character(ShortRead::id(randomSeqs[[i]])) %in% ids)]
  }, ids = as.character(ShortRead::id(trimmedSeqs)))
}

# Sequences have been trimmed and random sequnces collected (if desired). 
# Next step is to write to output file(s).
# For fasta format, this is as simple as writing out the sequences currently in
# the environment. For fastq format, the quality scores for the trimmed bases
# must be loaded and trimmed as well.

# Write sequence file.
write_seq_files(
  seqs = trimmedSeqs, 
  seqType = outType, 
  file = args$output,
  compress = args$compress)

# Write randomID file.
if(all(args$collectRandomIDs != FALSE)){
  if(length(randomSeqs) == 1){
    write_seq_files(
      seqs = randomSeqs[[1]],
      seqType = randomType,
      file = args$collectRandomIDs,
      compress = args$compress)
  }else{
    if(length(args$collectRandomIDs) != length(randomSeqs)){
      newFileName <- unlist(strsplit(args$collectRandomIDs[[1]], ".fa"))
      newNames <- paste0(
        newFileName[[1]], ".", 1:length(randomSeqs), ".", randomType)
      args$collectRandomIDs <- newNames
    }
    
    null <- mapply(
      write_seq_files,
      seqs = randomSeqs,
      file = args$collectRandomIDs,
      MoreArgs = list(
        seqType = randomType, 
        compress = args$compress)
    )
}}

cat("Script completed.\n")
q()
