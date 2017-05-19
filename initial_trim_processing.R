# Initial processing of read
data_dir <- "~/data/projects/guideseq_analysis/results/170411_second_sample_set_neg_strand/demultiplexed"
code_dir <- "~/dev/scripts/seqTrimR"
output_dir <- "~/data/projects/guideseq_analysis/results/170411_second_sample_set_neg_strand/trimmed"

if(!dir.exists(output_dir)){
  attempt <- try(system(paste0("mkdir ", output_dir)))
  if(attempt == 1) stop("Failed to make output directory.")
}

libs <- c("Biostrings", "GenomicRanges", "ShortRead", "dplyr", "pander",
          "stringr", "yaml", "parallel")
loaded <- suppressMessages(sapply(libs, require, character.only = TRUE))
if(all(loaded)){
  message("All dependancies loaded.")
}else{
  message("Check dependencies.")
  pandoc.table(loaded)
}

r1.lead <- c("")
r2.lead <- DNAStringSet(c("TTGAGTTGTCATATGTTAATAACGGTAT"))
r1.over <- DNAStringSet(reverseComplement(r2.lead), end = 15)
r2.over <- DNAStringSet(reverseComplement(DNAStringSet(
  c("ACACTCTTTCCCTACACGACGCTCTTCCGATCT"))), end = 15)
percentID <- 0.95
compress <- TRUE
min_length <- 30
filter <- TRUE
readNamePattern <- "[\\w:-]+"

banmat_path <- file.path(
  code_dir, "supporting_scripts/binary_ambiguous_nucleotide_scoring_matrix.R")
trim_lead <- file.path(code_dir, "supporting_scripts/trim_leading.R")
trim_over <- file.path(code_dir, "supporting_scripts/trim_overreading.R")

samples <- unique(str_extract(list.files(data_dir), "^[\\w]+"))
samples <- samples[!samples %in% c("ambiguous", "unassigned")]

buster <- makeCluster(8)

clusterExport(
  cl = buster, varlist = c("banmat_path", "trim_lead", "trim_over"))

process <- parLapply(
  buster,
  samples, 
  function(sample, r1.lead, r2.lead, r1.over, r2.over, 
           percentID, data_dir, output_dir, filter, min_length, readNamePattern, 
           compress){
    libs <- c("Biostrings", "GenomicRanges", "ShortRead", "dplyr", "pander",
              "stringr", "yaml", "parallel")
    loaded <- suppressMessages(sapply(libs, require, character.only = TRUE))
    stopifnot(all(loaded))
    
    source(banmat_path)
    source(trim_lead)
    source(trim_over)
    
    files <- list.files(data_dir, pattern = paste0(sample, ".r[1-2]?.fastq.gz"))

    # Trim R1 sequences
    r1 <- readFastq(
      file.path(data_dir, grep("r1.fastq.gz", files, value = TRUE)))
    r1.seqs <- sread(r1)
    r1.qual <- quality(r1)
    names(r1.seqs) <- names(r1.qual@quality) <- ShortRead::id(r1)
    r1.seqs <- trim_leading(r1.lead, r1.seqs, percentID, filter = filter)
    r1.seqs <- trim_overreading(r1.over, r1.seqs, percentID)

    # Trim R2 sequences
    r2 <- readFastq(
      file.path(data_dir, grep("r2.fastq.gz", files, value = TRUE)))
    r2.seqs <- sread(r2)
    r2.qual <- quality(r2)
    names(r2.seqs) <- names(r2.qual@quality) <- ShortRead::id(r2)
    r2.seqs <- trim_leading(r2.lead, r2.seqs, percentID, filter = filter)
    r2.seqs <- trim_overreading(r2.over, r2.seqs, percentID)

    # Remove reads below minimum threshold
    r1.seqs <- r1.seqs[width(r1.seqs) >= min_length]
    r2.seqs <- r2.seqs[width(r2.seqs) >= min_length]
    
    # Filter for only sequences with R1 and R2 reads
    r1.reads <- str_extract(names(r1.seqs), readNamePattern)
    r2.reads <- str_extract(names(r2.seqs), readNamePattern)
    common_reads <- intersect(r1.reads, r2.reads)
    r1.seqs <- r1.seqs[match(common_reads, r1.reads)]
    r2.seqs <- r2.seqs[match(common_reads, r2.reads)]

    # Convert trimming data to quality scores
    r1.qual@quality@ranges <- r1.seqs@ranges
    r2.qual@quality@ranges <- r2.seqs@ranges
        
  if(compress){
    writeXStringSet(
      r1.seqs, 
      file = file.path(output_dir, paste0(sample, ".r1.trimmed.fastq.gz")),
      format = "fastq",
      qualities = r1.qual@quality,
      compress = TRUE)
    writeXStringSet(
      r2.seqs, 
      file = file.path(output_dir, paste0(sample, ".r2.trimmed.fastq.gz")),
      format = "fastq",
      qualities = r2.qual@quality,
      compress = TRUE)
  }else{
    writeXStringSet(
      r1.seqs, 
      file = file.path(output_dir, paste0(sample, ".r1.trimmed.fastq")),
      format = "fastq",
      qualities = r1.qual@quality)
    writeXStringSet(
      r2.seqs, 
      file = file.path(output_dir, paste0(sample, ".r2.trimmed.fastq")),
      format = "fastq",
      qualities = r2.qual@quality)
  }
  sample}, r1.lead = r1.lead, r2.lead = r2.lead, r1.over = r1.over, 
  r2.over = r2.over, percentID = percentID, data_dir = data_dir, 
  output_dir = output_dir, filter = filter, min_length = min_length, 
  readNamePattern = readNamePattern, compress = compress)

stopCluster(buster)
