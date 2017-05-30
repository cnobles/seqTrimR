# seqTrimR
Trims 5' and 3' nucleotide sequences from paired-end reads based on designated sequences and/or quality scores.
--------------------------------------------------------------------------------
Usage:
```
Rscript seqTrim.R seqFile.fastq -o output.fastq -l leadingTrimSequence --maxMisMatch 3

Rscript seqTrim.R seqFile.fastq -o output.fastq -l leadingTrimSequence -r overReadingTrimSequence \
    --phasing 8 --leadMisMatch 2 --overMisMatch 1 --overMaxLength 20 --collectRandomIDs --compress --cores 10
```

##Arguments
**[seqFile]** Sequence file to trim, either fasta or fastq (gzip compression tolerated).
**[-h, --help]** Show help message and exit.
**[-o,--output]** Output file name.
**[-l, --leadTrimSeq]** Sequence to trim from 5' end of reads, or the leading sequence. See above for sequence flexibility.
**[-r, --overTrimSeq]** Sequence to trim from 3' end of reads, or the overreading sequence. See above for sequence flexibility.
**[--phasing]** Number of nucleotides to remove from 5' end of sequence before trimming. Default = 0.
**[--maxMisMatch]** Maximum allowable mismatches in leading or overreading trim sequences.
**[--leadMisMatch]** Maximum allowable mismatches in leading trim sequence. Default = 0.
**[--overMisMatch]** Maximum allowable mismatches in overreading trim sequence. Default = 0.
**[--overMaxLength]** Maximum length to consider of the overTrimSeq to use for alignments. See above for in-depth explanation of this feature.
**[--minSeqLength]** Minimum length of trimmed sequence. Any trimmed sequence with a length below this value will be filtered out. Default = 30.
**[--collectRandomIDs]** Option to collect random nucleotide sequences from trimmed portions. If used, provide an output file name.
**[--ignoreAmbiguousNts]**  Conversely, ambiguous nucleotides can be ignored from collection but still enforced in matching for trimming.
**[--noFiltering]** Will not filter reads based on leadTrimSeq, the default behavior.
**[--compress]** Output fastq/fasta files are gzipped.
**[-c, --cores]** Max cores to be used. If 0 (default), program will not utilize parallel processing.

##Dependencies
seqTrimR is coded in R, and was developed on v3.2.2, though it should run with earlier versions given the appropriate dependencies. The script uses 9 additional packages:
  * argparse
  * pander
  * stringr
  * ShortRead
  * Biostrings
  * BiocGenerics
  * IRanges
  * GenomicRanges
  * parallel (if multicore processing is desired)
