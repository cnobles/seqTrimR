# seqTrimR
Trims 5' and 3' nucleotide sequences from paired-end reads based on designated sequences and/or quality scores.
--------------------------------------------------------------------------------
Trimming diagram:
```
         leading trimming                    overreading trimming
         ----------CA                                    TCA-----
         ||||||||||||                                    ||||||||    
[seq] 5' ----------CAAGTC----------------------------TCCATCA-----3'

[result] AGTC----------------------------TCCA
```
Random ID collection diagram:
```
         leading trimming                    overreading trimming
         ----NNNN--CA                                    TCA-----
         ||||    ||||                                    ||||||||    
[seq] 5' ----CGCT--CAAGTC----------------------------TCCATCA-----3'

[result] AGTC----------------------------TCCA
[random] CGCT
```
Usage:
```
Rscript seqTrim.R seqFile.fastq -o output.fastq -l leadingTrimSequence 

Rscript seqTrim.R seqFile.fastq -o output.fastq -l leadingTrimSequence \
  -r overReadingTrimSequence --phasing 8 --leadMisMatch 2 --overMisMatch 1 \
  --overMaxLength 20 --collectRandomIDs --compress --cores 10
```
## Trimming sequence flexibility
The arguments **[leadTrimSeq]** and **[overTrimSeq]** take character string inputs, such as below.
```
# [leadTrimSeq] examples
ATGCCGTTAGCTATGC	    #Fixed nucleotide sequence
ATGCCGTTNNNNNNAGCTATGC	    #Random 6 nucleotide barcode embedded in sequence
ATGCCGNNNNTTAGNNNNCTATGC    #Dual 4 nucleotide barcodes embedded in sequence
# Note: Random nucleotides from within leading trim sequences containing random 
#   embeded nucleotides can be collected using the [collectRandomIDs] argument.

# [overTrimSeq] examples
GCTAACGTAC			  #Short 10 nucleotide fixed sequence
GCTAACGTACGTTTCAAGCTACGGACATGC    #Longer nucleotide fixed sequence
GCTAACGTACGTTTCNNNNNNCGGACATGC    #Longer nucleotide sequence with 
				  #  embedded random sequence
# Note: Longer overreading sequences can take longer to process, alternatively
#   using the argument [overMaxLength] will limit the number of nucleotides
#   from the beginning to improve performance. Using this argument does not 
#   change the allowed mismatch parameter, which will always be based on the
#   full length sequence.
```
Using the argument **[collectRandomIDs]** will collect the random sequences from within the leading trimming sequence whenever a match is made. Random sequences are returned in the order in which they appear from 5' to 3' (left to right). Multiple output file names can be given to name each random sequence file captured after the **[collectRandomIDs]** flag. 

Additionally, if random or ambiguous nucleotide sequences are to be used in matching, the argument **[ignoreAmbiguousNts]** can be used to ignore collection of random sequences. The pattern matching will follow the ambiguous nucleotide code.

## Arguments
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

## Dependencies
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
