#!/usr/bin/env bash
set -ev

# Test for leading and overreading trimming on R2 test sequences
Rscript seqTrim.R tests/testSeq-1.R2.fastq.gz -o tests/testSeq-1.R2.trim.fastq.gz \
  -l ACATATGACAACTCAATTAAACGCGAGC --leadMisMatch 3 \
  -r AGATCGGAAGAGCGTCGTGT --overMisMatch 4 --overMaxLength 20 --compress

# Test for only overreading trimming on R1 test sequences
Rscript seqTrim.R tests/testSeq-1.R1.fastq.gz -o tests/testSeq-1.R1.trim.fastq.gz \
    -r GCTCGCGTTTAATTGAGTTGTCATATGT --overMisMatch 4 --overMaxLength 20 --compress

# Check outputs for correct findings
test_R2_len=$(zcat tests/testSeq-1.R2.trim.fastq.gz | sed -n '2~4p' | wc -l)
test_R1_len=$(zcat tests/testSeq-1.R1.trim.fastq.gz | sed -n '2~4p' | wc -l)
if [ ! $test_R2_len = 50 ] | [ ! $test_R1_len = 50 ]; then
    exit 1
fi

# R2 test sequences
zcat tests/testSeq-1.R2.trim.fastq.gz | sed -n '2~4p' | head -n 5

# R1 test sequences
zcat tests/testSeq-1.R1.trim.fastq.gz | sed -n '2~4p' | head -n 5

# Cleanup
rm -f tests/*.trim.fastq*

echo "Passed all tests."
exit