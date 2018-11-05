#!/usr/bin/env bash
set -ev

# Test for leading and overreading trimming on R2 test sequences
Rscript seqTrim.R tests/testSeq-1.R2.fastq.gz -o tests/testSeq-1.R2.trim.fastq.gz \
  -l ACATATGACAACTCAATTAAACGCGAGC --leadMismatch 3 \
  -r AGATCGGAAGAGCGTCGTGT --overMismatch 4 --overMaxLength 20 \
  --stat tests/test.R2.stat.csv --compress

# Test for only overreading trimming on R1 test sequences
Rscript seqTrim.R tests/testSeq-1.R1.fastq.gz -o tests/testSeq-1.R1.trim.fastq.gz \
  -r GCTCGCGTTTAATTGAGTTGTCATATGT --overMismatch 4 --overMaxLength 20 \
  --stat tests/test.R1.stat.csv --compress

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

# Concatenate the stat files
cat tests/test.R1.stat.csv tests/test.R2.stat.csv

# Cleanup
rm -f tests/*.trim.fastq* tests/test.R?.stat.csv

echo "Passed all tests."
exit
