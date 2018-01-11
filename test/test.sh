Rscript seqTrim.R test/testSeq-1.R2.fastq.gz -o test/testSeq-1.R2.trim.fastq.gz \
  -l ACATATGACAACTCAATTAAACGCGAGC --leadMisMatch 3 \
  -r AGATCGGAAGAGCGTCGTGT --overMisMatch 4 --overMaxLength 20 --compress -c 2
