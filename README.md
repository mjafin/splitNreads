splitNreads
===========

Reads SAM format from stdin or SAM/BAM from first argument.

Parses any reads with splice junctions (N in CIGAR) and splits them into
multiple reads. 0x40 and 0x80 is set in FLAG for the split fragments to 
indicate "chimeric" alignments. Output is uncompressed BAM format that 
needs sorting and indexing.

Usage:

```samtools view -h my.bam | SplitNReads.py | samtools sort - out```

```splitNReads.py my.bam | samtools sort - out```

Requirements:
pysam
