splitNreads
===========

Reads SAM format from stdin or SAM/BAM from first argument.

Parses any reads with splice junctions (N in CIGAR) and splits them into
multiple reads. 0x40 and 0x80 is set in FLAG for the split fragments to 
indicate "chimeric" alignments. Output is uncompressed SAM format that 
needs parsing for read group (RG:A may be produced instead of RG:Z), 
sorting and indexing.

Usage:

```samtools view -h my.bam | SplitNReads.py | sed 's/RG:A:/RG:Z:/g' | samtools calmd -Su - /path/to/genome.fa | samtools sort - out```

```splitNReads.py my.bam | sed 's/RG:A:/RG:Z:/g' | samtools calmd -Su - /path/to/genome.fa | samtools sort - out```

Requirements:
pysam
