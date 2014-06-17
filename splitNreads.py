#!/usr/bin/env python
from __future__ import print_function
import sys
import pysam

"""
M 	BAM_CMATCH 	0
I 	BAM_CINS 	1
D 	BAM_CDEL 	2
N 	BAM_CREF_SKIP 	3
S 	BAM_CSOFT_CLIP 	4
H 	BAM_CHARD_CLIP 	5
P 	BAM_CPAD 	6
= 	BAM_CEQUAL 	7
X 	BAM_CDIFF 	8
"""
CIGAR_DICT = {0:{"seq_add": 1, "ref_add": 1}, # M
              1:{"seq_add": 1, "ref_add": 0}, # I
              2:{"seq_add": 0, "ref_add": 1}, # D
              3:{"seq_add": 0, "ref_add": 1}, # N
              4:{"seq_add": 1, "ref_add": 0}, # S
              5:{"seq_add": 0, "ref_add": 0}, # H
              6:{"seq_add": 0, "ref_add": 1}, # P
              7:{"seq_add": 1, "ref_add": 1}, # =
              8:{"seq_add": 1, "ref_add": 1}} # X

def main(argv):
    """Reads SAM format from stdin or SAM/BAM from first argument.
    Parses any reads with splice junctions (N in CIGAR) and splits them into
    multiple reads. 0x40 and 0x80 is set in FLAG for the split fragments to 
    indicate "chimeric" alignments. Output is uncompressed SAM format that 
    needs parsing for read group (RG:A may be produced instead of RG:Z), 
    sorting and indexing. 
    
    Usage:
    samtools view -h my.bam | SplitNReads.py | sed 's/RG:A:/RG:Z:/g' | samtools calmd -Su - /path/to/genome.fa | samtools sort - out
    splitNReads.py my.bam | sed 's/RG:A:/RG:Z:/g' | samtools calmd -Su - /path/to/genome.fa | samtools sort - out
    """
    # default read mode: sam
    read_mode = "r" # by default input is sam
    write_mode = "wh" # wbu" # write SAM
    out_filename = "-" # write to stdout
    if len(argv) < 2:
        # read from stdin
        in_filename = "-"
    else:
        # read from a file
        in_filename = argv[1]
        # check if input is binary (in which case output is too)
        if in_filename.endswith("bam"):
            read_mode = "rb"
    
    infile = pysam.Samfile( in_filename, read_mode )
    outfile = pysam.Samfile( out_filename, write_mode, template = infile )
    for s in infile: 
        if 0x4&s.flag or not "N" in s.cigarstring: # unaligned or no splice junctions
            changeMAPQ(my_read = s, new_qual = 60)
            outfile.write(s)
            continue
        # break CIGAR into N-divided blocks        
        cigar_array = breakdown(s.cigar, s.aend-s.pos)
        is_first = True
        for my_op in cigar_array:
            a = pysam.AlignedRead()
            copy_defaults(new = a, old = s) # copy fields that are not affected
            if is_first:
                is_first = False
                #a.flag = a.flag|0x40|0x80 # means: not first and not last segment
            else:
                #a.flag = a.flag|0x40|0x80|0x800 # means: not first and not last segment + supplementary
                a.flag = a.flag|0x800 # means: supplementary, keep other flags intact
            
            # base start location in current read fragment
            start_pos = my_op["seq_start"]
            # offset to base start location in current read
            offset = my_op["seq_offset"]
            # fix start position
            a.pos = s.pos + my_op["ref_offset"]
            # extract bases
            a.seq = s.seq[start_pos:start_pos+offset]
            # extract quality values
            a.qual = s.qual[start_pos:start_pos+offset]
            # update CIGAR
            a.cigar = my_op["cigar"]
            # fix relative start/end position of mate
            if s.tlen > 0:
                a.tlen = a.tlen - my_op["ref_offset"] # dist from this reads beginning to the mate's end
            elif s.tlen < 0 and s.pnext > 0:
                a.tlen = s.pnext - a.aend # dist from mate's beginning to this read's end
            changeMAPQ(my_read = a, new_qual = 60)
            outfile.write(a)
        
def changeMAPQ(my_read, new_qual=60):
    if my_read.mapq == 255: 
        my_read.mapq = new_qual # Gatk best practices..

def copy_defaults(new, old):
    new.qname = old.qname
    #new.seq  = "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
    new.flag = old.flag
    new.rname = old.rname
    #new.pos = old.pos
    new.mapq = old.mapq
    #new.cigar = ( (0,10), (2,1), (0,25) )
    #new.mrnm = 0
    new.rnext = old.rnext
    #new.mpos=199
    new.pnext = old.pnext
    #new.isize = 167
    new.tlen = old.tlen
    #new.qual="<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<"
    new.tags = old.tags[:]
    #import pdb; pdb.set_trace()

def breakdown(cigar, orig_astart_end_diff):
    # break down cigar string based on splice junctions
    cigar_array = []
    ref_offset = 0 # offset to reference start position
    seq_offset = 0 # read start base position offset
    temp_dict = {"ref_offset": ref_offset, "seq_start": 0, "cigar": tuple()}
    for c in cigar:
        if c[0] == 3: # is N?
            if check_cigar_has_content(temp_dict.get("cigar", ())):
                temp_dict["seq_offset"] = seq_offset
                temp_dict["cigar"] += ((5, orig_astart_end_diff-ref_offset),) # add hard clipping as per Gatk best practices
                cigar_array.append(temp_dict)
            ref_offset += c[1]
            temp_dict = {"ref_offset": ref_offset, "seq_start": seq_offset, 
                         "cigar": ((5, ref_offset),)} # add some hard clipping as per Gatk best practices
            continue
        # add to offsets depending on whether CIGAR operator increments one of
        # reference position or read base position
        seq_offset += CIGAR_DICT[c[0]]["seq_add"] * c[1]
        ref_offset += CIGAR_DICT[c[0]]["ref_add"] * c[1]
        temp_dict["cigar"] += (c,)
    if check_cigar_has_content(temp_dict.get("cigar", ())):
        temp_dict["seq_offset"] = seq_offset
        cigar_array.append(temp_dict)
    return cigar_array

def check_cigar_has_content(cigar):
    if len(cigar) < 1:
        return False
    for c in cigar:
        if c[0] in (0, 7, 8) and c[1] > 0: # is M, = or X
            return True
    return False
if __name__ == "__main__":
    main(sys.argv)
