#!/usr/bin/python

#########################################
# Author: Yunjiang Qiu <serein927@gmail.com>
# File: bam2juice.py
# Create Date: 2016-06-22 18:43:11
# updated by David Katz Feb 28, 2024 <davidkatz02@gmail.com>
#########################################

import sys
import argparse
import pysam
from bisect import bisect

def main():
    parser = argparse.ArgumentParser(description='transform sam file to juicer format')
    parser.add_argument("-s", "--sam", dest="sam", required=True, help="input sam file")
    parser.add_argument("-f", "--frag", dest="frag", required=True, help="fragment site file")
    args = parser.parse_args()

    frag = {}
    with open(args.frag, 'r') as f:
        for line in f:
            tmp = line.rstrip().split(' ')
            frag[tmp[0].replace("chr", "").replace("MT", "M")] = [int(x) for x in tmp[1:]]

    bam = pysam.AlignmentFile(args.sam, "r")
    reads = {}
    
    for read in bam:
        if read.query_name in reads:
            read1 = read
            read2 = reads.pop(read.query_name)

            try:
                frag1 = bisect(frag[bam.get_reference_name(read1.reference_id).replace("chr", "")], read1.reference_start)
                frag2 = bisect(frag[bam.get_reference_name(read2.reference_id).replace("chr", "")], read2.reference_start)
            except:
                sys.exit(f'Something wrong locating fragments for {bam.get_reference_name(read1.reference_id).replace("chr", "")} or {bam.get_reference_name(read2.reference_id).replace("chr", "")}.\n')

            first = 0
            if read1.reference_id > read2.reference_id:
                first = 1
            elif read1.reference_id == read2.reference_id:
                if frag1 > frag2:
                    first = 1
                elif frag1 == frag2:
                    if (read1.flag & 16) > (read2.flag & 16):
                        first = 1
                    elif (read1.flag & 16) == (read2.flag & 16):
                        if read1.reference_start > read2.reference_start:
                            first = 1

            if first == 1:
                print("%s" % ' '.join([str(i) for i in [(read2.flag & 16), bam.get_reference_name(read2.reference_id), read2.reference_start, frag2, (read1.flag & 16), bam.get_reference_name(read1.reference_id), read1.reference_start, frag1, read2.mapping_quality, read2.cigarstring, read2.query_sequence, read1.mapping_quality, read1.cigarstring, read1.query_sequence, read2.query_name, read1.query_name]]))
            else:
                print("%s" % ' '.join([str(i) for i in [(read1.flag & 16), bam.get_reference_name(read1.reference_id), read1.reference_start, frag1, (read2.flag & 16), bam.get_reference_name(read2.reference_id), read2.reference_start, frag2, read1.mapping_quality, read1.cigarstring, read1.query_sequence, read2.mapping_quality, read2.cigarstring, read2.query_sequence, read1.query_name, read2.query_name]]))

        else:
            reads[read.query_name] = read

if __name__ == "__main__":
    sys.exit(main())



