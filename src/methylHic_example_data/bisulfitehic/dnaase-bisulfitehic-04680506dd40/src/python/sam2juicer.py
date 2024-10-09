#!/usr/bin/python

#########################################
# Author: Yunjiang Qiu <serein927@gmail.com>
# File: bam2juice.py
# Create Date: 2016-06-22 18:43:11
#########################################

import sys
import argparse
import pysam
from _bisect import *

def main():
    parser = argparse.ArgumentParser(description='transfrom sam file to juicer format')
    parser.add_argument("-s", "--sam", dest="sam",required=True,help="input sam file")
    parser.add_argument("-f", "--frag", dest="frag",required=True,help="fragment site file")
    args = parser.parse_args()

    frag = {}
    with open(args.frag, 'r') as f:
        for line in f:
            tmp = line.rstrip().split(' ')
            frag[tmp[0].replace("chr","").replace("MT","M")] = [int(x) for x in tmp[1:]]

    bam = pysam.Samfile(args.sam,"r")
    reads = {}
    end_of_file = False
    read = bam.next()

    while not end_of_file:
        if read.qname in reads:
            read1 = read
            read2 = reads.pop(read.qname)

            try:
                frag1 = bisect(frag[bam.getrname(read1.tid).replace("chr","")], read1.pos)
                frag2 = bisect(frag[bam.getrname(read2.tid).replace("chr","")], read2.pos)
            except:
                sys.exit('Something wrong locating fragments for %s or %s.\n' % (bam.getrname(read1.tid).replace("chr",""), bam.getrname(read2.tid).replace("chr","")))

            first = 0
            if read1.tid > read2.tid:
                first = 1
            elif read1.tid == read2.tid:
                 if frag1 > frag2:
                     first = 1
                 elif frag1 == frag2:
                     if (read1.flag & 16) > (read2.flag & 16):
                         first = 1
                     elif (read1.flag & 16) == (read2.flag & 16):
                         if read1.pos > read2.pos:
                             first = 1

            if first == 1:
                print("%s" % ' '.join([str(i) for i in [(read2.flag & 16), bam.getrname(read2.tid), read2.pos, frag2, (read1.flag & 16), bam.getrname(read1.tid), read1.pos, frag1, read2.mapq, read2.cigarstring, read2.seq, read1.mapq, read1.cigarstring, read1.seq, read2.qname, read1.qname] ]))
            else:
                print("%s" % ' '.join([str(i) for i in [(read1.flag & 16), bam.getrname(read1.tid), read1.pos, frag1, (read2.flag & 16), bam.getrname(read2.tid), read2.pos, frag2, read1.mapq, read1.cigarstring, read1.seq, read2.mapq, read2.cigarstring, read2.seq, read1.qname, read2.qname] ]))

        else:
            reads[read.qname] = read

        try:
            read = bam.next()
        except:
            end_of_file = True

if __name__ == "__main__":
    sys.exit(main())
