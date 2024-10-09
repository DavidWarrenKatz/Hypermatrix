import pysam
import argparse, re
from bx.intervals.intersection import Intersecter, Interval


def pass_filter(read, min_mapq, min_insert_size):
    return (not (read.is_unmapped or read.mate_is_unmapped or read.is_secondary
                 or read.is_duplicate or read.is_qcfail or read.is_supplementary)) \
           and read.is_paired and read.mapping_quality >= min_mapq \
           and (abs(read.template_length)>=min_insert_size or abs(read.template_length)==0)


def split_regions(regions):
    chrom, loc = re.split(':', regions)
    start, end = re.split('-', loc)
    return [chrom, (int)(start), (int)(end)]


def interact_density_around_regions(in_bam, out_summary, from_region, plot_center, plot_bin_size, plot_boundary, min_mapq, min_insert_size):
    from_chr, from_start, from_end = split_regions(from_region)
    plot_chr, plot_start, plot_end = split_regions(plot_center)
    plot_cent_loc = (int)((plot_start + plot_end) / 2)
    to_intervals = Intersecter()
    to_intervals.add_interval(Interval(int(plot_cent_loc-plot_boundary), int(plot_cent_loc+plot_boundary), chrom=plot_chr))
    num_bins = (int)(2*plot_boundary/plot_bin_size)
    target_bins=[0]*num_bins
    with pysam.AlignmentFile(in_bam, "rb") as samfile:
        for read in samfile.fetch(contig=from_chr, start=from_start, stop=from_end):
            if pass_filter(read, min_mapq, min_insert_size):
                if read.next_reference_name == plot_chr:
                    mate_start = read.next_reference_start
                    if to_intervals.find(mate_start, mate_start+1):
                        bin_num_mate = (int)(num_bins/2) + (int)((mate_start-plot_cent_loc)/plot_bin_size)
                        target_bins[bin_num_mate]+=1
    print(sum(target_bins))
    with open(out_summary, 'w') as out:
        i = 0-(int)(num_bins/2)
        for s in target_bins:
            out.write(str(i) + '\t')
            i = i+1
        out.write('\n')
        for s in target_bins:
            out.write(str(s) + '\t')
        out.write('\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--in_bam')
    p.add_argument('--out_summary') ##header: bin location; 2nd row: the density of reads
    p.add_argument('--from_region', default='chr17:85611135-85828994')
    p.add_argument('--plot_center', default='chr1:115724771-115798945')
    p.add_argument('--plot_bin_size', type=int, default=50000)
    p.add_argument('--plot_boundary', type=int, default=5000000)
    p.add_argument('--min_mapq', type=int, default=30)
    p.add_argument('--min_insert_size', type=int, default=20000)
    args = p.parse_args()
    interact_density_around_regions(**vars(args))
