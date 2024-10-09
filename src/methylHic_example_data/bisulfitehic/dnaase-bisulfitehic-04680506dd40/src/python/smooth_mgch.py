import argparse, gzip
from bx.intervals.intersection import Intersecter, Interval
from scipy import interpolate
import numpy as np
import pysam
from scipy.interpolate import UnivariateSpline


def smooth_mgch_per_chrom(tbx, out, chrom, size, window, max_window, min_gch, column, s):
    column = column - 1
    intersecter = Intersecter()
    num = 0
    for row in tbx.fetch(chrom, 0, size):
        splitline = row.split("\t")
        chrom = splitline[0]
        start = int(splitline[1])
        end = int(splitline[2])
        v = float(splitline[column])
        interv_tmp = Interval(start, end, chrom=chrom,
                              value={'methylation': v})
        intersecter.add_interval(interv_tmp)
        num = num + 1
    print('with size of GCH:' + str(num))

    for line in tbx.fetch(chrom, 0, size):
        splitline = line.split("\t")
        chrom = splitline[0]
        start = int(splitline[1])
        end = int(splitline[2])
        v = float(splitline[column])
        sites_after = intersecter.after(start)
        if sites_after is not None and len(sites_after) > 0:
            site_after = sites_after[0]
            num_gch = 0
            methy_list = []
            range_list = []
            w_start_1 = max(0, start - window)
            w_end_1 = end
            w_start_2 = site_after.start
            w_end_2 = min(size, site_after.start + window)
            while (num_gch < min_gch) and (
                    abs(w_start_1 - int(start)) <= max_window and abs(w_end_2 - int(site_after.end)) <= max_window):
                for interval in intersecter.find(w_start_1, w_end_1):
                    s = interval.start
                    m = interval.value['methylation']
                    methy_list.append(m)
                    range_list.append(s)
                    num_gch = num_gch + 1
                for interval in intersecter.find(w_start_2, w_end_2):
                    s = interval.start
                    m = interval.value['methylation']
                    methy_list.append(m)
                    range_list.append(s)
                    num_gch = num_gch + 1
                w_end_1 = w_start_1
                w_start_1 = max(0, w_start_1 - window)
                w_start_2 = w_end_2
                w_end_2 = min(size, int(w_end_2 - 1) + window)
            methy_list = np.array(methy_list, dtype='float32')
            range_list = np.array(range_list, dtype='int32')
            # f = interpolate.interp1d(range_list, methy_list, kind=kind)
            sort_index = np.argsort(range_list)
            #print(sort_index)
            range_list = range_list[sort_index]
            methy_list = methy_list[sort_index]
            f = UnivariateSpline(range_list, methy_list, k=1, s=s)
            xn = range(start, site_after.start)
            yn = f(xn)
            for i in range(len(xn)):
                i_start = xn[i]
                i_end = i_start + 1
                i_y = max(0.0, yn[i])
                out.write(chrom + '\t' + str(i_start) + '\t' + str(i_end) + '\t' + str(round(i_y, 4)) + '\n')
        else:
            out.write(chrom + '\t' + str(start) + '\t' + str(end) + '\t' + str(round(v, 4)) + '\n')


def smooth_mgch(in_bed, out_bed, chrom_size, window, max_window, min_gch, column, s):
    chrom_dict = dict()
    with open(chrom_size, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            splitline = line.split("\t")
            chrom = splitline[0]
            size = int(splitline[1])
            chrom_dict[chrom] = size
    tbx = pysam.TabixFile(in_bed)
    if out_bed.endswith('.gz'):
        out = gzip.open(out_bed, 'wt')
    else:
        out = open(out_bed, 'w')
    for chrom in chrom_dict.keys():
        if chrom in tbx.contigs:
            print('chrom:' + chrom)
            size = chrom_dict[chrom]
            smooth_mgch_per_chrom(tbx, out, chrom, size, window, max_window, min_gch, column, s)
    out.close()


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--in_bed')  ##bgzip and indexed 6plus2 bed file, standard bed + methy (%, 0-100 scale) + coverage
    p.add_argument(
        '--out_bed')
    p.add_argument('--chrom_size')
    p.add_argument('--window', type=int, default=1000)  ##left 1000bp and right 1000bp
    p.add_argument('--max_window', type=int, default=10000)
    p.add_argument('--min_gch', type=int, default=100)  ##if less than 3 GCH, extend another 100bp, until +/-1000bp
    p.add_argument('--column', type=int, default=7)  ##the column number in the file to do the smoothing
    p.add_argument('--s', type=float, default=1)
    args = p.parse_args()
    smooth_mgch(**vars(
        args))  ##1d smooth GCH with adjacent +/-1kb window with at least 100 GCH
