import argparse
from bx.intervals.intersection import Intersecter, Interval
from scipy import stats


def binom_test_mgch(in_bed, out_bed, min_gch, min_point, window, max_window, global_acc):
    gch_intervals_dict = dict()
    num = 0
    with open(in_bed, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('track') or line.startswith('#'):
                continue
            chrom, start, end, _, score, strand, m, c = line.split("\t")
            interv_tmp = Interval(int(start), int(end), chrom=chrom, strand=strand,
                                  value={'score': int(score), 'methylation': float(m), 'coverage': int(c)})
            if chrom in gch_intervals_dict.keys():
                intersecter = gch_intervals_dict.get(chrom)
            else:
                intersecter = Intersecter()
            intersecter.add_interval(interv_tmp)
            gch_intervals_dict[chrom] = intersecter
            num = num + 1
    print('build interval tree with size:' + str(num))

    with open(out_bed, 'w') as out:
        with open(in_bed, 'r') as f:
            for line in f:
                line = line.rstrip('\n')
                if line.startswith('track') or line.startswith('#'):
                    continue
                chrom, start, end, _, score, strand, m, c = line.split("\t")
                methy_count = int(float(m) * int(c) / 100.0)
                intersecter = gch_intervals_dict.get(chrom)
                num_gch = 0
                c_window = 0
                m_window = 0
                w_start_1 = int(start) - window
                w_end_1 = int(start)
                w_start_2 = int(end)
                w_end_2 = int(end) + window
                while (num_gch < min_gch or c_window < min_point) and abs(w_start_1 - int(start)) <= max_window:
                    for interval in intersecter.find(w_start_1, w_end_1):
                        c_int = interval.value['coverage']
                        c_window = c_window + c_int
                        m_window = m_window + int(interval.value['methylation'] * c_int / 100.0)
                        num_gch = num_gch + 1
                    for interval in intersecter.find(w_start_2, w_end_2):
                        c_int = interval.value['coverage']
                        c_window = c_window + c_int
                        m_window = m_window + int(interval.value['methylation'] * c_int / 100.0)
                        num_gch = num_gch + 1
                    w_start_1 = w_start_1 - window
                    w_end_1 = w_start_1 + window
                    w_start_2 = w_start_2 + window
                    w_end_2 = w_start_2 + window
                if c_window > 0:
                    window_p = m_window / c_window
                else:
                    window_p = global_acc
                p_value = stats.binom_test(methy_count, n=int(c), p=window_p, alternative='greater')
                out.write(line + '\t' + str(round(p_value, 2)) + '\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--in_bed')
    p.add_argument('--out_bed')
    p.add_argument('--window', type=int, default=100)  ##left 100bp and right 100bp
    p.add_argument('--min_gch', type=int, default=3)  ##if less than 3 GCH, extend another 100bp, until +/-1000bp
    p.add_argument('--min_point', type=int,
                   default=10)  ##if less than 10 data point, extend another 100bp, until +/-1000bp
    p.add_argument('--max_window', type=int, default=1000)
    p.add_argument('--global_acc', type=float, default=0.285)
    args = p.parse_args()
    binom_test_mgch(**vars(args))
