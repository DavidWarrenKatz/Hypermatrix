import argparse, gzip
from bx.intervals.intersection import Intersecter, Interval
from scipy import stats
import numpy as np
import pysam
from fastbetabino import *

def fit_betabinomial(methy_list, cov_list): # a should be an array of the number of methylated point in each GCH
    alpha, beta = fit_alpha_beta(cov_list, methy_list)
    return (alpha, beta)


def betabinom_test_mgch_per_chrom(tbx, out, chrom, size, min_gch, min_point, window1, window2, binomial_only, two_sided):
    intersecter = Intersecter()
    num = 0
    c_global = 0
    m_global = 0
    global_methy = []
    global_m = []
    global_cov = []
    for row in tbx.fetch(chrom, 0, size):
        chrom, start, end, _, score, strand, m, c = row.split("\t")
        interv_tmp = Interval(int(start), int(end), chrom=chrom, strand=strand,
                              value={'score': int(score), 'methylation': float(m), 'coverage': int(c)})
        intersecter.add_interval(interv_tmp)
        c_global = c_global + int(c)
        m_global = m_global + int(float(m) * int(c) / 100.0)
        global_methy.append(int(float(m) * int(c) / 100.0))
        global_cov.append(int(c))
        global_m.append(float(m)/100.0)
        num = num + 1
    print('with size of GCH:' + str(num))
    global_methy = np.array(global_methy, dtype='int32')
    global_cov = np.array(global_cov, dtype='int32')
    global_m = np.array(global_m, dtype='float32')
    global_mean = np.nanmean(global_m)
    global_std = np.nanstd(global_m)
    if num < min_gch:
        print('too small number of GCH in this chrom, skip')
        return
    if not binomial_only:
        if num > 1000000:
            print('to speed up, random sampling 1000000 GCH in this chrom to estimate the distribution')
            s_index = np.random.choice(global_methy.shape[0], 1000000, replace=False)
            global_methy_s = global_methy[s_index]
            global_cov_s = global_cov[s_index]
        else:
            global_methy_s = global_methy
            global_cov_s = global_cov
        global_methy_alpha, global_methy_beta = fit_betabinomial(global_methy_s, global_cov_s)
        print('global alpha:' + str(global_methy_alpha))
        print('global beta:' + str(global_methy_beta))
        del global_methy_s
        del global_cov_s

    del global_methy
    del global_cov
    del global_m


    for line in tbx.fetch(chrom, 0, size):
        chrom, start, end, _, score, strand, m, c = line.split("\t")
        methy_count = int(float(m) * int(c) / 100.0)
        total_count = int(c)
        methy = float(m) / 100
        num_gch_1 = 0
        c_window_1 = 0
        m_window_1 = 0
        meth_list1 = []
        m_list1 = []
        cov_list1 = []
        w1_start_1 = max(int(start) - window1, 0)
        w1_end_1 = int(start)
        w1_start_2 = int(end)
        w1_end_2 = min(int(end) + window1, size)
        for interval in intersecter.find(w1_start_1, w1_end_1):
            c_int = interval.value['coverage']
            c_window_1 = c_window_1 + c_int
            m_window_1 = m_window_1 + int(interval.value['methylation'] * c_int / 100.0)
            meth_list1.append(int(interval.value['methylation'] * c_int / 100.0))
            cov_list1.append(c_int)
            m_list1.append(float(interval.value['methylation'])/100.0)
            num_gch_1 = num_gch_1 + 1
        for interval in intersecter.find(w1_start_2, w1_end_2):
            c_int = interval.value['coverage']
            c_window_1 = c_window_1 + c_int
            m_window_1 = m_window_1 + int(interval.value['methylation'] * c_int / 100.0)
            meth_list1.append(int(interval.value['methylation'] * c_int / 100.0))
            cov_list1.append(c_int)
            m_list1.append(float(interval.value['methylation']) / 100.0)
            num_gch_1 = num_gch_1 + 1
        meth_list1 = np.array(meth_list1, dtype='int32')
        cov_list1 = np.array(cov_list1, dtype='int32')
        m_list1 = np.array(m_list1, dtype='float32')
        local_mean1 = np.nanmean(m_list1)
        local_std1 = np.nanstd(m_list1)
        #if not binomial_only:
            #local_methy1_alpha, local_methy1_beta = fit_betabinomial(meth_list1, cov_list1)
        if c_window_1 > 0:
            local_fold1 = methy - (m_window_1 / c_window_1)
            zscore_w1 = (methy - local_mean1) / local_std1
        else:
            local_fold1 = methy - m_global/c_global
            zscore_w1 = (methy - global_mean) / global_std
        #print(local_methy1_alpha, local_methy1_beta)
        #print(meth_list1.shape, cov_list1.shape)
        num_gch_2 = 0
        c_window_2 = 0
        m_window_2 = 0
        meth_list2 = []
        m_list2 = []
        cov_list2 = []
        w2_start_1 = max(int(start) - window2, 0)
        w2_end_1 = int(start)
        w2_start_2 = int(end)
        w2_end_2 = min(int(end) + window2, size)
        for interval in intersecter.find(w2_start_1, w2_end_1):
            c_int = interval.value['coverage']
            c_window_2 = c_window_2 + c_int
            m_window_2 = m_window_2 + int(interval.value['methylation'] * c_int / 100.0)
            meth_list2.append(int(interval.value['methylation'] * c_int / 100.0))
            cov_list2.append(c_int)
            m_list2.append(float(interval.value['methylation']) / 100.0)
            num_gch_2 = num_gch_2 + 1
        for interval in intersecter.find(w2_start_2, w2_end_2):
            c_int = interval.value['coverage']
            c_window_2 = c_window_2 + c_int
            m_window_2 = m_window_2 + int(interval.value['methylation'] * c_int / 100.0)
            meth_list2.append(int(interval.value['methylation'] * c_int / 100.0))
            cov_list2.append(c_int)
            m_list2.append(float(interval.value['methylation']) / 100.0)
            num_gch_2 = num_gch_2 + 1
        meth_list2 = np.array(meth_list2, dtype='int32')
        cov_list2 = np.array(cov_list2, dtype='int32')
        m_list2 = np.array(m_list2, dtype='float32')
        local_mean2 = np.nanmean(m_list2)
        local_std2 = np.nanstd(m_list2)
        #if not binomial_only:
           # local_methy2_alpha, local_methy2_beta = fit_betabinomial(meth_list2, cov_list2)
        #print(local_methy2_alpha, local_methy2_beta)
        #print(meth_list2.shape, cov_list2.shape)
        if c_window_2 > 0:
            local_fold2 = methy - (m_window_2 / c_window_2)
            zscore_w2 = (methy - local_mean2) / local_std2
        else:
            local_fold2 = methy - (m_global / c_global)
            zscore_w2 = (methy - global_mean) / global_std
        global_fold = methy - (m_global / c_global)
        zscore_global = (methy - global_mean) / global_std

        if num_gch_1 >= min_gch and c_window_1 >= min_point:
            window1_p = m_window_1 / c_window_1
            #if not binomial_only:
                #beta_binom_p_w1 = 0 - np.log10(stats.betabinom.cdf(methy_count, total_count, local_methy1_alpha, local_methy1_beta))
        else:
            window1_p = m_global / c_global
            #if not binomial_only:
                #beta_binom_p_w1 = 0 - np.log10(stats.betabinom.cdf(methy_count, total_count, global_methy_alpha, global_methy_beta))
        if two_sided:
            binomial_p_w1 = 0 - np.log10(stats.binom_test(methy_count, n=total_count, p=window1_p, alternative='two-sided'))
        else:
            binomial_p_w1 = 0 - np.log10(stats.binom_test(methy_count, n=total_count, p=window1_p, alternative='less'))
        if num_gch_2 >= min_gch and c_window_2 >= min_point :
            window2_p = m_window_2 / c_window_2
            #if not binomial_only:
                #beta_binom_p_w2 = 0 - np.log10(stats.betabinom.cdf(methy_count, total_count, local_methy2_alpha, local_methy2_beta))
        else:
            window2_p = m_global / c_global
            #if not binomial_only:
                #beta_binom_p_w2 = 0 - np.log10(stats.betabinom.cdf(methy_count, total_count, global_methy_alpha, global_methy_beta))
        if two_sided:
            binomial_p_w2 = 0 - np.log(stats.binom_test(methy_count, n=total_count, p=window2_p, alternative='two-sided'))
        else:
            binomial_p_w2 = 0 - np.log(stats.binom_test(methy_count, n=total_count, p=window2_p, alternative='less'))
        if two_sided:
            binomial_p_global = 0 - np.log10(
                stats.binom_test(methy_count, n=total_count, p=m_global / c_global, alternative='two-sided'))
        else:
            binomial_p_global = 0 - np.log10(
                stats.binom_test(methy_count, n=total_count, p=m_global / c_global, alternative='less'))
        binomial_p_final = np.nanmax((binomial_p_w1, binomial_p_w2, binomial_p_global))


        if not binomial_only:
            beta_binom_p_global = 0 - np.log10(stats.betabinom.cdf(methy_count, total_count, global_methy_alpha, global_methy_beta))
            #beta_binom_p_final = np.nanmax((beta_binom_p_w1, beta_binom_p_w2, beta_binom_p_global))
            out.write(line + '\t' + str(round(local_fold1, 4)) + '\t' + str(round(local_fold2, 4)) + '\t' + str(
                round(global_fold, 4)) + '\t' + str(round(binomial_p_w1, 4)) + '\t' + str(
                round(binomial_p_w2, 4)) + '\t' + str(round(binomial_p_global, 4)) + '\t' + str(
                binomial_p_final) + '\t' + str(round(beta_binom_p_global, 4)) + '\t' + str(round(zscore_w1, 4)) + '\t' + str(round(zscore_w2, 4)) + '\t' + str(round(zscore_global, 4)) + '\n')
        else:
            out.write(line + '\t' + str(round(local_fold1, 4)) + '\t' + str(round(local_fold2, 4)) + '\t' + str(
                round(global_fold, 4)) + '\t' + str(round(binomial_p_w1, 4)) + '\t' + str(
                round(binomial_p_w2, 4)) + '\t' + str(round(binomial_p_global, 4)) + '\t' + str(
                binomial_p_final) + '\t' + str(round(zscore_w1, 4)) + '\t' + str(round(zscore_w2, 4)) + '\t' + str(round(zscore_global, 4)) + '\n')


def betabinom_test_mgch(in_bed, out_bed, chrom_size, min_gch, min_point, window1, window2, binomial_only, two_sided):
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
            betabinom_test_mgch_per_chrom(tbx, out, chrom, size, min_gch, min_point, window1, window2, binomial_only, two_sided)
    out.close()


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--in_bed')  ##bgzip and indexed 6plus2 bed file, standard bed + methy (%, 0-100 scale) + coverage
    p.add_argument(
        '--out_bed')  ## standard input bed, then log2 fold change (local bg 1, local bg 2, global bg, min), -log10p (local bg 1, local bg 2, global bg, max), z-score(local bg1, local bg2, global bg)
    p.add_argument('--chrom_size')  ## to get accessibility background in each chr; recommend only autosome + chrX
    p.add_argument('--window1', type=int, default=2500)  ##left 2500bp and right 2500bp
    p.add_argument('--window2', type=int, default=5000)  ##left 5000bp and right 5000bp
    p.add_argument('--min_gch', type=int, default=10)  ##requires at least 10 GCH to fit a good beta-binomial distribution
    p.add_argument('--min_point', type=int,
                   default=100)  ##requires at least 100 data points
    p.add_argument('--binomial_only', type=bool, default=False)
    p.add_argument('--two_sided', type=bool, default=False)
    args = p.parse_args()
    betabinom_test_mgch(**vars(
        args))  ##get beta-binomial test, -log10 pvalue, for regions greater methy, "greater", regions smaller methy, use "smaller"hypothesis testing. # maybe also get FDR value
