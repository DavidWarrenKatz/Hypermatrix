import pysam
import argparse


def mh_hic_filter(in_cram, out_summary):
    uniq_mapped_mapq30_nopcr_cis_1kb = 0
    uniq_mapped_mapq30_nopcr_cis_10kb = 0
    uniq_mapped_mapq30_nopcr_cis_20kb = 0
    long_cis_contact_per_chr_filter_flag = 0  #mark as 1 if in any chr, total long-range cis contacts of each chromosome <= x, for a chromosome of length × Mb
    pre_read=None
    sam_chrom_dict=dict()
    cis_chrom_dict=dict()
    with pysam.AlignmentFile(in_cram, "rc") as samfile:
        sam_header_dict = samfile.header.to_dict()['SQ']
        for item in sam_header_dict:
            chrom=item['SN']
            val = int(item['LN']/1000000)
            if (not chrom.startswith('chrUn')) and (not chrom.endswith('random')) and (not chrom=='chrM') and (not chrom=='chrX') and (not chrom=='chrY'):
                sam_chrom_dict[chrom]=val
                cis_chrom_dict[chrom]=0
        for read in samfile.fetch(until_eof=True):
            if pre_read is None:
                pre_read = read
                continue
            elif pre_read.query_name != read.query_name:
                pre_read = read
                continue
            elif pre_read.query_name == read.query_name:
                if read.is_paired and (not read.is_unmapped) and (not pre_read.is_unmapped) and (not read.is_secondary) and (not pre_read.is_secondary) \
                        and read.mapping_quality >= 30 and pre_read.mapping_quality >= 30 and (not read.is_duplicate) \
                    and read.reference_name == read.next_reference_name:
                    if abs(read.template_length) >= 1000:
                        uniq_mapped_mapq30_nopcr_cis_1kb = uniq_mapped_mapq30_nopcr_cis_1kb + 1
                        if abs(read.template_length) >= 10000:
                            uniq_mapped_mapq30_nopcr_cis_10kb = uniq_mapped_mapq30_nopcr_cis_10kb + 1
                            chrom = read.reference_name
                            if (not chrom.startswith('chrUn')) and (not chrom.endswith('random')) and (not chrom=='chrM') and (not chrom=='chrX') and (not chrom=='chrY'):
                                cis_chrom_dict[chrom] = cis_chrom_dict[chrom] + 1
                            if abs(read.template_length) >= 20000:
                                uniq_mapped_mapq30_nopcr_cis_20kb = uniq_mapped_mapq30_nopcr_cis_20kb + 1

    for chrom in sam_chrom_dict.keys():
        chr_len = sam_chrom_dict[chrom]
        frag_num = cis_chrom_dict[chrom]
        if frag_num <= chr_len:
            long_cis_contact_per_chr_filter_flag = 1
            break

    with open(out_summary, 'w') as out:
        line = "UniqMappedMapQ30NoPcrCisMore1kb:\t" + str(uniq_mapped_mapq30_nopcr_cis_1kb) \
               + "\nUniqMappedMapQ30NoPcrCisMore10kb:\t" + str(uniq_mapped_mapq30_nopcr_cis_10kb) \
               + "\nUniqMappedMapQ30NoPcrCisMore20kb:\t" + str(uniq_mapped_mapq30_nopcr_cis_20kb) \
               + "\nPerChromCisLongFilter:\t" + str(long_cis_contact_per_chr_filter_flag) + "\n"
        #for chrom in sam_chrom_dict.keys():
        #    frag_num = cis_chrom_dict[chrom]
        #    line = line + chrom + "\t" + str(sam_chrom_dict[chrom]) + "\t" + str(frag_num) + "\n"
        out.write(line)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--in_cram') ##required sorted by read_name. R1 and R2's read name should be the same.  not with /1 or /2 suffix.
    p.add_argument('--out_summary') ##only check pair-end reads. only use unique mapped reads with mapQ>=30, then count PCR dups, cis-, trans-,
    args = p.parse_args()
    mh_hic_filter(**vars(args))
