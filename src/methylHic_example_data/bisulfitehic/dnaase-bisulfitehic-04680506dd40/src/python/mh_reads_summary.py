import pysam
import argparse


def mh_reads_summary(in_cram, out_summary):
    total_reads = 0
    total_mapped = 0
    uniq_mapped = 0
    uniq_reads_to_lambda = 0
    uniq_pairs_to_lambda = 0
    uniq_mapped_to_target_species = 0
    uniq_paired_to_target_species = 0
    uniq_mapped_nopcr = 0
    uniq_mapped_nopcr_mapq1 = 0
    uniq_mapped_nopcr_mapq1_trans = 0
    uniq_mapped_nopcr_mapq1_cis = 0
    uniq_mapped_nopcr_mapq1_cis_1kb = 0
    uniq_mapped_nopcr_mapq1_cis_20kb = 0
    uniq_mapped_nopcr_mapq30 = 0
    uniq_mapped_nopcr_mapq30_trans = 0
    uniq_mapped_nopcr_mapq30_cis = 0
    uniq_mapped_nopcr_mapq30_cis_1kb = 0
    uniq_mapped_nopcr_mapq30_cis_20kb = 0
    with pysam.AlignmentFile(in_cram, "rc") as samfile:
        for read in samfile.fetch(until_eof=True):
            total_reads = total_reads + 1
            if not read.is_unmapped:
                total_mapped = total_mapped + 1
                if not read.is_secondary:
                    uniq_mapped = uniq_mapped + 1
                    if read.reference_name == "lambda":
                        uniq_reads_to_lambda = uniq_reads_to_lambda + 1
                        if read.is_paired and (not read.mate_is_unmapped) and read.next_reference_name == "lambda":
                            uniq_pairs_to_lambda = uniq_pairs_to_lambda + 1
                    else:
                        uniq_mapped_to_target_species = uniq_mapped_to_target_species + 1
                        if read.is_paired and (not read.mate_is_unmapped):
                            uniq_paired_to_target_species = uniq_paired_to_target_species + 1
                            if not read.is_duplicate:
                                uniq_mapped_nopcr = uniq_mapped_nopcr + 1
                                if read.mapping_quality >= 1:
                                    uniq_mapped_nopcr_mapq1 = uniq_mapped_nopcr_mapq1 + 1
                                    if read.reference_name != read.next_reference_name:
                                        uniq_mapped_nopcr_mapq1_trans = uniq_mapped_nopcr_mapq1_trans + 1
                                    else:
                                        uniq_mapped_nopcr_mapq1_cis = uniq_mapped_nopcr_mapq1_cis + 1
                                        if abs(read.template_length) >= 1000:
                                            uniq_mapped_nopcr_mapq1_cis_1kb = uniq_mapped_nopcr_mapq1_cis_1kb + 1
                                            if abs(read.template_length) >= 20000:
                                                uniq_mapped_nopcr_mapq1_cis_20kb = uniq_mapped_nopcr_mapq1_cis_20kb + 1
                                    if read.mapping_quality >= 30:
                                        uniq_mapped_nopcr_mapq30 = uniq_mapped_nopcr_mapq30 + 1
                                        if read.reference_name != read.next_reference_name:
                                            uniq_mapped_nopcr_mapq30_trans = uniq_mapped_nopcr_mapq30_trans + 1
                                        else:
                                            uniq_mapped_nopcr_mapq30_cis = uniq_mapped_nopcr_mapq30_cis + 1
                                            if abs(read.template_length) >= 1000:
                                                uniq_mapped_nopcr_mapq30_cis_1kb = uniq_mapped_nopcr_mapq30_cis_1kb + 1
                                                if abs(read.template_length) >= 20000:
                                                    uniq_mapped_nopcr_mapq30_cis_20kb = uniq_mapped_nopcr_mapq30_cis_20kb + 1

    with open(out_summary, 'w') as out:
        line = "TotalReads:\t" + str(total_reads) + "\t100%" \
               + "\nTotalMapped:\t" + str(total_mapped) + "\t" + str(total_mapped * 100.0 / total_reads) + "%" \
               + "\nUniqMapped:\t" + str(uniq_mapped) + "\t" + str(uniq_mapped * 100.0 / total_reads) + "%" \
               + "\nUniqMappedToLambda:\t" + str(uniq_reads_to_lambda) + "\t" + str(
            uniq_reads_to_lambda * 100.0 / total_reads) + "%" \
               + "\nUniqPairedToLambda:\t" + str(uniq_pairs_to_lambda) + "\t" + str(
            uniq_pairs_to_lambda * 100.0 / total_reads) + "%" \
               + "\nUniqMappedToTargetSpecies:\t" + str(uniq_mapped_to_target_species) + "\t" + str(
            uniq_mapped_to_target_species * 100.0 / total_reads) + "%" \
               + "\nUniqPairedToTargetSpecies:\t" + str(uniq_paired_to_target_species) + "\t" + str(
            uniq_paired_to_target_species * 100.0 / total_reads) + "%" \
               + "\nUniqMappedNoPcr:\t" + str(uniq_mapped_nopcr) + "\t" + str(
            uniq_mapped_nopcr * 100.0 / uniq_paired_to_target_species) + "%" \
               + "\nUniqMappedNoPcrMinMapQ1:\t" + str(uniq_mapped_nopcr_mapq1) + "\t" + str(
            uniq_mapped_nopcr_mapq1 * 100.0 / uniq_paired_to_target_species) + "%" \
               + "\nUniqMappedNoPcrMinMapQ1Trans:\t" + str(uniq_mapped_nopcr_mapq1_trans) + "\t" + str(
            uniq_mapped_nopcr_mapq1_trans * 100.0 / uniq_mapped_nopcr_mapq1) + "%" \
               + "\nUniqMappedNoPcrMinMapQ1Cis:\t" + str(uniq_mapped_nopcr_mapq1_cis) + "\t" + str(
            uniq_mapped_nopcr_mapq1_cis * 100.0 / uniq_mapped_nopcr_mapq1) + "%" \
               + "\nUniqMappedNoPcrMinMapQ1CisMore1kb:\t" + str(uniq_mapped_nopcr_mapq1_cis_1kb) + "\t" + str(
            uniq_mapped_nopcr_mapq1_cis_1kb * 100.0 / uniq_mapped_nopcr_mapq1_cis) + "%" \
               + "\nUniqMappedNoPcrMinMapQ1CisMore20kb:\t" + str(uniq_mapped_nopcr_mapq1_cis_20kb) + "\t" + str(
            uniq_mapped_nopcr_mapq1_cis_20kb * 100.0 / uniq_mapped_nopcr_mapq1_cis) + "%" \
               + "\nUniqMappedNoPcrMinMapQ30:\t" + str(uniq_mapped_nopcr_mapq30) + "\t" + str(
            uniq_mapped_nopcr_mapq30 * 100.0 / uniq_paired_to_target_species) \
               + "\nUniqMappedNoPcrMinMapQ30Trans:\t" + str(uniq_mapped_nopcr_mapq30_trans) + "\t" + str(
            uniq_mapped_nopcr_mapq30_trans * 100.0 / uniq_mapped_nopcr_mapq30) + "%" \
               + "\nUniqMappedNoPcrMinMapQ30Cis:\t" + str(uniq_mapped_nopcr_mapq30_cis) + "\t" + str(
            uniq_mapped_nopcr_mapq30_cis * 100.0 / uniq_mapped_nopcr_mapq30) + "%" \
               + "\nUniqMappedNoPcrMinMapQ30CisMore1kb:\t" + str(uniq_mapped_nopcr_mapq30_cis_1kb) + "\t" + str(
            uniq_mapped_nopcr_mapq30_cis_1kb * 100.0 / uniq_mapped_nopcr_mapq30_cis) + "%" \
               + "\nUniqMappedNoPcrMinMapQ30CisMore20kb:\t" + str(uniq_mapped_nopcr_mapq30_cis_20kb) + "\t" + str(
            uniq_mapped_nopcr_mapq30_cis_20kb * 100.0 / uniq_mapped_nopcr_mapq30_cis) + "%\n"
        out.write(line)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--in_cram')
    p.add_argument('--out_summary')
    args = p.parse_args()
    mh_reads_summary(**vars(args))
