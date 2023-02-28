#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023 Ryan L. Collins, Riaz Gillani and the Van Allen Laboratory
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Final cleanup of pediatric WGS SV callsets
"""


filters_to_info = [('BOTHSIDES_SUPPORT', 
                    'Variant has read-level support for both sides of breakpoint'),
                   ('HIGH_SR_BACKGROUND', 
                    'High number of SR splits in background samples indicating messy region'),
                   ('PESR_GT_OVERDISPERSION', 
                    'High PESR dispersion')]


import argparse
import pysam
from sys import stdin, stdout


def check_duplicate(record, records_seen):
    """
    Checks if a record's coordinates and SVTYPE have already been observed
    Returns boolean (and updated records_seen if False)
    """

    svtype = record.info.get('SVTYPE')
    chrom = record.chrom
    start = record.pos
    end = record.stop

    is_dup = False
    if svtype not in records_seen.keys():
        records_seen[svtype] = {chrom : {start : set([end])}}
    elif chrom not in records_seen[svtype].keys():
        records_seen[svtype][chrom] = {start : set([end])}
    elif start not in records_seen[svtype][chrom]:
        records_seen[svtype][chrom][start] = set([end])
    elif end not in records_seen[svtype][chrom][start]:
        records_seen[svtype][chrom][start].add(end)
    else:
        is_dup = True

    return is_dup, records_seen


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf_in', help='input .vcf')
    parser.add_argument('vcf_out', help='output .vcf')
    args = parser.parse_args()

    # Open connection to input vcf
    if args.vcf_in in '- stdin /dev/stdin'.split():
        invcf = pysam.VariantFile(stdin)
    else:
        invcf = pysam.VariantFile(args.vcf_in)
    samples = [s for s in invcf.header.samples]

    mf_infos = [k for k in invcf.header.info.keys() \
                if k.startswith('MALE_') or k.startswith('FEMALE_')]

    # Reformat header
    header = invcf.header
    for key, descrip in filters_to_info:
        if key in header.filters.keys():
            header.filters.remove_header(key)
            header.add_line('##INFO=<ID={},Number=0,Type=Flag,Description="{}">'.format(key, descrip))
    for key in mf_infos:
        header.info.remove_header(key)
    header.add_line('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">')

    # Open connection to output vcf
    if args.vcf_out in '- stdout /dev/stdout':
        outvcf = pysam.VariantFile(stdout, 'w', header=header)
    else:
        outvcf = pysam.VariantFile(args.vcf_out, 'w', header=header)

    # Iterate over records in invcf, clean up, and write to outvcf
    records_seen = {}
    for record in invcf.fetch():

        # Skip duplicated records (defined as those with identical chrom/pos/end/svtype)
        is_dup, records_seen = check_duplicate(record, records_seen)
        if is_dup:
            continue

        # Label small MCNVs as UNRESOLVED
        if record.info.get('SVTYPE') == 'CNV' \
        and record.info.get('SVLEN', 0) < 5000:
            record.filter.add('UNRESOLVED')

        # Label large, common, complex SVs as UNRESOLVED
        if record.info.get('SVTYPE') == 'CPX' \
        and record.info.get('SVLEN', 0) >= 50000 \
        and record.info.get('AF', 0) >= 0.1:
            record.filter.add('UNRESOLVED')

        # Apply stricter NCR filter for artifact deletion peak
        if record.info.get('SVTYPE') == "DEL" \
        and record.info.get('SVLEN', 0) > 400 \
        and record.info.get('SVLEN', 0) < 1000 \
        and record.info.get('AC', (0, ))[0] / record.info.get('AN', 1) < 0.1:
            if record.info.get('PCRMINUS_NCR', 0) > 0.01:
                record.filter.add('HIGH_PCRMINUS_NOCALL_RATE')

        # Clear all MALE/FEMALE AF annotations
        for key in mf_infos:
            if key in record.info.keys():
                record.info.pop(key)

        # Relocate certain FILTERs to INFOs
        for key, descrip in filters_to_info:
            if key in record.filter.keys():
                record.info[key] = True
        old_filts = set(record.filter.keys()).difference(set([k for k, v in filters_to_info]))
        if len(old_filts) == 0:
            old_filts = set(['PASS'])
        record.filter.clear()
        for filt in old_filts:
            record.filter.add(filt)

        # Write to outvcf
        outvcf.write(record)

    # Close connection to output file to clear buffer
    outvcf.close()


if __name__ == '__main__':
    main()

