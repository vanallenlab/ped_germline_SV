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
import numpy as np
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


def calc_ncr(record):
    """
    Calculate no-call rate for a record's genotypes
    """

    total = 0
    nocalls = 0
    for sinfo in record.samples.values():
        total += 1

        if all([a is None for a in sinfo.get('GT', (None, ))]):
            nocalls += 1

    if total > 0:
        return nocalls / total
    else:
        return None


def recalibrate_qual(record):
    """
    Recalibrate record QUAL (quality)
    Defined as median GQ among all non-reference GTs
    """

    gqs = []
    for sinfo in record.samples.values():
        a = [a for a in sinfo.get('GT', (None, None)) if a is not None]
        if np.nansum(a) > 0:
            gqs.append(sinfo.get('GQ', None))

    if len(gqs) > 0:
        return np.nanmedian(gqs)
    else:
        return 0


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
    # header.filters.remove_header('HIGH_PCRMINUS_NOCALL_RATE')
    header.add_line('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">')
    header.add_line('##INFO=<ID=OLD_ID,Number=1,Type=String,Description="Original GATK-SV variant ID before polishing">')

    # Open connection to output vcf
    if args.vcf_out in '- stdout /dev/stdout':
        outvcf = pysam.VariantFile(stdout, 'w', header=header)
    else:
        outvcf = pysam.VariantFile(args.vcf_out, 'w', header=header)

    # Iterate over records in invcf, clean up, and write to outvcf
    svtype_counter = {}
    records_seen = {}
    for record in invcf.fetch():

        # Skip duplicated records (defined as those with identical chrom/pos/end/svtype)
        is_dup, records_seen = check_duplicate(record, records_seen)
        if is_dup:
            continue

        # Get reused record info
        svtype = record.info.get('SVTYPE')
        svlen = record.info.get('SVLEN', 0)

        # Skip empty records
        if record.info.get('AC', (1, ))[0] == 0 \
        and (svtype != 'CNV' or 'MULTIALLELIC' not in record.filter.keys()):
            continue

        # Label small MCNVs as UNRESOLVED
        if (svtype == 'CNV' or 'MULTIALLELIC' in record.filter.keys()) \
        and svlen < 5000:
            record.filter.add('UNRESOLVED')

        # Label very large, common, unbalanced SVs as UNRESOLVED
        if svtype in 'DEL DUP INS CPX'.split() \
        and svlen >= 500000 \
        and record.info.get('AF', (0, ))[0] >= 0.05:
            record.filter.add('UNRESOLVED')

        # Label common CTX as UNRESOLVED
        if svtype == 'CTX' \
        and record.info.get('AF', (0, ))[0] >= 0.01:
            record.filter.add('UNRESOLVED')

        # Apply stricter NCR filter for artifact deletion peak
        try:
            if svtype == "DEL" \
            and svlen > 400 \
            and svlen < 1000 \
            and record.info.get('AC', (0, ))[0] / record.info.get('AN', 1) < 0.1:
                if record.info.get('PCRMINUS_NCR', 0) > 0.01:
                    record.filter.add('HIGH_NCR')
        except:
            import pdb; pdb.set_trace()

        # Recompute NCR
        if svtype != 'CNV' and 'MULTIALLELIC' not in record.filter.keys():
            record.info['NCR'] = calc_ncr(record)

        # Clear old high NCR FILTER and reannotate based on updated NCR
        if 'HIGH_NCR' in record.filter.keys():
            original_filters = [k for k in record.filter.keys()]
            record.filter.clear()
            for k in original_filters:
                if k in 'HIGH_NCR HIGH_PCRMINUS_NOCALL_RATE'.split():
                    record.filter.add(k)
        if record.info.get('NCR', 0) >= 0.1:
            record.filter.add('HIGH_NCR')

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

        # Recalibrate QUAL score
        if svtype != 'CNV' and 'MULTIALLELIC' not in record.filter.keys():
            record.qual = recalibrate_qual(record)

        # Rename record
        if record.chrom not in svtype_counter.keys():
            svtype_counter[record.chrom] = {}
        if svtype not in svtype_counter[record.chrom].keys():
            svtype_counter[record.chrom][svtype] = 0
        svtype_counter[record.chrom][svtype] += 1
        new_id = '_'.join(['PedSV.v2.0', svtype, record.chrom,
                           str(svtype_counter[record.chrom][svtype])])
        record.info['OLD_ID'] = record.id
        record.id = new_id

        # Write to outvcf
        outvcf.write(record)

    # Close connection to output file to clear buffer
    outvcf.close()


if __name__ == '__main__':
    main()

