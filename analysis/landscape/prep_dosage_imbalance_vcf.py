#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2024 Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Convert a VCF into a simple sites-only VCF with deletion and duplication imbalance annotations
"""


import argparse
import pysam
from sys import stdout


def parse_cpx_intervals(record):
    """
    Compute the sum of deleted and duplicated nucleotides for a complex SV
    """

    idel, idup = 0, 0

    for interval in record.info.get('CPX_INTERVALS', ''):
        itype = interval.split('_')[0]
        if itype in 'DEL DUP'.split():
            coords = interval.split(':')[1]
            pos = sorted([int(k) for k in coords.split('-')], reverse=True)
            size = abs(pos[0] - pos[1])
            if itype == 'DEL':
                idel += size
            elif itype == 'DUP':
                idup += size

    if idel > 0:
        record.info['DEL_IMBALANCE'] = idel
    if idup > 0:
        record.info['DUP_IMBALANCE'] = idup

    return record


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcfin', help='input .vcf')
    parser.add_argument('--vcf-out', help='output .vcf [default: stdout]', 
                        default='stdout')
    args = parser.parse_args()

    # Open connection to input VCF
    invcf = pysam.VariantFile(args.vcfin)
    
    # Prepare header for output VCF
    header = invcf.header
    header.add_line('##INFO=<ID=DEL_IMBALANCE,Number=1,Type=Integer,Description="Nucleotides deleted.">')
    header.add_line('##INFO=<ID=DUP_IMBALANCE,Number=1,Type=Integer,Description="Nucleotides duplicated.">')

    # Open connection to output VCF
    if args.vcf_out in '- /dev/stdout stdout':
        outfile = stdout
    else:
        outfile = args.vcf_out
    outvcf = pysam.VariantFile(outfile, 'w', header=header)

    # Process each record in input VCF
    for record in invcf.fetch():

        # Compute genomic imbalance annotations per record
        if record.info.get('SVTYPE', '') == 'DEL':
            record.info['DEL_IMBALANCE'] = record.info.get('SVLEN', 0)
        elif record.info.get('SVTYPE', '') == 'DUP':
            record.info['DUP_IMBALANCE'] = record.info.get('SVLEN', 0)
        elif record.info.get('SVTYPE', '') == 'CPX':
            record = parse_cpx_intervals(record)
        else:
            continue

        # Strip all other info
        for key in record.info.keys():
            if key in 'CHROM2 POS2 END SVTYPE SVLEN DEL_IMBALANCE DUP_IMBALANCE AC AF'.split():
                continue
            else:
                record.info.pop(key)

        # Write to output VCF
        outvcf.write(record)

    # Close connection to output VCF
    outvcf.close()


if __name__ == '__main__':
    main()
