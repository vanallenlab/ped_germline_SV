#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023 Ryan L. Collins and the Van Allen Laboratory
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Subset and rename records in a VCF
"""


import argparse
import csv
import pysam


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--vcf-in', required=True, help='input .vcf')
    parser.add_argument('--vid-map', required=True, help='two-columns .tsv ' + 
                        'mapping old to new variant IDs')
    parser.add_argument('--vcf-out', required=True, help='output .vcf')
    args = parser.parse_args()

    # Open connections to input and output VCFs
    invcf = pysam.VariantFile(args.vcf_in)
    outvcf = pysam.VariantFile(args.vcf_out, mode='w', header=invcf.header)

    # Load variant ID map
    vid_map = {}
    with open(args.vid_map) as tsvin:
        for old, new in csv.reader(tsvin, delimiter='\t'):
            vid_map[old] = new

    # Iterate over input VCF and subset/rename records
    for record in invcf.fetch():
        if record.id not in vid_map.keys():
            continue
        record.id = vid_map[record.id]
        outvcf.write(record)

    # Clear buffer
    outvcf.close()


if __name__ == '__main__':
    main()
