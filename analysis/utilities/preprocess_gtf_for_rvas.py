#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023 Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Preprocess a GTF for single-chromosome RVAS
"""


import argparse
import pybedtools as pbt
from sys import stdout


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('gtf', help='input .bed [can be bgzipped]')
    parser.add_argument('-c', '--chromosome', help='chromosome to extract ' +
                        '[default: keep all genes]')
    parser.add_argument('-x', '--exclusion-bed', help='.bed of regions to exclude')
    parser.add_argument('-f', '--exclusion-frac', help='max fraction of gene covered ' +
                        'by --exclusion-bed before being excluded [default: 0.5]',
                        type=float, default=0.5)
    parser.add_argument('-o', '--outfile', help='output .tsv [default: stdout]',
                        default='stdout')
    args = parser.parse_args()

    # Open connection to gtf
    gtf = pbt.BedTool(args.gtf)

    # Subset to single chromosome, if optioned
    if args.chromosome is not None:
        gtf = gtf.filter(lambda x: x.chrom == args.chromosome)

    # Transform gtf into simple four-column bed
    gstrs = ''
    for feat in gtf.filter(lambda x: x[2] == 'gene'):
        gstrs += '{}\t{}\t{}\t{}\n'.format(feat.chrom, feat.start, feat.end,
                                           feat.attrs['gene_name'])
    genes = pbt.BedTool(gstrs, from_string=True).sort()

    # Filter genes based on exclusion bed, if optioned
    if args.exclusion_bed is not None:
        genes = genes.coverage(args.exclusion_bed).\
                      filter(lambda x: float(x[-1]) <= args.exclusion_frac)

    # Write results to output file
    if args.outfile in '- /dev/stdout stdout':
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')
    for feat in genes:
        outfile.write('\t'.join(feat.fields[:4]) + '\n')
    outfile.close()


if __name__ == '__main__':
    main()
