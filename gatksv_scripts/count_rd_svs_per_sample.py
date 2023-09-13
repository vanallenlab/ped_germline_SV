#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023 Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Count number of read depth-based CNVs per sample from an svtk vcf2bed output
"""


import argparse
import gzip
import pandas as pd
from sys import stdout


cnv_types = 'DEL DUP'.split()


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('bed', help='input .bed [can be bgzipped]')
    parser.add_argument('-o', '--outfile', help='output .tsv [default: stdout]',
                        default='stdout')
    parser.add_argument('-s', '--samples-list', help='list of samples to evaluate')
    args = parser.parse_args()

    # Build collector for results
    counts = {cnv : {} for cnv in cnv_types}
    if args.samples_list is not None:
        with open(args.samples_list) as slin:
            for line in slin.readlines():
                sample = line.rstrip()
                for cnv in cnv_types:
                    if sample not in counts[cnv].keys():
                        counts[cnv][sample] = 0

    # Process BED in one shot
    with gzip.open(args.bed, 'rt') as bedin:
        # Get index of header columns
        col_idxs = {col : i for i, col in enumerate(bedin.readline().rstrip().split('\t'))}

        # Iterate over BED records and only count those meeting RD criteria
        for rec in bedin.readlines():
            rvals = {k : v for k, v in zip(col_idxs.keys(), rec.rstrip().split('\t'))}
            svtype = rvals.get('SVTYPE')
            if (rvals.get('ALGORITHMS') == 'depth' \
                or rvals.get('EVIDENCE') in 'RD BAF,RD'.split()) \
            and svtype in cnv_types:
                for sid in rvals.get('samples').split(','):
                    if sid not in counts[svtype].keys():
                        counts[svtype][sid] = 0
                    counts[svtype][sid] += 1

    # Write results to output file
    if args.outfile in '- /dev/stdout stdout':
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')
    outfile.write('#sample\tcnv_type\tcount\n')
    for cnv in counts.keys():
        for sid, k in counts[cnv].items():
            outfile.write('{}\t{}\t{}\n'.format(sid, cnv, k))
    outfile.close()


if __name__ == '__main__':
    main()
