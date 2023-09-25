#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023 Ryan L. Collins, Riaz Gillani, Jett Crowdis and the Van Allen Laboratory
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Parse hg38 alternate haplotype alignments vs. primary assembly

Expects input schema as defined by UCSC:
http://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=map&hgta_track=altSeqLiftOverPsl&hgta_table=altSeqLiftOverPsl&hgta_doSchema=describe+table+schema
"""


import argparse
import pandas as pd
import pybedtools as pbt
from sys import stdout


columns = 'bin matches misMatches repMatches nCount qNumInsert qBaseInsert ' + \
          'tNumInsert tBaseInsert strand qName qSize qStart qEnd tName tSize ' + \
          'tStart tEnd blockCount blockSizes qStarts tStarts'
columns = columns.split()
primary_contigs = ['chr' + str(i+1) for i in range(22)] + ['chrX', 'chrY']


def parse_alignment(aln_info):
    """
    Parse a single row in an alignment table

    Return a pbt.BedTool of regions
    """

    chrom = aln_info.tName
    starts = [int(x) for x in aln_info.tStarts.split(',') if x != '']
    sizes = [int(x) for x in aln_info.blockSizes.split(',') if x != '']
    ends = [x + k for x, k in zip(starts, sizes)]
    intervals = ['\t'.join([chrom, str(a), str(b)]) for a, b in zip(starts, ends)]

    return pbt.BedTool('\n'.join(intervals), from_string=True)


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('tsv', help='input .tsv from UCSC')
    parser.add_argument('-o', '--outfile', help='output BED [can be gzipped]', 
                        default='stdout')
    parser.add_argument('-g', '--genome', help='genome file [for sorting]')
    args = parser.parse_args()

    # Load tsv as pd.DataFrame
    df = pd.read_csv(args.tsv, names=columns, sep='\t',
                     usecols='tName strand tStarts blockSizes'.split())

    # Restrict to primary contigs
    df = df[df.tName.isin(primary_contigs)]

    # Parse each row
    bts = df.apply(parse_alignment, axis=1)
    bed = bts[0].cat(*bts[1:], postmerge=True)
    if args.genome is not None:
        bed = bed.sort(g=args.genome)

    # Write to output file
    if args.outfile in '- stdout /dev/stdout'.split():
        for line in bed:
            stdout.write('\t'.join(line.fields) + '\n')
    else:
        bed.saveas(args.outfile)


if __name__ == '__main__':
    main()

