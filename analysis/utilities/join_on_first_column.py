#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2024 Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Join two or more .tsvs on their first column
"""


import argparse
import pandas as pd
from sys import stdout


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('tsvs', nargs='+', help='.tsv files to join')
    parser.add_argument('--join', help='how should files be joined? [default: inner]', 
                        choices='left right inner outer cross'.split(),
                        default='inner')
    parser.add_argument('--tsv-out', help='output .tsv [default: stdout]',
                        default='stdout')
    args = parser.parse_args()

    # Seed data.frame is first tsv
    m = pd.read_csv(args.tsvs[0], sep='\t')

    # Sequentially join each other .tsv passed as a positional argument
    if len(args.tsvs[0]) > 1:
        for fin in args.tsvs[1:]:
            other = pd.read_csv(fin, sep='\t')
            m = m.merge(other, how=args.join, left_on=m.columns[0], 
                        right_on=other.columns[0])

    # Write results to output file
    if args.tsv_out in '- /dev/stdout stdout':
        outfile = stdout
    else:
        outfile = args.tsv_out
    m.to_csv(outfile, sep='\t', na_rep='NA', index=False)


if __name__ == '__main__':
    main()
