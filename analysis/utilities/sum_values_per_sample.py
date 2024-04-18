#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2024 Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Sum a list of values per sample
"""


import argparse
from sys import stdin, stdout


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--tsv-in', help='input .tsv [default: stdin]', 
                        default='stdin')
    parser.add_argument('--tsv-out', help='output .tsv [default: stdout]',
                        default='stdout')
    args = parser.parse_args()

    # Prep collection dictionary
    res = {}

    # Open connection to input file
    if args.tsv_in in '- /dev/stdin stdin':
        infile = stdin
    else:
        infile = open(args.tsv_in)

    # Read each line from input file
    for line in infile.readlines():
        sample, value = line.rstrip().split('\t')
        if sample not in res.keys():
            res[sample] = 0
        res[sample] += float(value)

    # Write results to output file
    if args.tsv_out in '- /dev/stdout stdout':
        outfile = stdout
    else:
        outfile = open(args.tsv_out, 'w')
    for sample, value in sorted(res.items()):
        outfile.write('\t'.join([sample, str(value)]) + '\n')
    outfile.close()


if __name__ == '__main__':
    main()
