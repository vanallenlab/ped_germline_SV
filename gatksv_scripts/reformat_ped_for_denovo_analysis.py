#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2024 Ryan L. Collins and the Van Allen Laboratory
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Reformat a .ped file for the Talkowski Lab de novo SV pipeline
"""


import argparse
import csv


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--ped-in', required=True, help='input .ped')
    parser.add_argument('--keep-samples', required=True, help='.txt list of ' +
                        'samples to retain.')
    parser.add_argument('--ped-out', required=True, help='output .ped')
    args = parser.parse_args()

    # Read list of eligible samples
    with open(args.keep_samples) as fin:
        keepers = [sid.rstrip() for sid in fin.readlines()]

    # Open connections to input and output .ped files
    inped = open(args.ped_in)
    outped = open(args.ped_out, 'w')
    header = 'FamID IndividualID FatherID MotherID Gender Affected'.split()
    outped.write('\t'.join(header) + '\n')

    # Iterate over input .ped file and only retain lines that are complete trios
    for famid, proid, faid, moid, sex, aff in csv.reader(inped, delimiter='\t'):
        if proid not in keepers \
        or faid not in keepers \
        or moid not in keepers:
            continue
        outped.write('\t'.join([famid, proid, faid, moid, sex, aff]) + '\n')
        outped.write('\t'.join([famid, faid, '0', '0', '1', '1']) + '\n')
        outped.write('\t'.join([famid, moid, '0', '0', '2', '1']) + '\n')

    # Clear buffer
    outped.close()


if __name__ == '__main__':
    main()
