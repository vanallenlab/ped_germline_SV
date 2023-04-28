#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023 Ryan L. Collins, Riaz Gillani and the Van Allen Laboratory
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Downsample control samples to ancestry-match individual pediatric cancer subtypes
"""


import argparse
import numpy as np
import pandas as pd
import random
import string
from optimize_sample_inclusion_lists import balance_ancestries
from scipy.stats import chisquare


ascii_idx = {l : i + 1 for i, l in enumerate(string.ascii_lowercase)}


def melt_meta(meta, keepers):
    """
    Melt metadata pd.DataFrame into dict format expected by balance_ancestries()
    """

    def __build_entry(sdat):
        if sdat.study_phase == 'case_control':
            pheno = sdat.disease
        else:
            if sdat.disease == 'control' \
            and sdat.proband == 'No':
                pheno = 'control'
            elif sdat.disease != 'control' \
            and sdat.proband == 'Yes':
                pheno = sdat.disease
            else:
                pheno = None

        return {'pop' : sdat.ancestry_inferred_by_SVs,
                'pheno' : pheno,
                'validation' : sdat.study_phase == 'case_control'}


    sample_dict = {}
    for idx, sdat in meta.loc[meta.sample_id.isin(keepers), :].iterrows():
         sample_dict[sdat.sample_id] = __build_entry(sdat)

    return sample_dict


def rebalance_cancer(samples, cancer):
    """
    Conduct ancestry matching separately for case:control and trio arms 
    """

    seed = int(''.join([str(ascii_idx[l]) for l in cancer]))

    n_ctrl_start = len([s for s, v in samples.items() if v['pheno'] == 'control'])

    print('\n\nNow pruning ' + cancer + ':')
    keepers = []

    if len([v for v in samples.values() if v['pheno'] == cancer and not v['validation']]) > 0:
        print('\nTrio cohort summary:')
        trio_samples = \
            balance_ancestries(samples.copy(), 
                               require={'pheno' : [cancer, 'control'], 
                                        'validation' : False}, 
                               seed=seed, abs_fold=False, 
                               allow_drop_cases=False, 
                               return_keepers=True, verbose=True)
        keepers += [s for s, v in trio_samples.items() if v ['pheno'] == 'control']

    if len([v for v in samples.values() if v['pheno'] == cancer and v['validation']]) > 0:
        print('\nValidation cohort summary:')
        validation_samples = \
            balance_ancestries(samples.copy(), 
                               require={'pheno' : [cancer, 'control'], 
                                        'validation' : True}, 
                               seed=seed, abs_fold=False, 
                               allow_drop_cases=False, 
                               return_keepers=True, verbose=True)
        keepers += [s for s, v in validation_samples.items() if v ['pheno'] == 'control']
    
    n_ctrl_end = len(keepers)

    print('\nRetained {:,} of {:,} controls for {}'.format(n_ctrl_end, n_ctrl_start, cancer))

    return keepers


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--metadata', required=True, help='sample metadata .tsv')
    parser.add_argument('--keep-samples', help='list of IDs to retain during ' +
                        'ancestry matching default: keep all samples')
    parser.add_argument('--outfile', required=True, help='output .tsv')
    args = parser.parse_args()

    # Load sample metadata
    meta = pd.read_csv(args.metadata, sep='\t').\
              rename(columns={'entity:sample_id' : 'sample_id'})

    # Load list of samples to keep
    if args.keep_samples is None:
        keepers = set(meta.sample_id.values.tolist())
    else:
        with open(args.keep_samples) as keep_fin:
            keepers = set([x.rstrip() for x in keep_fin.readlines()])

    # Melt sample metadata to dict for compliance with balance_ancestries()
    samples = melt_meta(meta, keepers)

    # Rebalance each disease and cohort separately
    control_ids = {}
    for cancer in [x for x in meta.disease.unique() if x != 'control']:
        control_ids[cancer] = rebalance_cancer(samples, cancer)

    # Add columns to meta indicating control status for each disease
    for cancer, c_ids in control_ids.items():
        meta[cancer + '_control'] = meta.sample_id.isin(c_ids)
    
    # Write to outfile
    meta.rename(columns={'sample_id' : 'entity:sample_id'}).\
         to_csv(args.outfile, sep='\t', index=False)

if __name__ == '__main__':
    main()
