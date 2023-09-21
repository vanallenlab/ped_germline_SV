#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023 Ryan L. Collins, Riaz Gillani, Jett Crowdis and the Van Allen Laboratory
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
from scipy.stats import chisquare


ascii_idx = {l : i + 1 for i, l in enumerate(string.ascii_lowercase)}


def melt_meta(meta, keepers):
    """
    Melt metadata pd.DataFrame into dict format expected by balance_ancestries()
    """

    def __build_entry(sdat):
        if sdat.disease == 'control' \
        and (sdat.proband == 'No' or sdat.isnull().proband):
            pheno = 'control'
        elif sdat.disease != 'control' \
        and (sdat.proband == 'Yes' or sdat.isnull().proband):
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


def remove_sample(samples, sid):
    """
    Remove a sample from samples dict while accounting for relatives / "complete_trio" label
    """

    samples.pop(sid)

    for oid, values in samples.items():
        if 'trio_members' not in values.keys():
            continue
        if sid in values['trio_members']:
            samples[rel]['complete_trio'] = False

    return samples


def balance_ancestries(samples, require={}, seed=2023, abs_fold=True,
                       allow_drop_cases=True, return_keepers=False, verbose=False):
    """
    Downsample controls to balance ancestry distributions between cases and controls
    """

    # Subset samples based on requirements, if optioned
    samples_subset = samples.copy()
    for sid, vals in samples.items():
        for key, value in require.items():
            if isinstance(value, list):
                keep = vals[key] in value
            else:
                keep = vals[key] == value
            if not keep:
                samples_subset.pop(sid)
                break

    def __count_pops(samples):
        pops = {}
        for vals in samples.values():
            pop = vals['pop']
            if pop not in pops.keys():
                pops[pop] = {True : 0, False : 0}
            is_case = (vals['pheno'] != 'control')
            pops[pop][is_case] += 1
        return pd.DataFrame.from_dict(pops)

    def __compare_pops(counts):
        k_case = counts.loc[True, :].values
        k_ctrl = counts.loc[False, :].values
        ratio = sum(k_case) / sum(k_ctrl)
        x2 = chisquare(k_case, ratio * k_ctrl)
        return x2.pvalue

    def __get_fold(counts, abs_fold=True):
        fold = counts.loc[False, :] / counts.loc[True, :]
        if abs_fold:
            order = np.log2(fold).abs().sort_values(ascending=False).index
        else:
            order = np.log2(fold).sort_values(ascending=False).index
        return fold[order]

    counts = __count_pops(samples_subset)
    if verbose:
        print('Samples prior to pruning:')
        print(counts)

    # Drop all samples from a population that has strictly zero cases or controls
    sids_to_prune = set()
    pop_is_quasi_empty = (counts == 0).any()
    for pop in pop_is_quasi_empty.index[pop_is_quasi_empty].tolist():
        for sid, vals in samples_subset.items():
            if vals['pop'] == pop:
                sids_to_prune.add(sid)
    for sid in sids_to_prune:
        sample_subset = remove_sample(samples_subset, sid)

    counts = counts.loc[:, ~pop_is_quasi_empty]
    pval = __compare_pops(counts)
    random.seed(seed)
    while pval < 0.05 or np.isnan(pval):
        fold = __get_fold(counts, abs_fold=abs_fold)
        pop = fold.index[0]
        if fold[pop] >= 1 or not allow_drop_cases:
            pop_sids = [sid for sid, vals in samples_subset.items() \
                        if vals['pop'] == pop and vals['pheno'] == 'control']
        else:
            pop_sids = [sid for sid, vals in samples_subset.items() \
                        if vals['pop'] == pop and vals['pheno'] != 'control']
        loser = random.choice(pop_sids)
        sids_to_prune.add(loser)
        sample_subset = remove_sample(samples_subset, loser)
        counts = __count_pops(samples_subset)
        pval = __compare_pops(counts)

    if verbose:
        print('Samples after pruning:')
        print(counts)

    if return_keepers:
        return samples_subset

    for sid in sids_to_prune:
        samples = remove_sample(samples, sid)

    return samples


def rebalance_cancer(samples, cancer):
    """
    Conduct ancestry matching separately for case:control and trio arms 
    """

    seed = int(''.join([str(ascii_idx[l]) for l in cancer]))

    n_ctrl_start = len([s for s, v in samples.items() if v['pheno'] == 'control'])

    print('\n\nNow pruning ' + cancer + ':')
    keepers = []

    def _pheno_match(svals, pheno):
        if pheno == 'pancan':
            return svals['pheno'] in 'ewing neuroblastoma osteosarcoma'.split()
        else:
            return svals['pheno'] == pheno

    if cancer == 'pancan':
        req_phenos = 'ewing neuroblastoma osteosarcoma control'.split()
    else:
        req_phenos = [cancer, 'control']

    if len([v for v in samples.values() if _pheno_match(v, cancer) and v['validation']]) > 0:
        print('\nDiscovery cohort summary:')
        validation_samples = \
            balance_ancestries(samples.copy(), 
                               require={'pheno' : req_phenos, 
                                        'validation' : True}, 
                               seed=seed, abs_fold=False, 
                               allow_drop_cases=False, 
                               return_keepers=True, verbose=True)
        keepers += [s for s, v in validation_samples.items() if v ['pheno'] == 'control']

    if len([v for v in samples.values() if _pheno_match(v, cancer) and not v['validation']]) > 0:
        print('\nTrio/replication cohort summary:')
        trio_samples = \
            balance_ancestries(samples.copy(), 
                               require={'pheno' : req_phenos, 
                                        'validation' : False}, 
                               seed=seed, abs_fold=False, 
                               allow_drop_cases=False, 
                               return_keepers=True, verbose=True)
        keepers += [s for s, v in trio_samples.items() if v ['pheno'] == 'control']
    
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
    cancers = ['pancan'] + [x for x in meta.disease.unique() if x != 'control']
    for cancer in cancers:
        control_ids[cancer] = rebalance_cancer(samples, cancer)

    # Add columns to meta indicating control status for each disease
    for cancer, c_ids in control_ids.items():
        meta[cancer + '_control'] = meta.sample_id.isin(c_ids)
    
    # Write to outfile
    meta.rename(columns={'sample_id' : 'entity:sample_id'}).\
         to_csv(args.outfile, sep='\t', index=False, na_rep='NA')

if __name__ == '__main__':
    main()
