#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023 Ryan L. Collins, Riaz Gillani and the Van Allen Laboratory
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Optimize lists of samples to include for final pediatric SV callsets for downstream analysis
"""


import argparse
import csv
import numpy as np
import pandas as pd
import random
from scipy.stats import chisquare


default_entry = {'in_trio_vcf' : False,
                 'trio_vcf_id' : None,
                 'in_validation_vcf' : False,
                 'validation_vcf_id' : None,
                 'was_trio' : False,
                 'trio_members' : [],
                 'complete_trio' : False,
                 'member' : 'proband',
                 '1000G' : False,
                 'pop' : 'OTH',
                 'pheno' : None}


def _clean_sid(sid):
    """
    Remove bcftools artifact from sample ID
    """

    if sid.startswith('2:'):
        sid = sid[2:]

    return sid


def summarize_cohort(samples, stage='Progress', verbose=False):
    """
    Collect and report a snapshot of sample membership
    """

    n_all, n_trio, n_trio_trios, n_validation, n_validation_trios, n_both, n_both_trios = [0] * 7
    for key, vals in samples.items():
        n_all += 1
        if vals['in_trio_vcf']:
            n_trio += 1
            if vals['complete_trio'] and vals['member'] == 'proband':
                n_trio_trios += 1
        if vals['in_validation_vcf']:
            n_validation += 1
            if vals['complete_trio'] and vals['member'] == 'proband':
                n_validation_trios += 1
        if vals['in_trio_vcf'] and vals['in_validation_vcf']:
            n_both += 1
            if vals['complete_trio'] and vals['member'] == 'proband':
                n_both_trios += 1

    if verbose:
        print('\n{}. Summary:'.format(stage))
        print('  * Found a total of {:,} unique sample IDs'.format(n_all))
        print('  * {:,} samples are present in trio VCF'.format(n_trio))
        print('  * {:,} complete trios are present in trio VCF'.format(n_trio_trios))
        print('  * {:,} samples are present in validation VCF'.format(n_validation))
        print('  * {:,} complete trios are present in validation VCF'.format(n_validation_trios))
        print('  * {:,} samples are present in both VCFs'.format(n_both))
        print('  * {:,} complete trios are present in both VCFs'.format(n_both_trios))
    else:
        print('\t'.join([stage] + \
                        ['{:,}'.format(x) for x in \
                         [n_all, n_trio, n_trio_trios, n_validation, 
                          n_validation_trios, n_both, n_both_trios]]))


def load_sample_universe(trio_samples_in, validation_samples_in, ped_in):
    """
    Build a dict of all sample IDs with their callset membership
    """

    # Load all samples in trio VCF
    with open(trio_samples_in) as fin:
        samples = {_clean_sid(sid.rstrip()) : default_entry.copy() \
                   for sid in fin.readlines()}
    for sid, values in samples.items():
        samples[sid]['in_trio_vcf'] = True
        samples[sid]['trio_vcf_id'] = sid

    # Add all samples in validation VCF
    with open(validation_samples_in) as fin:
        for line in fin.readlines():
            sid = _clean_sid(line.rstrip())
            if sid not  in samples.keys():
                samples[sid] = default_entry.copy()
            samples[sid]['in_validation_vcf'] = True
            samples[sid]['validation_vcf_id'] = sid

    # Add trio information
    with open(ped_in) as fin:
        for pro, fa, mo in csv.reader(fin, delimiter='\t'):
            members = [pro, fa, mo]
            if not any([sid == '0' for sid in members]):
                for sid in members:
                    if sid in samples.keys():
                        samples[sid]['was_trio'] = True
                        samples[sid]['trio_members'] = members
                for sid in [fa, mo]:
                    if sid in samples.keys():
                        samples[sid]['member'] = 'parent'
            if all([sid in samples.keys() for sid in members]):
                for sid in members:
                    samples[sid]['complete_trio'] = True

    # Summarize all loaded samples before any pruning
    summarize_cohort(samples, stage='Before pruning')

    return samples


def _unify_values(vals1, vals2):
    """
    Take the union of values for two samples
    """

    new_entry = default_entry.copy()

    for key in vals1.keys():
        subvals = [v[key] for v in [vals1, vals2]]
        if any([isinstance(sv, bool) for sv in subvals]):
            newval = any(subvals)
        elif all([sv is None for sv in subvals]):
            newval = None
        else:
            newval = [sv for sv in subvals if sv is not None][0]

        new_entry[key] = newval
    
    return new_entry


def remove_sample(samples, sid):
    """
    Remove a sample from samples dict while accounting for relatives / "complete_trio" label
    """

    samples.pop(sid)

    rels = [oid for oid, values in samples.items() if sid in values['trio_members']]
    for rel in rels:
        samples[rel]['complete_trio'] = False

    return samples


def add_cryptic_duplicates(samples, relatives_in):
    """
    Annotate cryptic duplicate samples based on genetic relatedness
    Note: this explicitly ignores pairs of samples with identical IDs 
    since such cases are obvious duplicates
    """

    pruned = set()

    for pairstr, rel in csv.reader(open(relatives_in), delimiter='\t'):
        if rel != 'duplicates':
            continue

        sid1, sid2 = pair = [_clean_sid(x) for x in pairstr.split('|')]

        if not all(sid in samples.keys() for sid in pair):
            continue
        if sid1 == sid2:
            continue

        # Check to make sure neither sample ID has been pruned already
        # If so, can skip processing this pair
        if any([sid in pruned for sid in pair]):
            continue

        # If both IDs are present in the trio cohort,
        # keep the one from a complete trio (if either are), then
        # defer to keeping affected individuals where possible, then
        # arbitrarily choose the first of the pair
        elif all([samples[sid]['in_trio_vcf'] for sid in pair]):
            if any([samples[sid]['complete_trio'] for sid in pair]):
                keep_id = [sid for sid in pair if samples[sid]['complete_trio']][0]
            elif any([samples[sid]['pheno'] != 'control' for sid in pair]):
                keep_id = [sid for sid in pair if samples[sid]['pheno'] != 'control'][0]
            else:
                keep_id = sid1
            remove_id = [sid for sid in pair if sid != keep_id][0]
            samples = remove_sample(samples, remove_id)
            pruned.add(remove_id)

        # Otherwise, just collapse all information onto a single entry and drop the second
        else:
            new_entry = _unify_values(samples[sid1], samples[sid2])
            samples[sid1] = new_entry
            samples = remove_sample(samples, sid2)
            pruned.add(sid2)

    return samples


def prune_duplicates(samples):
    """
    Prune duplicate samples based on the following priority:
    1. If the sample is from the 1000 Genomes Project, prune from the case/control cohort;
    2. If the sample is part of a complete trio, prune from the case/control cohort;
    3. If the sample is a parent from an incomplete trio, prune from the trio cohort;
    4. If the sample is a proband from an incomplete trio, prune from the case/control cohort;
    5. If the sample has non-European ancestry, prune from the case/control cohort;
    6. Otherwise, prune from the trio cohort
    """

    for sid, vals in samples.items():

        # Skip all samples that don't appear in both cohorts
        if not all([vals[key] for key in 'in_trio_vcf in_validation_vcf'.split()]):
            continue

        if vals['1000G']:
            vals['in_validation_vcf'] = False

        elif vals['complete_trio']:
            vals['in_validation_vcf'] = False

        elif vals['member'] != 'proband' \
        and vals['was_trio'] \
        and not vals['complete_trio']:
            vals['in_trio_vcf'] = False

        elif vals['member'] == 'proband' \
        and vals['was_trio'] \
        and not vals['complete_trio']:
            vals['in_validation_vcf'] = False

        elif vals['pop'] != 'EUR':
            vals['in_validation_vcf'] = False

        else:
            vals['in_trio_vcf'] = False

        samples[sid] = vals

    return samples


def prune_related(samples, kinship_in, require={}, max_kinship=1/8):
    """
    Identify and prune related individuals based on kinship coefficient
    """

    # Build list of eligible samples
    elig_samps = set()
    for sid, vals in samples.items():
        elig = True
        for key, value in require.items():
            if vals[key] != value:
                elig = False
                break
        if elig:
            elig_samps.add(sid)

    # Load kinship data and subset to eligible samples surpassing max_kinship
    kdf = pd.read_csv(kinship_in, sep='\t').rename(columns={'#sample_1' : 'sample_1'})
    kdf['sample_1'] = kdf.sample_1.map(_clean_sid)
    kdf['sample_2'] = kdf.sample_2.map(_clean_sid)
    kdf = kdf.loc[kdf.sample_1.isin(elig_samps) & \
                  kdf.sample_2.isin(elig_samps) & \
                  (kdf.sample_1 != kdf.sample_2) & \
                  (kdf.kin >= max_kinship), :]

    # Prune one sample at a time to maximize number of unrelated samples retained
    def __count_links(kdf):
        return pd.Series(np.concatenate([kdf.sample_1.values, kdf.sample_2.values])).value_counts()
    links = __count_links(kdf)
    while len(links) > 0:
        prune_id = links.index[0]
        samples = remove_sample(samples, prune_id)
        kdf = kdf.loc[(kdf.sample_1 != prune_id) & (kdf.sample_2 != prune_id)]
        links = __count_links(kdf)

    return samples


def balance_ancestries(samples, require={}, seed=2023):
    """
    Downsample controls to balance ancestry distributions between cases and controls
    """

    # Subset samples based on requirements, if optioned
    samples_subset = samples.copy()
    for sid, vals in samples.items():
        for key, value in require.items():
            if vals[key] != value:
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

    def __get_fold(counts):
        fold = counts.loc[False, :] / counts.loc[True, :]
        order = np.log2(fold).abs().sort_values(ascending=False).index
        return fold[order]

    counts = __count_pops(samples_subset)
    pval = __compare_pops(counts)
    random.seed(seed)
    sids_to_prune = set()
    while pval < 0.05 or np.isnan(pval):
        fold = __get_fold(counts)
        pop = fold.index[0]
        if fold[pop] > 1:
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

    for sid in sids_to_prune:
        samples = remove_sample(samples, sid)

    return samples


def write_sample_list(samples, require, out_path):
    """
    Write a list of sample IDs meeting criteria to out_path
    """

    fout = open(out_path, 'w')

    for sid, vals in samples.items():
        for key, value in require.items():
            if vals[key] == value:
                fout.write(sid + '\n')

    fout.close()


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--kinship', required=True, help='relatedness metrics .tsv')
    parser.add_argument('--predicted-relatives', required=True, 
                        help='.tsv with sample pairs and predicted relationship')
    parser.add_argument('--pedfile', required=True, help='Cohort .ped file ' +
                        'truncated to only include columns 2-4 (sample IDs)')
    parser.add_argument('--samples-from-1000g', required=True, help='List of IDs ' +
                        'for samples from the 1000 Genomes Project')
    parser.add_argument('--ancestry-labels', required=True,
                        help='two-column .tsv mapping sample IDs to ancestry')
    parser.add_argument('--phenotype-labels', required=True,
                        help='two-column .tsv mapping sample IDs to phenotype')
    parser.add_argument('--trio-callset-samples', required=True,
                        help='list of IDs present in the trio callset')
    parser.add_argument('--validation-callset-samples', required=True,
                        help='list of IDs present in the trio callset')
    parser.add_argument('--out-prefix', required=True, help='output file prefix')
    args = parser.parse_args()

    # Build dict of all samples to consider and annotate with trio membership
    samples = load_sample_universe(args.trio_callset_samples,
                                   args.validation_callset_samples,
                                   args.pedfile)

    # Annotate with 1000G membership
    for line in open(args.samples_from_1000g).readlines():
        sid = line.rstrip()
        if sid in samples.keys():
            samples[sid]['1000G'] = True

    # Annotate with ancestry
    with open(args.ancestry_labels) as fin:
        for sid, ancestry in csv.reader(fin, delimiter='\t'):
            if sid in samples.keys():
                samples[sid]['pop'] = ancestry

    # Annotate with phenotype
    with open(args.phenotype_labels) as fin:
        for sid, pheno in csv.reader(fin, delimiter='\t'):
            if sid in samples.keys():
                samples[sid]['pheno'] = pheno

    # Prune known duplicate samples
    samples = prune_duplicates(samples)
    summarize_cohort(samples, 
                     stage='Obvious duplicates removed based on identical sample IDs')

    # Identify & prune remaining cryptic duplicate samples
    samples = add_cryptic_duplicates(samples, args.predicted_relatives)
    samples = prune_duplicates(samples)
    summarize_cohort(samples, 
                     stage='Cryptic duplicates removed based on genetic relatedness')

    # Prune related individuals in the validation arm
    samples = prune_related(samples, args.kinship, 
                            require={'in_validation_vcf' : True})
    summarize_cohort(samples, 
                     stage='Related samples removed from validation cohort')

    # Downsample controls from validation cohort as needed to balance 
    # ancestry distributions among cases and controls
    samples = balance_ancestries(samples, require={'in_validation_vcf' : True})
    summarize_cohort(samples, stage='Ancestries rebalanced in case/control arm')

    # Write out samples to be retained in each cohort
    write_sample_list(samples, {'in_trio_vcf' : True},
                      args.out_prefix + '.trio_cohort_final_samples.list')
    write_sample_list(samples, {'in_validation_vcf' : True},
                      args.out_prefix + '.validation_cohort_final_samples.list')



if __name__ == '__main__':
    main()

