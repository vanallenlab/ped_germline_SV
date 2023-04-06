#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023 Ryan L. Collins, Riaz Gillani and the Van Allen Laboratory
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Annotate AFs in pediatric cancer SV callsets
"""


import argparse
import csv
import numpy as np
import pandas as pd
import pysam
import time
from copy import deepcopy
from itertools import product


af_entries = {'AN' : '##INFO=<ID={}AN,Number=1,Type=Integer,Description="Total number of alleles genotyped in {} samples (biallelic sites only).">',
              'AC' : '##INFO=<ID={}AC,Number=A,Type=Integer,Description="Number of non-reference alleles observed in {} samples (biallelic sites only).">',
              'AF' : '##INFO=<ID={}AF,Number=A,Type=Float,Description="Allele frequency in {} samples (biallelic sites only).">',
              'N_BI_GENOS' : '##INFO=<ID={}N_BI_GENOS,Number=1,Type=Integer,Description="Total number of {} samples with complete genotypes (biallelic sites only).">',
              'N_HOMREF' : '##INFO=<ID={}N_HOMREF,Number=1,Type=Integer,Description="Number of {} samples with homozygous reference genotypes (biallelic sites only).">',
              'N_HET' : '##INFO=<ID={}N_HET,Number=1,Type=Integer,Description="Number of {} samples with heterozygous genotypes (biallelic sites only).">',
              'H_HOMALT' : '##INFO=<ID={}N_HOMALT,Number=1,Type=Integer,Description="Number of {} samples with homozygous alternate genotypes (biallelic sites only).">',
              'FREQ_HOMREF' : '##INFO=<ID={}FREQ_HOMREF,Number=1,Type=Float,Description="Homozygous reference genotype frequency in {} samples (biallelic sites only).">',
              'FREQ_HET' : '##INFO=<ID={}FREQ_HET,Number=1,Type=Float,Description="Heterozygous genotype frequency in {} samples (biallelic sites only).">',
              'FREQ_HOMALT' : '##INFO=<ID={}FREQ_HOMALT,Number=1,Type=Float,Description="Homozygous alternate genotype frequency in {} samples (biallelic sites only).">',
              'CN_NUMBER' : '##INFO=<ID={}CN_NUMBER,Number=1,Type=Integer,Description="Total number of {} samples with estimated copy numbers (multiallelic CNVs only).">',
              'CN_COUNT' : '##INFO=<ID={}CN_COUNT,Number=.,Type=Integer,Description="Number of {} samples observed at each copy state, starting from CN=0 (multiallelic CNVs only).">',
              'CN_FREQ' : '##INFO=<ID={}CN_FREQ,Number=.,Type=Float,Description="Frequency of {} samples observed at each copy state, starting from CN=0 (multiallelic CNVs only).">',
              'CN_NONDIPLOID_COUNT' : '##INFO=<ID={}CN_NONDIPLOID_COUNT,Number=1,Type=Integer,Description="Number of {} samples with non-diploid copy states (multiallelic CNVs only).">',
              'CN_NONDIPLOID_FREQ' : '##INFO=<ID={}CN_NONDIPLOID_FREQ,Number=1,Type=Float,Description="Frequency of {} samples with non-diploid copy states (multiallelic CNVs only).">'}
popmax_entries = {'AF' : '##INFO=<ID=POPMAX_{}AF,Number=A,Type=Float,Description="Maximum allele frequency in {} samples across populations (biallelic sites only).">',
                  'FREQ_HOMREF' : '##INFO=<ID=POPMAX_{}FREQ_HOMREF,Number=1,Type=Float,Description="Maximum frequency of homozygous reference genotypes in {} samples across populations (biallelic sites only).">',
                  'FREQ_HET' : '##INFO=<ID=POPMAX_{}FREQ_HET,Number=1,Type=Float,Description="Maximum frequency of heterozygous genotypes in {} samples across populations (biallelic sites only).">',
                  'FREQ_HOMALT' : '##INFO=<ID=POPMAX_{}FREQ_HOMALT,Number=1,Type=Float,Description="Maximum frequency of homozygous alternate genotype frequencies in {} samples across populations (biallelic sites only).">',
                  'CN_FREQ' : '##INFO=<ID=POPMAX_{}CN_FREQ,Number=.,Type=Float,Description="Maximum frequency of each copy state in {} samples across populations, starting from CN=0 (multiallelic CNVs only).">',
                  'CN_NONDIPLOID_FREQ' : '##INFO=<ID=POPMAX_{}CN_NONDIPLOID_FREQ,Number=1,Type=Float,Description="Maximum frequency of non-diploid copy states in {} across populations (multiallelic CNVs only).">'}


def load_labels(sample_info, labels_in, key):
    """
    Add sample labels to sample_info
    """

    with open(labels_in) as fin:
        for sid, val in csv.reader(fin, delimiter='\t'):
            if sid in sample_info.keys():
                sample_info[sid][key] = val

    return sample_info


def expand_categories(svals):
    """
    Expand all possible combinations of designations for a single sample's labels
    """

    none_vals = tuple([None] * len(svals))

    return product(*pd.DataFrame([svals, none_vals]).transpose().values.tolist())


def enumerate_categories(sample_info, drop_pop=False):
    """
    Define all combinations of sample labels with at least one sample to be annotated
    """

    keys = list(sample_info[list(sample_info.keys())[0]].keys())
    if drop_pop:
        keys = [k for k in keys if k != 'ancestry']

    none_vals = tuple([None] * len(keys))
    
    combos = set([none_vals])

    for sid, vals in sample_info.items():
        svals = tuple([vals[k] for k in keys])
        combos.add(svals)
        for ncombo in expand_categories(svals):
            combos.add(ncombo)

    def _combo_sorter(combo):
        sort_vals = list(combo)
        for i in range(len(sort_vals)):
            if sort_vals[i] is None:
                sort_vals[i] = 'AAAA'
        return sort_vals

    return sorted(combos, key=_combo_sorter)


def get_combo_prefix(combo, prefix_type='info'):
    """
    Format the INFO prefix for a given combination of sample labels
    """

    cvals = [c for c in combo if c is not None]

    if len(cvals) == 0:
        key_prefix = ''
        descrip_filler = 'all'
    else:
        key_prefix = '_'.join(cvals) + '_'
        descrip_filler = ' '.join(cvals)

    if prefix_type == 'info':
        return key_prefix
    elif prefix_type == 'description':
        return descrip_filler


def reformat_header(header_in, categories, categories_noPop):
    """
    Reformat the input VCF header for the output VCF
    """

    # Add all category AF values to header
    for combo in categories:

        key_prefix = get_combo_prefix(combo, prefix_type='info')
        descrip_filler = get_combo_prefix(combo, prefix_type='description')
        
        for key, descrip in af_entries.items():
            header_in.add_line(descrip.format(key_prefix, descrip_filler))

    # Add popmax values to header
    for combo in categories_noPop:

        key_prefix = get_combo_prefix(combo, prefix_type='info')
        descrip_filler = get_combo_prefix(combo, prefix_type='description')
        
        for key, descrip in popmax_entries.items():
            header_in.add_line(descrip.format(key_prefix, descrip_filler))


def is_multiallelic(record):
    """
    Reports whether a pysam.VariantRecord is an mCNV
    """

    return len(record.alleles) > 2 \
           or 'MULTIALLELIC' in record.filter.keys() \
           or record.info.get('SVTYPE') in 'CNV MCNV'.split()


def categorize_gt(gt):
    """
    Categorize genotype as HOMREF, HET, or HOMALT
    """

    if gt[0] == gt[1]:
        if gt[0] == 0:
            return 'HOMREF'
        else:
            return 'HOMALT'
    else:
        return 'HET'


def annotate_freqs(record, sample_info, categories, categories_noPop, pops=[]):
    """
    Add frequency annotations to a pysam.VariantRecord
    """

    is_mcnv = is_multiallelic(record)

    if is_mcnv:
        empty_counter = {k : 0 for k in 'CN_NUMBER CN_NONDIPLOID_COUNT'.split()}
        max_cn = np.nanmax([v['CN'] for s, v in record.samples.items() if v['CN'] is not None])
        empty_counter['CN_COUNT'] = [0] * (max_cn + 1)
    else:
        empty_counter = {k : 0 for k in 'AN AC N_BI_GENOS N_HOMREF N_HET N_HOMALT'.split()}

    freqs = {categ : deepcopy(empty_counter) for categ in categories}

    # Populate an allele counter for each combination
    for sid, svals in record.samples.items():

        if is_mcnv:
            cn = svals['CN']
            if cn is None:
                continue
            if np.isnan(cn):
                continue
            for combo in expand_categories(sample_info[sid].values()):
                freqs[combo]['CN_NUMBER'] += 1
                freqs[combo]['CN_COUNT'][cn] += 1
                if cn != 2:
                    freqs[combo]['CN_NONDIPLOID_COUNT'] += 1

        else:
            called_alleles = [a is not None for a in svals['GT']]
            if not any(called_alleles):
                continue
            bi_geno = int(all(called_alleles))
            an = len([a for a in svals['GT'] if a is not None])
            ac = len([a for a in svals['GT'] if a is not None and a > 0])
            gt_label = categorize_gt(svals['GT'])

            for combo in expand_categories(sample_info[sid].values()):
                freqs[combo]['AN'] += an
                freqs[combo]['AC'] += ac
                freqs[combo]['N_BI_GENOS'] += bi_geno
                freqs[combo]['N_' + gt_label] += 1

    # Compute all frequencies
    for combo in categories:

        if is_mcnv:
            CN_NUMBER = freqs[combo]['CN_NUMBER']
            if CN_NUMBER > 0:
                CN_FREQ = [k / CN_NUMBER for k in freqs[combo]['CN_COUNT']]
                CN_NONDIPLOID_FREQ = freqs[combo]['CN_NONDIPLOID_COUNT'] / CN_NUMBER
            else:
                CN_FREQ = [None] * (max_cn + 1)
                CN_NONDIPLOID_FREQ = None
            freqs[combo]['CN_FREQ'] = CN_FREQ
            freqs[combo]['CN_NONDIPLOID_FREQ'] = CN_NONDIPLOID_FREQ

        else:
            AN = freqs[combo]['AN']
            if AN > 0:
                AF = freqs[combo]['AC'] / AN
            else:
                AF = None
            freqs[combo]['AF'] = AF
            N_BI = freqs[combo]['N_BI_GENOS']
            if N_BI > 0:
                FREQ_HOMREF = freqs[combo]['N_HOMREF'] / N_BI
                FREQ_HET = freqs[combo]['N_HET'] / N_BI
                FREQ_HOMALT = freqs[combo]['N_HOMALT'] / N_BI
            else:
                FREQ_HOMREF, FREQ_HET, FREQ_HOMALT = [None] * 3
            freqs[combo]['FREQ_HOMREF'] = FREQ_HOMREF
            freqs[combo]['FREQ_HET'] = FREQ_HET
            freqs[combo]['FREQ_HOMALT'] = FREQ_HOMALT

        # Write frequencies to INFO
        prefix = get_combo_prefix(combo, prefix_type='info')
        for key, value in freqs[combo].items():
            record.info[prefix + key] = value

    # Add POPMAX frequencies if optioned
    for combo in categories_noPop:
        prefixes = [get_combo_prefix(tuple([pop] + list(combo)), prefix_type='info') \
                    for pop in pops]
        if is_mcnv:
            popmax_keys = 'CN_FREQ CN_NONDIPLOID_FREQ'.split()
        else:
            popmax_keys = 'AF FREQ_HOMREF FREQ_HET FREQ_HOMALT'.split()
        for key in popmax_keys:
            terms = [t for t in [p + key for p in prefixes] if t in record.info.keys()]
            vals = np.array([record.info.get(t, 0) for t in terms], dtype=float)
            if key == 'CN_FREQ':
                val = tuple(np.apply_along_axis(np.nanmax, 0, vals).tolist())
            else:
                vals = vals[vals != np.array(None)]
                if len(vals) > 0:
                    val = np.nanmax(vals)
                else:
                    val = None
            record.info['POPMAX_' + get_combo_prefix(combo, prefix_type='info') + key] = val

    return record


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcfin', help='input vcf [can be bgzipped]')
    parser.add_argument('vcfout', help='output vcf [can be bgzipped]')
    parser.add_argument('--ancestry-labels', help='.tsv mapping sample IDs ' +
                        'to ancestries')
    parser.add_argument('--disease-labels', help='.tsv mapping sample IDs ' +
                        'to disease status')
    parser.add_argument('--family-labels', help='.tsv mapping sample IDs ' +
                        'to family membership (parent or proband)')
    parser.add_argument('--keep-samples', help='list of sample IDs to retain ' +
                        '[default: keep all samples]')
    parser.add_argument('--pace', help='Report pace of progress', 
                        action='store_true')
    args = parser.parse_args()

    # Open connection to input vcf
    vcf_in = pysam.VariantFile(args.vcfin)

    # Load list of samples to keep
    out_samples = set(vcf_in.header.samples)
    if args.keep_samples is not None:
        with open(args.keep_samples) as fin:
            keep_samples = set([s.rstrip() for s in fin.readlines()])
        out_samples = out_samples.intersection(keep_samples)

    # Subset input VCF to list of samples to keep
    vcf_in.subset_samples(out_samples)

    # Load dictionary of sample metadata to annotate
    sample_info = {s : {} for s in out_samples}
    if args.ancestry_labels is not None:
        sample_info = load_labels(sample_info, args.ancestry_labels, 'ancestry')
        pops = set([v['ancestry'] for v in sample_info.values()])
    else:
        pops = []
    if args.disease_labels is not None:
        sample_info = load_labels(sample_info, args.disease_labels, 'phenotype')
    if args.family_labels is not None:
        sample_info = load_labels(sample_info, args.family_labels, 'membership')

    # Enumerate all categories to annotate
    categories = enumerate_categories(sample_info)
    if args.ancestry_labels is not None:
        categories_noPop = enumerate_categories(sample_info, drop_pop=True)
    else:
        categories_noPop = []

    # Reformat header
    reformat_header(vcf_in.header, categories, categories_noPop)

    # Open connection to output VCF
    vcf_out = pysam.VariantFile(args.vcfout, 'w', header=vcf_in.header)

    # Annotate frequencies for each variant
    k = 0
    start = time.time()
    prev = start
    for record in vcf_in.fetch():
        record = annotate_freqs(record, sample_info, categories, 
                                categories_noPop, pops)
        vcf_out.write(record)
        k += 1
        if k % 10 == 0 and args.pace:
            now = time.time()
            print('Progress: annotated {:,} variants in {:.1f} seconds (last 10 variants: {:.1f} secs)'.format(k, now - start, now - prev))
            prev = now

    # Close output VCF to clear buffer
    vcf_out.close()


if __name__ == '__main__':
    main()

