#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023 Ryan L. Collins, Riaz Gillani and the Van Allen Laboratory
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Combine frequency information for both case/control & trio cohorts for PedSV study
"""


import argparse
import copy
import csv
import numpy as np
import pandas as pd
import pysam
import re
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
              'CN_NONREF_COUNT' : '##INFO=<ID={}CN_NONREF_COUNT,Number=1,Type=Integer,Description="Number of {} samples with non-diploid copy states (multiallelic CNVs only).">',
              'CN_NONREF_FREQ' : '##INFO=<ID={}CN_NONREF_FREQ,Number=1,Type=Float,Description="Frequency of samples from {} with non-diploid copy states (multiallelic CNVs only).">'}
popmax_entries = {'AF' : '##INFO=<ID={}POPMAX_AF,Number=A,Type=Float,Description="Maximum allele frequency in samples from {} across populations (biallelic sites only and minimum {} samples per population).">',
                  'FREQ_HOMREF' : '##INFO=<ID={}POPMAX_FREQ_HOMREF,Number=1,Type=Float,Description="Maximum frequency of homozygous reference genotypes in samples from {} across populations (biallelic sites only and minimum {} samples per population).">',
                  'FREQ_HET' : '##INFO=<ID={}POPMAX_FREQ_HET,Number=1,Type=Float,Description="Maximum frequency of heterozygous genotypes in samples from {} across populations (biallelic sites only and minimum {} samples per population).">',
                  'FREQ_HOMALT' : '##INFO=<ID={}POPMAX_FREQ_HOMALT,Number=1,Type=Float,Description="Maximum frequency of homozygous alternate genotype frequencies in samples from {} across populations (biallelic sites only and minimum {} samples per population).">',
                  'CN_FREQ' : '##INFO=<ID={}POPMAX_CN_FREQ,Number=.,Type=Float,Description="Maximum frequency of each copy state in {} samples across populations, starting from CN=0 (multiallelic CNVs only).">',
                  'CN_NONREF_FREQ' : '##INFO=<ID={}POPMAX_CN_NONREF_FREQ,Number=1,Type=Float,Description="Maximum frequency of non-diploid copy states in {} across populations (multiallelic CNVs only).">'}


def get_categories(header):
    """
    Infer which categories need to be combined between case_control and trio cohorts
    """

    cc_afs = [h for h in header.info.keys() if 'case_control' in h and h.endswith('_AF')]
    cc_cats = [re.sub('_AF$', '', re.sub('^case_control_', '', h)) for h in cc_afs]

    trio_afs = [h for h in header.info.keys() if 'trio' in h and h.endswith('_AF')]
    trio_cats = [re.sub('_AF$', '', re.sub('^trio_', '', h)) for h in trio_afs]
    
    cats = [c + '_' for c in set(cc_cats).intersection(set(trio_cats))]

    return [''] + cats


def reformat_header(header, categories, min_n_popmax):
    """
    Reformat the input VCF header for the output VCF
    """

    # Add pooled categories for each cohort to header
    for cohort in 'trio case_control'.split():
        for key, descrip in af_entries.items():
            key_prefix = cohort + '_'
            descrip_filler = re.sub('_', ' ', cohort) + ' cohort'
            header.add_line(descrip.format(key_prefix, descrip_filler))
            for sex in 'MALE FEMALE'.split():
                sex_prefix = key_prefix + sex + '_'
                sex_descrip = ' '.join([descrip_filler, sex])
                header.add_line(descrip.format(sex_prefix, sex_descrip))

    # Add all category AF values to header
    for key_prefix in categories:

        descrip_filler = re.sub(' $', '', ' '.join(['all'] + key_prefix.split('_')))
        
        for key, descrip in af_entries.items():
            header.add_line(descrip.format(key_prefix, descrip_filler))

    # Remove all existing POPMAX annotations
    old_popmax_keys = [k for k in header.info.keys() if 'POPMAX' in k]
    for key in old_popmax_keys:
        header.info.remove_header(key)

    # Add one POPMAX field for each cohort, gnomAD, and study-wide pooled
    for prefix in [''] + 'case_control trio gnomad_v2.1_sv'.split():

        if prefix == '':
            descrip_filler = 'both cohorts'
        elif prefix == 'gnomad_v2.1_sv':
            descrip_filler = 'gnomAD-SV v2.1'
        else:
            descrip_filler = 'the ' + ' '.join(prefix.split('_')) + ' cohort'
        
        for key, descrip in popmax_entries.items():
            header.add_line(descrip.format(re.sub('^_', '', re.sub('[_]+', '_', prefix + '_')), 
                                           descrip_filler, min_n_popmax))

    return header


def is_multiallelic(record):
    """
    Reports whether a pysam.VariantRecord is an mCNV
    """

    return len(record.alleles) > 2 \
           or 'MULTIALLELIC' in record.filter.keys() \
           or record.info.get('SVTYPE') in 'CNV MCNV'.split()


def pool_freqs(record, sub_prefixes, new_prefix):
    """
    Pool a specified set of frequency annotations for a pysam.VariantRecord
    """

    is_mcnv = is_multiallelic(record)

    # Compute pooled frequencies
    if is_mcnv:
        
        CN_n, NONDIP = 0, 0
        max_cn = np.max([len(record.info[p + 'CN_COUNT']) for p in sub_prefixes])
        CN_k = [0] * max_cn

        for prefix in sub_prefixes:
            CN_n += record.info.get('_'.join([prefix + 'CN_NUMBER']), 0)
            for i in range(max_cn):
                cn_k_tmp = record.info.get('_'.join([prefix + 'CN_COUNT']))
                if len(cn_k_tmp) - 1 >= i:
                    CN_k[i] += cn_k_tmp[i]
                    if i != 2:
                        NONDIP += cn_k_tmp[i]

        record.info[new_prefix + 'CN_NUMBER'] = CN_n
        record.info[new_prefix + 'CN_COUNT'] = CN_k
        record.info[new_prefix + 'CN_NONREF_COUNT'] = NONDIP
        if CN_n > 0:
            record.info[new_prefix + 'CN_FREQ'] = tuple(np.array(CN_k) / CN_n)
            record.info[new_prefix + 'CN_NONREF_FREQ'] = NONDIP / CN_n
        else:
            record.info[new_prefix + 'CN_FREQ'] = tuple([np.NaN] * max_cn)
            record.info[new_prefix + 'CN_NONREF_FREQ'] = np.NaN

    else:

        AN, AC, n_bi, n_homref, n_het, n_homalt = 0, 0, 0, 0, 0, 0
        
        for prefix in sub_prefixes:
            AN += record.info.get('_'.join([prefix + 'AN']), 0)
            AC += record.info.get('_'.join([prefix + 'AC']), (0, ))[0]
            n_bi += record.info.get('_'.join([prefix + 'N_BI_GENOS']), 0)
            n_homref += record.info.get('_'.join([prefix + 'N_HOMREF']), 0)
            n_het += record.info.get('_'.join([prefix + 'N_HET']), 0)
            n_homalt += record.info.get('_'.join([prefix + 'N_HOMALT']), 0)
        
        record.info[new_prefix + 'AN'] = AN
        record.info[new_prefix + 'AC'] = (AC, )
        if AN > 0:
            AF = (AC / AN, )
        else:
            AF = np.NaN
        record.info[new_prefix + 'AF'] = AF
        record.info[new_prefix + 'N_BI_GENOS'] = n_bi
        record.info[new_prefix + 'N_HOMREF'] = n_homref
        record.info[new_prefix + 'N_HET'] = n_het
        record.info[new_prefix + 'N_HOMALT'] = n_homalt
        if n_bi > 0:
            freq_homref = n_homref / n_bi
            freq_het = n_het / n_bi
            freq_homalt = n_homalt / n_bi
        else:
            freq_homref, freq_het, freq_homalt = np.NaN, np.NaN, np.NaN
        record.info[new_prefix + 'FREQ_HOMREF'] = freq_homref
        record.info[new_prefix + 'FREQ_HET'] = freq_het
        record.info[new_prefix + 'FREQ_HOMALT'] = freq_homalt

    return record


def pool_freqs_per_cohort(record, cohort):
    """
    Wrapper function to pool all frequencies per trio or case/control cohort
    """

    # Infer populations to collapse
    field_prefixes = set([re.sub('AN$', '', k) \
                          for k in record.header.info.keys() \
                          if k.startswith(cohort + '_') \
                          and k.endswith('_AN') \
                          and 'MALE' not in k])
    field_prefixes.remove(cohort + '_')

    # Pool frequency information across populations for this cohort
    try:
        record = pool_freqs(record, field_prefixes, cohort + '_')
    except:
        import pdb; pdb.set_trace()

    # Pool sex-specific frequencies
    for sex in 'MALE FEMALE'.split():
        sex_field_prefixes = [p + sex + '_' for p in field_prefixes]
        record = pool_freqs(record, sex_field_prefixes, '_'.join([cohort, sex]) + '_')

    return record


def pool_freqs_for_study(record, categories):
    """
    Wrapper function to pool all frequencies across both trio & case/control cohorts
    """

    for categ in categories:
        field_prefixes = ['_'.join([c, categ]) for c in 'trio case_control'.split()]
        record = pool_freqs(record, field_prefixes, categ)

    return record


def add_popmax(record, prefix, all_pops, min_n):
    """
    Add popmax annotations for a single cohort
    """

    afs = []
    f_homrefs = []
    f_hets = []
    f_homalts = []
    f_cn_all = []
    f_nondip = []

    if is_multiallelic(record) and 'gnomad' not in prefix:

        for pop in all_pops:
            
            CN_n_field = prefix  + pop + '_CN_NUMBER'
            
            if CN_n_field not in record.info.keys():
                continue
            if record.info[CN_n_field] < min_n:
                continue

            f_cn_all.append(record.info[prefix + pop + '_CN_FREQ'])
            f_nondip.append(record.info[prefix + pop + '_CN_NONREF_FREQ'])

        if len(f_cn_all) > 0:
            max_cn = np.nanmax([len(t) for t in f_cn_all])
            f_cn = []
            for i in range(max_cn):
                f_cn_tmp = []
                for t in f_cn_all:
                    if len(t) < i + 1:
                        continue
                    f_cn_tmp.append(t[i])
                f_cn.append(np.nanmax(f_cn_tmp, initial=0))
            record.info[prefix + 'POPMAX_CN_FREQ'] = f_cn
        record.info[prefix + 'POPMAX_CN_NONREF_FREQ'] = np.nanmax(f_nondip, initial=np.NaN)
        
    else:

        for pop in all_pops:
        
            AN_field = prefix + pop + '_AN'
            
            if 'gnomad' in prefix:
                if prefix + pop + '_AF' not in record.info.keys():
                    continue
                afs.append(record.info[prefix + pop + '_AF'])

            else:
                if AN_field not in record.info.keys():
                    continue
                if record.info[AN_field] < 2 * min_n:
                    continue

                afs.append(record.info[prefix + pop + '_AF'][0])
                f_homrefs.append(record.info[prefix + pop + '_FREQ_HOMREF'])
                f_hets.append(record.info[prefix + pop + '_FREQ_HET'])
                f_homalts.append(record.info[prefix + pop + '_FREQ_HOMALT'])

        if len(afs) > 0:
            record.info[prefix + 'POPMAX_AF'] = np.nanmax(afs)
        if 'gnomad' not in prefix:
            if len(f_homrefs) > 0:
                record.info[prefix + 'POPMAX_FREQ_HOMREF'] = np.nanmax(f_homrefs)
            if len(f_hets) > 0:
                record.info[prefix + 'POPMAX_FREQ_HET'] = np.nanmax(f_hets)
            if len(f_homalts) > 0:
                record.info[prefix + 'POPMAX_FREQ_HOMALT'] = np.nanmax(f_homalts)
    
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
    parser.add_argument('--pace', help='Report pace of progress', 
                        action='store_true')
    parser.add_argument('--popmax-pops', help='Comma-separated list of population ' + 
                        'abbreviations to consider when updating POPMAX ' +
                        '[default: AMR,AFR,EAS,EUR,SAS]',
                        default='AMR,AFR,EAS,EUR,SAS', type=str)
    parser.add_argument('--min-n-popmax', default=100, type=int,
                        help='Minimum number of genotyped individuals to include ' +
                        'a population in POPMAX')
    args = parser.parse_args()

    # Open connection to input vcf
    vcf_in = pysam.VariantFile(args.vcfin)

    # Find prefixes of categories to be filled for pooled AF information
    categories = get_categories(vcf_in.header)

    # Reformat header
    header = reformat_header(vcf_in.header, categories, args.min_n_popmax)

    # Open connection to output VCF
    vcf_out = pysam.VariantFile(args.vcfout, 'w', header=header)

    # Add pooled frequencies for each variant
    k = 0
    start = time.time()
    prev = start
    for record in vcf_in.fetch():

        # First, add cross-ancestry frequencies for each of the trio & case/control cohorts
        for cohort in 'trio case_control'.split():
            record = pool_freqs_per_cohort(record, cohort)
        
        # Second, pool frequencies across both trio & case/control cohorts
        record = pool_freqs_for_study(record, categories)

        # Finally, update POPMAX annotations for overall study, each cohort individually, and gnomAD
        for prefix in [''] + 'case_control_ trio_ gnomad_v2.1_sv_'.split():
            record = add_popmax(record, prefix, set(args.popmax_pops.split(',')), 
                                args.min_n_popmax)

        vcf_out.write(record)

        # Check pace if optioned
        k += 1
        if k % 10 == 0 and args.pace:
            now = time.time()
            print('Progress: annotated {:,} variants in {:.1f} seconds (last 10 variants: {:.1f} secs)'.format(k, now - start, now - prev))
            prev = now

    # Close output VCF to clear buffer
    vcf_out.close()


if __name__ == '__main__':
    main()

