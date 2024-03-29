#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023 Ryan L. Collins, Riaz Gillani, Jett Crowdis and the Van Allen Laboratory
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Second of two passes for final cleanup of pediatric WGS SV callsets
Part 2: after outlier sample exclusion
"""


new_infos = ['##INFO=<ID=OLD_ID,Number=1,Type=String,Description="Original ' + \
             'GATK-SV variant ID before polishing">',
             '##INFO=<ID=FAILED_COHORT_COMPARISONS,Number=.,Type=String,Description=' + \
             '"Pairs of cohorts with significantly different frequencies">']
new_filts = ['##FILTER=<ID=MANUAL_FAIL,Description="This variant failed ' +
             'post hoc manual review and should not be trusted.">',
             '##FILTER=<ID=INTERCOHORT_HETEROGENEITY,Description="This variant ' + \
             'was genotyped at significantly different frequencies between at least ' + \
             'one pair of cohorts. Only applied to rare SVs (AF<5%). See ' + \
             'INFO:FAILED_COHORT_COMPARISONS for details.">']
infos_to_rm = ['NCR_TMP']



import argparse
import numpy as np
import pandas as pd
import pybedtools as pbt
import pysam
from scipy.stats import fisher_exact
from sys import stdin, stdout


def is_multiallelic(record):
    """
    Reports whether a pysam.VariantRecord is an mCNV
    """

    return len(record.alleles) > 2 \
           or 'MULTIALLELIC' in record.filter.keys() \
           or record.info.get('SVTYPE') in 'CNV MCNV'.split()


def is_depth_only(record):
    """
    Checks whether a record is a read depth-only CNV
    """

    if record.info['ALGORITHMS'] == ('depth', ) \
    and 'PE' not in record.info['EVIDENCE'] \
    and 'SR' not in record.info['EVIDENCE']:
        return True
    else:
        return False


def clean_rd_genos(record, bad_rd_samples):
    """
    Clean up RD genotypes for depth-only records
    """

    # Count proportion of non-ref GTs contributed by bad_rd_samples
    all_nonref = 0
    bad_nonref = 0
    for sid, sdat in record.samples.items():
        GT = sdat['GT']
        if all([a is None for a in GT]):
            continue
        if any([a > 0 for a in GT if a is not None]):
            all_nonref += 1
            if sid in bad_rd_samples \
            and sdat['EV'] == ('RD', ):
                bad_nonref += 1

    # Tag sample as non-PASS if at least half of all original non-ref GTs were from bad_rd_samples
    if all_nonref > 0:
        if bad_nonref / all_nonref >= 0.5:
            record.filter.add('UNRELIABLE_RD_GENOTYPES')
        
    # Mask RD-only genotypes for bad_rd_samples
    for sid in bad_rd_samples:
        if record.samples[sid]['EV'] == ('RD', ):
            record.samples[sid]['GT'] = (None, None)

    return record


def update_af(record):
    """
    Update AC, AN, and AF for a single record
    """

    ac, an = 0, 0
    af = None
    for sdat in record.samples.values():
        GT = [a for a in sdat['GT'] if a is not None]
        an += len(GT)
        ac += len([a for a in GT if a > 0])
    if an > 0:
        af = ac / an
    record.info['AC'] = ac
    record.info['AN'] = an
    record.info['AF'] = af

    return record


def is_artifact_deletion(record):
    """
    Checks whether a record is a deletion in the artifact zone
    """

    svtype = record.info['SVTYPE']
    svlen = record.info['SVLEN']
    if svtype == "DEL" \
    and svlen > 400 \
    and svlen < 1000 \
    and record.info.get('AC', (0, ))[0] / record.info.get('AN', 1) < 0.05:
        return True
    else:
        return False


def calc_ncr(record, exclude_samples=[]):
    """
    Calculate no-call rate for a record's genotypes
    """

    total = 0
    nocalls = 0
    for sid, sinfo in record.samples.items():
        
        if sid in exclude_samples:
            continue
        
        total += 1

        if all([a is None for a in sinfo.get('GT', (None, ))]):
            nocalls += 1

    if total > 0:
        return nocalls / total
    else:
        return None


def recalibrate_qual(record):
    """
    Recalibrate record QUAL (quality)
    Defined as median GQ among all non-reference GTs
    """

    gqs = []
    for sinfo in record.samples.values():
        a = [a for a in sinfo.get('GT', (None, None)) if a is not None]
        if np.nansum(a) > 0:
            gqs.append(sinfo.get('GQ', None))

    if len(gqs) > 0:
        return np.nanmedian(gqs)
    else:
        return 0


def count_by_cohort(record, cohort_map):
    """
    Collect AC & AN per cohort per ancestry
    """

    counts = {}

    # Iterate over all samples
    for sid, sdat in record.samples.items():
        cohort, ancestry = cohort_map.get(sid, {'cohort' : None, 'ancestry' : None}).values()

        # Only keep samples with defined cohort and ancestry
        if cohort is None or ancestry is None:
            continue
        if cohort not in counts.keys():
            counts[cohort] = {}
        if ancestry not in counts[cohort].keys():
            counts[cohort][ancestry] = {'AN' : 0, 'AC' : 0}

        # Parse sample genotype and update counter
        GT = sdat['GT']
        if all([a is None for a in GT]):
            continue
        sAN = len([a for a in GT if a is not None])
        sAC = len([a for a in GT if a is not None and a > 0])
        counts[cohort][ancestry]['AN'] += sAN
        counts[cohort][ancestry]['AC'] += sAC

    return counts


def compare_cohorts(record, counts, cohort1, cohort2, rare_af = 0.05, 
                    rare_pval = 0.05, common_pval = 0.05):
    """
    Test for frequency differences between two cohorts for a single record
    """

    # Format results for cohort1 and cohort2 as pd.DataFrames
    d1 = pd.DataFrame(counts[cohort1]).T
    d2 = pd.DataFrame(counts[cohort2]).T

    # Don't process records with no non-ref observations
    if all(d1.AC == 0) and all(d2.AC == 0):
        return record

    # Subset to populations represented in both cohort
    pops = list(set(d1.index.tolist()).intersection(set(d2.index.tolist())))
    d1.columns = ['c1' + k for k in d1.columns]
    d2.columns = ['c2' + k for k in d2.columns]
    dm = pd.concat([d1.loc[pops, :], d2.loc[pops, :]], axis=1)

    # Find population with greatest AC
    test_pop = dm.loc[:, 'c1AC c2AC'.split()].sum(axis=1).\
                  sort_values(ascending=False).index[0]

    # Compute Fisher's exact test for most informative population
    c1AN, c1AC, c2AN, c2AC = dm.loc[test_pop, :].values
    fisher_p = fisher_exact(np.array([[c1AN - c1AC, c2AN - c2AC], [c1AC, c2AC]]))[1]

    # Decide on pass/fail based on variant frequency
    AF = record.info.get('AF')[0]
    if AF < rare_af:
        cutoff = rare_pval
    else:
        cutoff = common_pval

    if fisher_p < cutoff:
        record.filter.add('INTERCOHORT_HETEROGENEITY')
        fail_label = cohort1 + '_vs_' + cohort2
        if 'FAILED_COHORT_COMPARISONS' in record.info.keys():
            record.info['FAILED_COHORT_COMPARISONS'] += (fail_label, )
        else:
            record.info['FAILED_COHORT_COMPARISONS'] = [fail_label]
    
    return record


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf_in', help='input .vcf')    
    parser.add_argument('vcf_out', help='output .vcf')
    parser.add_argument('--bad-rd-samples', help='list of samples with ' +
                        'unreliable RD genotypes')
    parser.add_argument('--fail-variants', help='list of variant IDs to be ' +
                        'marked as having failed manual review')
    parser.add_argument('--version-number', help='callset version number for ' +
                        'tagging variant IDs.', type=str)
    parser.add_argument('--sample-metadata', help='.tsv with sample metadata. ' +
                        'First column must be sample ID. Must include columns ' +
                        'named "sex", "cohort", "ancestry", "proband". Other ' +
                        'columns will be ignored.')
    args = parser.parse_args()

    # Open connection to input vcf
    if args.vcf_in in '- stdin /dev/stdin'.split():
        invcf = pysam.VariantFile(stdin)
    else:
        invcf = pysam.VariantFile(args.vcf_in)
    samples = [s for s in invcf.header.samples]

    mf_infos = [k for k in invcf.header.info.keys() \
                if k.startswith('MALE_') or k.startswith('FEMALE_')]

    # Reformat header
    header = invcf.header
    for key in mf_infos:
        header.info.remove_header(key)
    for info in new_infos:
        header.add_line(info)
    for filt in new_filts:
        header.add_line(filt)

    # Load sample metadata as pd.DataFrame and subset to samples present in VCF header
    md = pd.read_csv(args.sample_metadata, sep='\t')
    md.set_index(md[md.columns[0]], drop=True, inplace=True)
    md = md.loc[samples, :]

    # For convenience, reassign all ICGC osteo samples to StJude
    # since these are the only two cohorts that required realignment
    md.loc[md.cohort == 'ICGC', 'cohort'] = 'StJude'

    # Make list of male/female samples for handling sex chromosome NCRs    
    male_ids = set(md.index[md.sex == 'MALE'].tolist())
    female_ids = set(md.index[md.sex == 'FEMALE'].tolist())

    # Build sample-keyed dict for cohort & ancestry mappings
    # Do not include probands from trios nor samples lacking ancestry assignments
    cohort_map = md.loc[(md.proband != 'Yes') & (~md.ancestry.isna()), 
                        'cohort ancestry'.split()].\
                    to_dict(orient='index')

    # Read list of samples with unreliable RD genotypes and intersect with VCF header
    bad_rd_samples = set()
    if args.bad_rd_samples is not None:
        with open(args.bad_rd_samples) as fin:
            bad_rd_samples = set([s.rstrip() for s in fin.readlines()])
            bad_rd_samples = bad_rd_samples.intersection(set(header.samples))
    if len(bad_rd_samples) > 0:
        rd_filt = '##FILTER=<ID=UNRELIABLE_RD_GENOTYPES,Description="This ' + \
                  'variant is enriched for non-reference GTs in unreliable ' + \
                  'samples and is therefore less reliable overall.">'
        header.add_line(rd_filt)

    # Load list of variant IDs to manually fail, if optioned
    if args.fail_variants is not None:
        with open(args.fail_variants) as fin:
            fail_vids = [l.rstrip() for l in fin.readlines()]
    else:
        fail_vids = []

    # Open connection to output vcf
    if args.vcf_out in '- stdout /dev/stdout':
        outvcf = pysam.VariantFile(stdout, 'w', header=header)
    else:
        outvcf = pysam.VariantFile(args.vcf_out, 'w', header=header)

    # Iterate over records in invcf, clean up, and write to outvcf
    svtype_counter = {}
    for record in invcf.fetch():

        # Get reused record info
        svtype = record.info.get('SVTYPE')
        svlen = record.info.get('SVLEN', 0)

        # Clear unnecessary INFOs
        for key in infos_to_rm:
            if key in record.info.keys():
                record.info.pop(key)

        # Check if record should be marked as manual fail
        if record.id in fail_vids:
            record.filter.add('MANUAL_FAIL')

        # Clear all MALE/FEMALE AF annotations
        for key in mf_infos:
            if key in record.info.keys():
                record.info.pop(key)

        # Clean up RD genotypes
        if svtype in 'DEL DUP'.split() \
        and len(bad_rd_samples) > 0:
            record = clean_rd_genos(record, bad_rd_samples)
                
        # Update AC/AN/AF
        if not is_multiallelic(record):
            record = update_af(record)

        # Skip empty records
        if record.info.get('AC', (1, ))[0] == 0 \
        and not is_multiallelic(record):
            continue

        # Recompute NCR
        if not is_multiallelic(record):
            if record.chrom == 'chrX':
                record.info['NCR'] = calc_ncr(record, exclude_samples=male_ids)
            elif record.chrom == 'chrY':
                record.info['NCR'] = calc_ncr(record, exclude_samples=female_ids)
            else:
                record.info['NCR'] = calc_ncr(record)

        # Clear old high NCR FILTER and reannotate based on updated NCR
        original_filters = [k for k in record.filter.keys()]
        record.filter.clear()
        for k in original_filters:
            if k not in 'HIGH_NCR HIGH_PCRMINUS_NOCALL_RATE'.split():
                record.filter.add(k)
        if is_artifact_deletion(record):
            if record.info.get('NCR', 0) > 1/250:
                record.filter.add('HIGH_NCR')
        elif record.info.get('NCR', 0) >= 0.04:
            record.filter.add('HIGH_NCR')

        # Recalibrate QUAL score
        if is_multiallelic(record):
            record.qual = 99
        else:
            record.qual = recalibrate_qual(record)

        # Compare frequencies between case cohorts (ICGC+StJude vs. GMKF) 
        # and, separately, between control cohorts
        cpairs = [['StJude', 'GMKF'], ['Topmed_MESA', 'Topmed_BIOME']]
        if not is_multiallelic(record):
            if record.info.get('AF', (1,))[0] < 0.05:
                ac_by_cohort = count_by_cohort(record, cohort_map)
                for cpair in cpairs:
                    record = compare_cohorts(record, ac_by_cohort, 
                                             cohort1=cpair[0], cohort2=cpair[1])

        # Rename record
        if record.chrom not in svtype_counter.keys():
            svtype_counter[record.chrom] = {}
        if svtype not in svtype_counter[record.chrom].keys():
            svtype_counter[record.chrom][svtype] = 0
        svtype_counter[record.chrom][svtype] += 1
        if args.version_number is None:
            id_prefix = 'PedSV'
        else:
            id_prefix = 'PedSV.' + str(args.version_number)
        new_id = '_'.join([id_prefix, svtype, record.chrom,
                           str(svtype_counter[record.chrom][svtype])])
        record.info['OLD_ID'] = record.id
        record.id = new_id

        # Write to outvcf
        outvcf.write(record)

    # Close connection to output file to clear buffer
    outvcf.close()


if __name__ == '__main__':
    main()

