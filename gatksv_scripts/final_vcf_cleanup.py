#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023 Ryan L. Collins, Riaz Gillani and the Van Allen Laboratory
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Final cleanup of pediatric WGS SV callsets
"""


filters_to_info = [('BOTHSIDES_SUPPORT', 
                    'Variant has read-level support for both sides of breakpoint'),
                   ('HIGH_SR_BACKGROUND', 
                    'High number of SR splits in background samples indicating messy region'),
                   ('PESR_GT_OVERDISPERSION', 
                    'High PESR dispersion')]
new_filts = ['##FILTER=<ID=UNRELIABLE_RD_GENOTYPES,Description="This variant is ' + \
             'enriched for non-reference GTs in unreliable samples and is ' + \
             'therefore less reliable overall.">',
             '##FILTER=<ID=HG38_ALT_LOCUS,Description="This variant is at ' + \
             'least half covered by loci with alternate contigs in hg38.">',
             '##FILTER=<ID=MANUAL_FAIL,Description="This variant failed ' +
             'post hoc manual review and should not be trusted.">']



import argparse
import csv
import numpy as np
import pybedtools as pbt
import pysam
from sys import stdin, stdout


def check_duplicate(record, records_seen):
    """
    Checks if a record's coordinates and SVTYPE have already been observed
    Returns boolean (and updated records_seen if False)
    """

    svtype = record.info.get('SVTYPE')
    chrom = record.chrom
    start = record.pos
    end = record.stop

    is_dup = False
    if svtype not in records_seen.keys():
        records_seen[svtype] = {chrom : {start : set([end])}}
    elif chrom not in records_seen[svtype].keys():
        records_seen[svtype][chrom] = {start : set([end])}
    elif start not in records_seen[svtype][chrom]:
        records_seen[svtype][chrom][start] = set([end])
    elif end not in records_seen[svtype][chrom][start]:
        records_seen[svtype][chrom][start].add(end)
    else:
        is_dup = True

    return is_dup, records_seen


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
    revised_ac = 0
    for sid, sdat in record.samples.items():
        GT = sdat['GT']
        if None in GT:
            continue
        if any([a > 0 for a in GT]):
            all_nonref += 1
            if sid in bad_rd_samples \
            and sdat['EV'] == ('RD', ):
                bad_nonref += 1
            else:
                revised_ac += sum(GT)

    # Tag sample as non-PASS if at least half of all original non-ref GTs were from bad_rd_samples
    if all_nonref > 0:
        if bad_nonref / all_nonref >= 0.5:
            record.filter.add('UNRELIABLE_RD_GENOTYPES')
        
    # Mask RD-only genotypes for bad_rd_samples
    for sid in bad_rd_samples:
        if record.samples[sid]['EV'] == ('RD', ):
            record.samples[sid]['GT'] = (None, None)

    # Update AC
    record.info['AC'] = revised_ac

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
    parser.add_argument('--alt-loci-bed', help='.bed with coordinates of loci ' +
                        'overlapping alternate contigs')
    parser.add_argument('--alt-loci-frac', type=float, default=0.2, 
                        help='maximum fraction of overlap permitted with ' + 
                        '--alt-loci-frac before labeling as non-PASS')
    parser.add_argument('--sample-sex', help='two-column .tsv mapping sample IDs ' +
                        'to MALE or FEMALE. Used for handling sex chromosome NCRs.')
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
    for key, descrip in filters_to_info:
        if key in header.filters.keys():
            header.filters.remove_header(key)
            header.add_line('##INFO=<ID={},Number=0,Type=Flag,Description="{}">'.format(key, descrip))
    for key in mf_infos:
        header.info.remove_header(key)
    # header.filters.remove_header('HIGH_PCRMINUS_NOCALL_RATE')
    header.add_line('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">')
    header.add_line('##INFO=<ID=OLD_ID,Number=1,Type=String,Description="Original GATK-SV variant ID before polishing">')
    for filt in new_filts:
        header.add_line(filt)

    # Open connection to output vcf
    if args.vcf_out in '- stdout /dev/stdout':
        outvcf = pysam.VariantFile(stdout, 'w', header=header)
    else:
        outvcf = pysam.VariantFile(args.vcf_out, 'w', header=header)

    # Read list of samples with unreliable RD genotypes and intersect with VCF header
    bad_rd_samples = set()
    if args.bad_rd_samples is not None:
        with open(args.bad_rd_samples) as fin:
            bad_rd_samples = set([s.rstrip() for s in fin.readlines()])
            bad_rd_samples = bad_rd_samples.intersection(set(header.samples))

    # Read table of sample sexes and partition into male/female
    male_ids, female_ids = set(), set()
    if args.sample_sex is not None:
        with open(args.sample_sex) as tsvin:
            for sid, sex in csv.reader(tsvin, delimiter='\t'):
                if sex.upper() == 'MALE':
                    male_ids.add(sid)
                elif sex.upper() == 'FEMALE':
                    female_ids.add(sid)

    # Load alt contig bed, if optioned
    if args.alt_loci_bed is not None:
        alt_bt = pbt.BedTool(args.alt_loci_bed)
    else:
        alt_bt = None

    # Load list of variant IDs to manually fail, if optioned
    if args.fail_variants is not None:
        with open(args.fail_variants) as fin:
            fail_vids = [l.rstrip() for l in fin.readlines()]
    else:
        fail_vids = []

    # Iterate over records in invcf, clean up, and write to outvcf
    svtype_counter = {}
    records_seen = {}
    for record in invcf.fetch():

        # Skip duplicated records (defined as those with identical chrom/pos/end/svtype)
        is_dup, records_seen = check_duplicate(record, records_seen)
        if is_dup:
            continue

        # Get reused record info
        svtype = record.info.get('SVTYPE')
        svlen = record.info.get('SVLEN', 0)


        ### First, do simple variant-level evaluations ###

        # Check if record should be marked as manual fail
        if record.id in fail_vids:
            record.filter.add('MANUAL_FAIL')

        # Label small MCNVs as UNRESOLVED
        if is_multiallelic(record) and svlen < 5000:
            record.filter.add('UNRESOLVED')

        # Label very large, common, unbalanced SVs as UNRESOLVED
        if svtype in 'DEL DUP INS CPX'.split() \
        and svlen >= 500000 \
        and record.info.get('AF', (0, ))[0] >= 0.05:
            record.filter.add('UNRESOLVED')

        # Label common CTX as UNRESOLVED
        if svtype == 'CTX' \
        and record.info.get('AF', (0, ))[0] >= 0.01:
            record.filter.add('UNRESOLVED')

        # Clear all MALE/FEMALE AF annotations
        for key in mf_infos:
            if key in record.info.keys():
                record.info.pop(key)


        ### Second, edit genotypes as needed ###

        # Clean up RD genotypes
        if svtype in 'DEL DUP'.split() \
        and len(bad_rd_samples) > 0:
            record = clean_rd_genos(record, bad_rd_samples)


        ### Third, perform variant-level operations that depend on genotypes ###

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
        if is_depth_only(record):
            if record.info.get('NCR', 0) >= 0.08:
                record.filter.add('HIGH_NCR')
        elif is_artifact_deletion(record):
            if record.info.get('NCR', 0) > 1/250:
                record.filter.add('HIGH_NCR')
        elif record.info.get('NCR', 0) >= 0.1:
            record.filter.add('HIGH_NCR')

        # Relocate certain FILTERs to INFOs
        for key, descrip in filters_to_info:
            if key in record.filter.keys():
                record.info[key] = True
        old_filts = set(record.filter.keys()).difference(set([k for k, v in filters_to_info]))
        if len(old_filts) == 0:
            old_filts = set(['PASS'])
        record.filter.clear()
        for filt in old_filts:
            record.filter.add(filt)

        # Recalibrate QUAL score
        if is_multiallelic(record):
            record.qual = 99
        else:
            record.qual = recalibrate_qual(record)

        # Check for alt contig coverage if optioned
        if alt_bt is not None:
            cstr = '{}\t{}\t{}\n'.format(record.chrom, record.start, record.stop)
            cov = pbt.BedTool(cstr, from_string=True).coverage(alt_bt)[0][-1]
            if float(cov) > args.alt_loci_frac:
                record.filter.add('HG38_ALT_LOCUS')

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

