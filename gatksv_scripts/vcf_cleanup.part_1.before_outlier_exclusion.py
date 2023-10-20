#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023 Ryan L. Collins, Riaz Gillani, Jett Crowdis and the Van Allen Laboratory
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
First of two passes for final cleanup of pediatric WGS SV callsets
Part 1: prior to outlier sample exclusion
"""


filters_to_info = [('BOTHSIDES_SUPPORT', 
                    'Variant has read-level support for both sides of breakpoint'),
                   ('HIGH_SR_BACKGROUND', 
                    'High number of SR splits in background samples indicating messy region'),
                   ('PESR_GT_OVERDISPERSION', 
                    'High PESR dispersion')]
new_infos = ['##INFO=<ID=HG38_REF_PATCH_LOCUS,Number=0,Type=Flag,Description="This ' + \
             'variant is at least 20% covered by reference fix patch loci contigs">',
             '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
             '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">',
             '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">',
             '##INFO=<ID=NCR_TMP,Number=1,Type=Float,Description="Provisional rate of no-call GTs prior to final polishing">']
new_filts = ['##FILTER=<ID=HG38_ALT_LOCUS,Description="This variant is at ' + \
             'least 20% covered by loci with alternate contigs">',
             '##FILTER=<ID=PREDICTED_GRIP_JXN,Description="This variant is ' + \
             'predicted to mark a splice junction for a gene retrocopy ' + \
             'insertion event and should not be evaluated as a canonical deletion.">',
             '##FILTER=<ID=LOW_SL_MEAN,Description="This variant had an ' + \
             'abundance of low-quality genotypes prior to any filtering, ' + \
             'which usually indicates a lower-quality SV locus.">']



import argparse
import numpy as np
import pandas as pd
import pybedtools as pbt
import pysam
from sys import stdin, stdout


def make_intron_bed(gtf_in):
    """
    Build a map of all introns and return as pbt.BedTool
    """

    tx_info = {}
    prev_chrom = None
    intron_strs = ''
    ic_fmt = '{}\t{}\t{}\t{}\n'

    for record in pbt.BedTool(gtf_in):
        if record.fields[2] != 'exon':
            continue

        txid = record.attrs.get('transcript_id', None)
        if txid is None:
            continue

        # Check if we have moved to a new chromosome, in which case all old info
        # should be cleared from memory
        if prev_chrom is not None:
            if record.chrom != prev_chrom:
                tx_info = {}
        prev_chrom = record.chrom


        # If this transcript has been seen before, infer intron coordinates and add to intron_strs
        if txid in tx_info.keys():
            intron_coords = ic_fmt.format(record.chrom, tx_info[txid], 
                                          record.start, txid)
            intron_strs += intron_coords

        # Update latest exon junction
        tx_info[txid] = record.end

    return pbt.BedTool(intron_strs, from_string=True)


def intron_check(record, introns, ro=0.90):
    """
    Check if a record has reciprocal overlap with any introns
    """

    svc_fmt = '{}\t{}\t{}\n'
    rec_int = svc_fmt.format(record.chrom, record.start, record.stop)
    rec_bt = pbt.BedTool(rec_int, from_string=True)
    if len(introns.intersect(rec_bt, r=True, f=ro)) > 0:
        return True
    else:
        return False


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


def update_af(record):
    """
    Update AC, AN, AF, and NCR for a single record
    """

    ac, an = 0, 0
    af = None
    n_samples = len(record.samples.keys())
    for sdat in record.samples.values():
        GT = [a for a in sdat['GT'] if a is not None]
        an += len(GT)
        ac += len([a for a in GT if a > 0])
    if an > 0:
        af = ac / an
    record.info['AC'] = ac
    record.info['AN'] = an
    record.info['AF'] = af
    record.info['NCR_TMP'] = 1 - (an / (2 * n_samples))

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


def is_manta_andor_wham(record):
    """
    Check whether an SV was contributed by Manta and/or wham (but no other algorithms)
    """

    algs = sorted(record.info.get('ALGORITHMS'))

    if algs == ['manta'] \
    or algs == ['manta', 'wham'] \
    or algs == ['wham']:
        return True
    else:
        return False


def mask_gts_by_ogq(record):
    """
    Nullify GTs that have only PE- or SR- evidence and OGQ = 0
    """

    for sid, sdat in record.samples.items():
        OGQ = sdat.get('OGQ', 99)
        if OGQ > 0:
            continue
        EV = sdat.get('EV', tuple())
        if EV == ('SR',) or EV == ('PE',) and OGQ:
            record.samples[sid]['GT'] = (None, None)
    
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
    parser.add_argument('--gtf', required=True, help='input .gtf')
    parser.add_argument('--alt-loci-bed', help='.bed with coordinates of loci ' +
                        'with alternate contigs; these will be masked as non-PASS')
    parser.add_argument('--ref-patch-loci-bed', help='.bed with coordinates of loci ' +
                        'where the reference has been patched; these will be ' + 
                        'annotated as such in INFO but will not have FILTER changed.')
    parser.add_argument('--exclude-loci-frac', type=float, default=0.2, 
                        help='maximum fraction of overlap permitted with ' + 
                        '--alt-loci-bed or --ref-patch-loci-bed before being ' + 
                        'marked as overlapping')
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
    for info in new_infos:
        header.add_line(info)
    for filt in new_filts:
        header.add_line(filt)

    # Build map of all introns
    introns = make_intron_bed(args.gtf)

    # Load alt contig bed, if optioned
    if args.alt_loci_bed is not None:
        alt_bt = pbt.BedTool(args.alt_loci_bed)
    else:
        alt_bt = None

    # Load ref fix patch bed, if optioned
    if args.ref_patch_loci_bed is not None:
        patch_bt = pbt.BedTool(args.ref_patch_loci_bed)
    else:
        patch_bt = None

    # Open connection to output vcf
    if args.vcf_out in '- stdout /dev/stdout':
        outvcf = pysam.VariantFile(stdout, 'w', header=header)
    else:
        outvcf = pysam.VariantFile(args.vcf_out, 'w', header=header)

    # Iterate over records in invcf, clean up, and write to outvcf
    records_seen = {}
    for record in invcf.fetch():

        # Skip duplicated records (defined as those with identical chrom/pos/end/svtype)
        is_dup, records_seen = check_duplicate(record, records_seen)
        if is_dup:
            continue

        # Get reused record info
        svtype = record.info.get('SVTYPE')
        svlen = record.info.get('SVLEN', 0)

        # Label small MCNVs as UNRESOLVED
        if is_multiallelic(record) and svlen < 5000:
            record.filter.add('UNRESOLVED')

        # Mask GTs with OGQ = 0 for rare-ish variants that are 
        # from Manta and/or Wham and sample has single form of evidence
        if is_manta_andor_wham(record) \
        and record.info.get('AF', [1])[0] < 0.05 \
        and not is_multiallelic(record):
            record = mask_gts_by_ogq(record)

        # Update AC/AN/AF/NCR
        if not is_multiallelic(record):
            record = update_af(record)

        # Skip empty records
        if record.info.get('AC', (1, ))[0] == 0 \
        and not is_multiallelic(record):
            continue

        # Label very large, common, unbalanced SVs as UNRESOLVED
        if svtype in 'DEL DUP INS CPX'.split() \
        and svlen >= 500000 \
        and record.info.get('AF', (0, ))[0] >= 0.05:
            record.filter.add('UNRESOLVED')

        # Label common CTX as UNRESOLVED
        if svtype == 'CTX' \
        and record.info.get('AF', (0, ))[0] >= 0.01:
            record.filter.add('UNRESOLVED')

        # Label suspected GRIP junctions
        if record.info['SVTYPE'] == 'DEL':
            if 'RD' not in record.info['EVIDENCE']:
                if intron_check(record, introns):
                    record.filter.add('PREDICTED_GRIP_JXN')

        # Apply very targeted FILTER based on SL_MEAN + other factors
        SLMean = record.info.get('SL_MEAN', 100)
        SLMax = record.info.get('SL_MAX', 100)
        if SLMean is not None and SLMax is not None:
            if SLMean < 0 \
            and is_manta_andor_wham(record) \
            and record.info.get('SVLEN', 10e10) < 1000 \
            and not 'BAF' in record.info.get('EVIDENCE', tuple()) \
            and record.info.get('AF', [1])[0] < 0.05 \
            and record.info.get('NCR_TMP', 0) > 1 / 1000 \
            and not record.info.get('PESR_GT_OVERDISPERSION') \
            and SLMax < 75:
                record.filter.add('LOW_SL_MEAN')
        
        # Clear all MALE/FEMALE AF annotations
        for key in mf_infos:
            if key in record.info.keys():
                record.info.pop(key)

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

        # Check for alt contig coverage if optioned
        if alt_bt is not None:
            cstr = '{}\t{}\t{}\n'.format(record.chrom, record.start, record.stop)
            cov = pbt.BedTool(cstr, from_string=True).coverage(alt_bt)[0][-1]
            if float(cov) > args.exclude_loci_frac:
                record.filter.add('HG38_ALT_LOCUS')

        # Check for ref fix patch coverage if optioned
        if patch_bt is not None:
            cstr = '{}\t{}\t{}\n'.format(record.chrom, record.start, record.stop)
            cov = pbt.BedTool(cstr, from_string=True).coverage(patch_bt)[0][-1]
            if float(cov) > args.exclude_loci_frac:
                record.info['HG38_REF_PATCH_LOCUS'] = True

        # Write to outvcf
        outvcf.write(record)

    # Close connection to output file to clear buffer
    outvcf.close()


if __name__ == '__main__':
    main()

