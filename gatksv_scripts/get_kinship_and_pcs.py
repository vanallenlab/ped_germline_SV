#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023 Ryan L. Collins and the Van Allen Laboratory
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Development code to run KING and PCA in Hail
"""


import argparse
import hail as hl
import pandas as pd


# Set global constants
autosomes = ['chr' + str(k+1) for k in range(22)]


def write_pcs(scores, outfile):
    """
    Write PCs to output file
    """
    pca_df = scores.to_pandas().\
                    rename(columns = {'s' : '#sample'}).\
                    set_index('#sample', drop=True)
    pca_df = pd.DataFrame(pca_df['scores'].to_list(), index=pca_df.index)
    pca_df.columns = ['PC' + str(k+1) for k in range(pca_df.shape[1])]
    pca_df.to_csv(outfile, sep='\t')


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--vcf-in', required=True, help='input VCF')
    parser.add_argument('--pcs-out', required=True, help='path to output .tsv ' + 
                        'for genetic principal components')
    parser.add_argument('--kinship-out', required=True, help='path to output .tsv ' + 
                        'for kinship/relatedness metrics')
    args = parser.parse_args()

    # Initialize Hail
    hl.init(default_reference='GRCh38')

    # Read VCF, convert to Hail MatrixTable, and re-read the MatrixTable
    mt_fname = args.vcf_in.replace('.vcf.gz', '.mt')
    hl.import_vcf(args.vcf_in, force_bgz=True).write(mt_fname, overwrite=True)
    mt = hl.read_matrix_table(mt_fname)

    # Compute allele frequencies
    mt = hl.variant_qc(mt)

    # Subset to  high-quality, biallelic, well-genotyped, autosomal variants
    mt = mt.filter_rows((mt.filters.size() == 0) & \
                        (mt.qual > 10) & \
                        (mt.variant_qc.call_rate > 0.99) & \
                        (mt.variant_qc.p_value_hwe > 10e-6) & \
                        (mt.variant_qc.AF[0] >= 0.001) & \
                        (mt.variant_qc.AF[1] >= 0.001) & \
                        (mt.variant_qc.AF[0] <= 0.999) & \
                        (mt.variant_qc.AF[1] <= 0.999), 
                        keep=True)
    mt = hl.filter_intervals(mt, [hl.parse_locus_interval(x, reference_genome='GRCh38') for x in autosomes])

    # Generate top 20 PCs based on common variants
    eigenvalues, scores, loadings = hl.hwe_normalized_pca(mt.GT, k=20)

    # Write PCs to outfile
    write_pcs(scores, args.pcs_out)

    # Generate kinship coefficients with PC-relate
    rel = hl.pc_relate(mt.GT, min_individual_maf=0.001, 
                       scores_expr=scores[mt.col_key].scores,
                       min_kinship=1/16)

    # Write kinship coefficients to output file
    rel.to_pandas().\
        rename(columns = {'i.s' : '#sample_1', 'j.s' : 'sample_2'}).\
        sort_values('kin', ascending=False).\
        to_csv(args.kinship_out, sep='\t', index=False)


if __name__ == '__main__':
    main()
