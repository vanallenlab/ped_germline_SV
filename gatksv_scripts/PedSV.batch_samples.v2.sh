#!/usr/bin/env bash

# Code to batch PedSV v2 callset prior to GATK-SV multi-sample steps
# Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# June 2023

export CODEDIR=/Users/ryan/Desktop/Collins/VanAllen/pancancer_wgs/pancan_germline_wgs
export BASEDIR=/Users/ryan/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab
export DATADIR=/Users/ryan/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_v2_callset_generation/v2_mega_batching

# Run batching script
$CODEDIR/scripts/gatksv_helpers/make_batches.py \
  --match-on batching_sex \
  --match-on batching_pheno \
  --batch-by wgd_score \
  --batch-by median_coverage \
  --global-qc-cutoffs $BASEDIR/ped_germline_SV/gatksv_scripts/refs/global_qc_thresholds.json \
  --custom-qc-fail-samples \
    <( cat \
        $DATADIR/PedSV.v2.discordant_sex_info.samples.list \
        $DATADIR/PedSV.v2.ambiguous_sex_ploidy.samples.list ) \
  --batch-qc-cutoffs $BASEDIR/ped_germline_SV/gatksv_scripts/refs/batch_qc_thresholds.json \
  --prefix PedSV.v2 \
  --outfile $DATADIR/gatk_sv_pediatric_cancers_combined_sample_data_model_updated_links_6_2_23_with_annotations_for_batching.w_batch_assignments.tsv.gz \
  --batch-names-tsv $DATADIR/PedSV.v2.sample_set_entity.tsv \
  --batch-membership-tsv $DATADIR/PedSV.v2.sample_set_membership.tsv \
  --logfile $DATADIR/PedSV.v2.batching_and_qc.log \
  $DATADIR/gatk_sv_pediatric_cancers_combined_sample_data_model_updated_links_6_2_23_with_annotations_for_batching.txt


