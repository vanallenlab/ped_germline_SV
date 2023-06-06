#!/usr/bin/env bash

# Code to batch PedSV v2 callset prior to GATK-SV multi-sample steps
# Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# June 2023

export CODEDIR=/Users/ryan/Desktop/Collins/VanAllen/pancancer_wgs/pancan_germline_wgs
export DATADIR=/Users/ryan/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_v2_callset_generation/v2_mega_batching

# Run batching script
$CODEDIR/scripts/gatksv_helpers/make_batches.py \
  --match-on batching_sex \
  --match-on batching_pheno \
  --batch-by wgd_score \
  --batch-by median_coverage \
  --prefix $DATADIR/PedSV.v2 \
  $DATADIR/gatk_sv_pediatric_cancers_combined_sample_data_model_updated_links_6_2_23_with_annotations_for_batching.txt
