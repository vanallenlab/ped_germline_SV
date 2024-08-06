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


# Titrate different MAD cutoffs to evaluate impact for reviewer response
for global in $( seq 5 11 ); do
  # Update global cutoffs
  sed "s/8\*MAD/$global\*MAD/g" \
    $BASEDIR/ped_germline_SV/gatksv_scripts/refs/global_qc_thresholds.json \
  > ~/scratch/global_qc_thresholds.test.json

  # Update batch-specific cutoffs
  batchspec=$(( $global - 1 ))
  sed "s/7\*MAD/$batchspec\*MAD/g" \
    $BASEDIR/ped_germline_SV/gatksv_scripts/refs/batch_qc_thresholds.json \
  > ~/scratch/batch_qc_thresholds.test.json

  # Re-batch samples
  $CODEDIR/scripts/gatksv_helpers/make_batches.py \
    --match-on batching_sex \
    --match-on batching_pheno \
    --batch-by wgd_score \
    --batch-by median_coverage \
    --global-qc-cutoffs ~/scratch/global_qc_thresholds.test.json \
    --custom-qc-fail-samples \
      <( cat \
          $DATADIR/PedSV.v2.discordant_sex_info.samples.list \
          $DATADIR/PedSV.v2.ambiguous_sex_ploidy.samples.list ) \
    --batch-qc-cutoffs ~/scratch/batch_qc_thresholds.test.json \
    --prefix "PedSV.v2" \
    -o /dev/null \
    -l /dev/null \
    --batch-membership-tsv ~/scratch/batching_test.assignments.tsv \
    $DATADIR/gatk_sv_pediatric_cancers_combined_sample_data_model_updated_links_6_2_23_with_annotations_for_batching.txt

  # Print total number of starting samples
  sed '1d' \
    $DATADIR/gatk_sv_pediatric_cancers_combined_sample_data_model_updated_links_6_2_23_with_annotations_for_batching.txt \
  | wc -l

  # Count samples retained
  cat ~/scratch/batching_test.assignments.tsv | wc -l

  # Compare list of samples to those in final v2.5.3 cohort
  sed '1d' ~/scratch/batching_test.assignments.tsv | cut -f2 \
  | fgrep -wvf - $DATADIR/../v2.5.3/PedSV.v2.5.3.all_samples_in_final_vcf.list \
  | wc -l
done | paste - - -
