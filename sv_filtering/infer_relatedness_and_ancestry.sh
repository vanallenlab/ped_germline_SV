#!/usr/bin/env bash

# Study of Germline SVs in Pediatric Cancers
# Copyright (c) 2023-Present, the Van Allen laboratory at Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Infer relatedness and ancestry for both cohorts in pediatric SV study

# Set paths and global parameters
export merged_vcf="gs://val-ped-nonterra-hail/PedSV_merged_cohorts.merged.vcf.gz"
export WRKDIR=~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab
export CODEDIR=$WRKDIR/ped_germline_SV/

# Switch to pediatric SV billing project
old_gcloud_project=$( gcloud config list project | fgrep "project = " | cut -f2 -d\= | gsed 's/\ //g' )
gcloud config set project vanallen-nih-ped-nonterra

### Merge genotypes from matching sites across two cohorts
# NOTE: This needs to be run in Terra using the MergeSvsBetweenCohorts.wdl method

# Deploy Hail dataproc cluster
hailctl dataproc start \
  --region us-central1 \
  --zone us-central1-a \
  --bucket val-ped-nonterra-hail \
  --num-secondary-workers 40 \
  --packages pandas \
  --max-idle 60m \
  pedsv

# Generate PCs and kinship metrics using Hail
hailctl dataproc submit \
  --region us-central1 \
  pedsv \
  $CODEDIR/gatksv_scripts/get_kinship_and_pcs.py \
  --vcf-in $merged_vcf \
  --pcs-out gs://val-ped-nonterra-hail/PedSV.merged.PCs.tsv.gz \
  --kinship-out gs://val-ped-nonterra-hail/PedSV.merged.kinship.tsv.gz

# Shut down Hail dataproc cluster
hailctl dataproc stop pedsv --region us-central1 

# Switch back to old billing project
gcloud config set project $old_gcloud_project

# Download outputs
for subdir in data data/ancestry_and_relatedness; do
  if [ ! -e $WRKDIR/$subdir ]; then
    mkdir $WRKDIR/$subdir
  fi
done
gsutil -m cp \
  gs://val-ped-nonterra-hail/PedSV.merged.*.tsv.gz \
  $WRKDIR/data/ancestry_and_relatedness/

# Assign ancestries using top 4 PCs
if [ -e $WRKDIR/data/ancestry_and_relatedness/merged ]; then
  rm -rf $WRKDIR/data/ancestry_and_relatedness/merged
fi
mkdir $WRKDIR/data/ancestry_and_relatedness/merged
${CODEDIR}/gatksv_scripts/assign_ancestry.R \
  --PCs $WRKDIR/data/ancestry_and_relatedness/PedSV.merged.PCs.tsv.gz \
  --training-labels $WRKDIR/data/ancestry_and_relatedness/1000G_HGDP_MESA_training_labels.tsv.gz \
  --testing-labels ${WRKDIR}/data/ancestry_and_relatedness/PedSV.SNV_ancestry_labels.tsv.gz \
  --out-prefix $WRKDIR/data/ancestry_and_relatedness/merged/PedSV.merged \
  --use-N-PCs 4 \
  --min-probability 0 \
  --plot

# Visualize relationship labels from kinship metrics
if [ -e $WRKDIR/data/ancestry_and_relatedness/merged ]; then
  rm -rf $WRKDIR/data/ancestry_and_relatedness/merged
fi
mkdir $WRKDIR/data/ancestry_and_relatedness/merged
${CODEDIR}/gatksv_scripts/infer_relatives.R \
  --metrics $WRKDIR/data/ancestry_and_relatedness/PedSV.merged.kinship.tsv.gz \
  --training-labels $WRKDIR/data/ancestry_and_relatedness/PedSV.merged.known_relationships.tsv.gz \
  --out-prefix $WRKDIR/data/ancestry_and_relatedness/merged/PedSV.merged \
  --plot

# Prep files for final sample inclusion determination
# Combine ped files
gsutil -m cat \
  gs://fc-eefa28f4-a75d-4bbf-80d7-f37f98f4bb19/trio_cohort_ped_file_11_24_22_updated_for_sex_aneuploidies.ped.txt \
  gs://fc-47352478-ddb0-4040-802e-fc3349c50fdb/cohort_ped_file_gatk_sv_pediatric_cancers_validation_cases_with_selected_discovery_trios_2_16_23_no_header.ped.txt \
| cut -f2-4 | sort -Vk1,1 -k2,2V -k3,3V | uniq \
> $WRKDIR/data/ancestry_and_relatedness/PedSV.merged.ped_trunc
# Write list of 1000G samples
zcat ${WRKDIR}/data/ancestry_and_relatedness/1000G_HGDP_training_labels.tsv.gz \
| cut -f1 | fgrep -v "#" \
> $WRKDIR/data/ancestry_and_relatedness/1000G_samples.list
# Get lists of sample IDs present in each VCF
gsutil -m cp \
  gs://fc-eefa28f4-a75d-4bbf-80d7-f37f98f4bb19/submissions/020a7b29-eb96-4c5a-9c8f-140bc52b0340/MergeSvsBetweenCohorts/b7301f39-1666-4ad2-bede-959e1aa7395d/call-GetSamples/shard-0/cacheCopy/minGQ_v7_FDR2pct_NCR10pct.no_outliers.cleaned.samples.list \
  gs://fc-eefa28f4-a75d-4bbf-80d7-f37f98f4bb19/submissions/020a7b29-eb96-4c5a-9c8f-140bc52b0340/MergeSvsBetweenCohorts/b7301f39-1666-4ad2-bede-959e1aa7395d/call-GetSamples/shard-1/cacheCopy/minGQ_v1_trioCutoffs_FDR2pct_NCR10pct.no_outliers.cleaned.samples.list \
  $WRKDIR/data/ancestry_and_relatedness/

# Prune sample inclusion lists based on relatedness
${CODEDIR}/gatksv_scripts/optimize_sample_inclusion_lists.py \
  --kinship $WRKDIR/data/ancestry_and_relatedness/PedSV.merged.kinship.tsv.gz \
  --predicted-relatives $WRKDIR/data/ancestry_and_relatedness/merged/PedSV.merged.relatedness_labels.tsv \
  --pedfile $WRKDIR/data/ancestry_and_relatedness/PedSV.merged.ped_trunc \
  --samples-from-1000g $WRKDIR/data/ancestry_and_relatedness/1000G_samples.list \
  --ancestry-labels $WRKDIR/data/ancestry_and_relatedness/merged/PedSV.merged.ancestry_labels.tsv \
  --phenotype-labels $WRKDIR/data/PedSV.all_samples.phenotype_labels.tsv \
  --trio-callset-samples $WRKDIR/data/ancestry_and_relatedness/minGQ_v7_FDR2pct_NCR10pct.no_outliers.cleaned.samples.list \
  --validation-callset-samples $WRKDIR/data/ancestry_and_relatedness/minGQ_v1_trioCutoffs_FDR2pct_NCR10pct.no_outliers.cleaned.samples.list \
  --out-prefix $WRKDIR/data/ancestry_and_relatedness/PedSV.v1

