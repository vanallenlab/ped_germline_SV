#!/usr/bin/env bash

# Study of Germline SVs in Pediatric Cancers
# Copyright (c) 2023-Present, the Van Allen laboratory at Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Infer relatedness and ancestry for both cohorts in pediatric SV study

# Set paths and global parameters
export polished_vcf="gs://val-ped-nonterra-hail/PedSV.v2.1.polished_noCodingOutliers.vcf.gz"
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
  --vcf-in $polished_vcf \
  --pcs-out gs://val-ped-nonterra-hail/PedSV.v2.1.polished.PCs.tsv.gz \
  --kinship-out gs://val-ped-nonterra-hail/PedSV.v2.1.polished.kinship.tsv.gz

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
  gs://val-ped-nonterra-hail/PedSV.v2.1.polished.*.tsv.gz \
  $WRKDIR/data/ancestry_and_relatedness/

# Assign ancestries using top 10 PCs
if [ -e $WRKDIR/data/ancestry_and_relatedness/v2.1 ]; then
  rm -rf $WRKDIR/data/ancestry_and_relatedness/v2.1
fi
mkdir $WRKDIR/data/ancestry_and_relatedness/v2.1
${CODEDIR}/gatksv_scripts/assign_ancestry.R \
  --PCs $WRKDIR/data/ancestry_and_relatedness/PedSV.v2.1.polished.PCs.tsv.gz \
  --training-labels $WRKDIR/data/ancestry_and_relatedness/1000G_HGDP_MESA_training_labels.tsv.gz \
  --testing-labels ${WRKDIR}/data/ancestry_and_relatedness/PedSV.SNV_ancestry_labels.tsv.gz \
  --out-prefix $WRKDIR/data/ancestry_and_relatedness/v2.1/PedSV.v2.1.polished \
  --use-N-PCs 10 \
  --min-probability 0 \
  --plot

# Generate list of related samples to prune from analysis-ready VCFs
$CODEDIR/gatksv_scripts/parse_kinship.R \
  $WRKDIR/data/ancestry_and_relatedness/PedSV.v2.1.polished.kinship.tsv.gz \
  $WRKDIR/PedSV_v2_callset_generation/v2_mega_batching/gatk_sv_pediatric_cancers_combined_sample_data_model_updated_links_6_2_23_with_annotations_for_batching.w_batch_assignments.tsv.gz \
> $WRKDIR/data/ancestry_and_relatedness/v2.1/PedSV.v2.1.samples_to_exclude_for_analysis_vcfs.list


