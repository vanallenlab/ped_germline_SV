#!/usr/bin/env bash

# Study of Germline SVs in Pediatric Cancers
# Copyright (c) 2023-Present, the Van Allen laboratory at Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Infer relatedness and ancestry for both cohorts in pediatric SV study

# Set paths and global parameters
export trio_vcf="gs://fc-eefa28f4-a75d-4bbf-80d7-f37f98f4bb19/minGQ_v7_FDR2pct_NCR10pct.no_outliers.cleaned.vcf.gz"
export validation_vcf="gs://fc-47352478-ddb0-4040-802e-fc3349c50fdb/minGQ_v1_trioCutoffs_FDR2pct_NCR10pct.no_outliers.cleaned.vcf.gz"
export WRKDIR=~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab
export CODEDIR=$WRKDIR/ped_germline_SV/

# Switch to pediatric SV billing project
old_gcloud_project=$( gcloud config list project | fgrep "project = " | cut -f2 -d\= | gsed 's/\ //g' )
gcloud config set project vanallen-nih-ped

# Deploy Hail dataproc cluster
hailctl dataproc start \
  --region us-east1 \
  --num-secondary-workers 25 \
  --packages pandas \
  --max-idle 60m \
  pedsv

# Generate PCs and kinship metrics using Hail
for vcf in $trio_vcf $validation_vcf; do
  hailctl dataproc submit pedsv \
    $CODEDIR/gatksv_scripts/get_kinship_and_pcs.py \
    --vcf-in $vcf \
    --pcs-out $( echo $vcf | gsed 's/.vcf.gz$/.PCs.tsv.gz/g' ) \
    --kinship-out $( echo $vcf | gsed 's/.vcf.gz$/.kinship.tsv.gz/g' )
done

# Shut down Hail dataproc cluster
hailctl dataproc stop pedsv

# Switch back to old billing project
gcloud config set project $old_gcloud_project

# Download outputs
for subdir in data data/ancestry_and_relatedness; do
  if [ ! -e $WRKDIR/$subdir ]; then
    mkdir $WRKDIR/$subdir
  fi
done
for vcf in $trio_vcf $validation_vcf; do
  echo "$( echo $vcf | gsed 's/.vcf.gz$/.PCs.tsv.gz/g' )"
  echo "$( echo $vcf | gsed 's/.vcf.gz$/.kinship.tsv.gz/g' )"
done | gsutil -m cp -I $WRKDIR/data/ancestry_and_relatedness/

# Assign ancestries using top 4 PCs
for cohort in trio validation; do
  if [ -e $WRKDIR/data/ancestry_and_relatedness/$cohort ]; then
    rm -rf $WRKDIR/data/ancestry_and_relatedness/$cohort
  fi
  mkdir $WRKDIR/data/ancestry_and_relatedness/$cohort
done
${CODEDIR}/gatksv_scripts/assign_ancestry.R \
  --PCs $WRKDIR/data/ancestry_and_relatedness/minGQ_v7_FDR2pct_NCR10pct.no_outliers.cleaned.PCs.tsv.gz \
  --training-labels $WRKDIR/data/ancestry_and_relatedness/1000G_HGDP_training_labels.tsv.gz \
  --out-prefix $WRKDIR/data/ancestry_and_relatedness/trio/PedSV.trio_cohort \
  --use-N-PCs 4 \
  --min-probability 0.5 \
  --plot
${CODEDIR}/gatksv_scripts/assign_ancestry.R \
  --PCs $WRKDIR/data/ancestry_and_relatedness/minGQ_v1_trioCutoffs_FDR2pct_NCR10pct.no_outliers.cleaned.PCs.tsv.gz \
  --training-labels $WRKDIR/data/ancestry_and_relatedness/1000G_HGDP_training_labels.tsv.gz \
  --out-prefix $WRKDIR/data/ancestry_and_relatedness/validation/PedSV.validation_cohort \
  --use-N-PCs 4 \
  --min-probability 0.5 \
  --plot


# Assign pairwise relationships
# TODO: implement this
