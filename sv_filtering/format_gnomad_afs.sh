#!/usr/bin/env bash

# Study of Germline SVs in Pediatric Cancers
# Copyright (c) 2023-Present, the Van Allen laboratory at Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Format gnomaAD v4.1 non-cancer frequencies for annotation

# Set up working directory
WRKDIR=`mktemp -d`
cd $WRKDIR

# Download gnomAD non-neuro controls index
gsutil -m cp \
  gs://gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.non_neuro_controls.sites.vcf.gz.tbi \
  ./

# Authenticate GCP for remote htslib streaming
export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

# Write gnomAD VCF header for reference
bcftools view -h \
  gs://gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.non_neuro_controls.sites.vcf.gz \
> $WRKDIR/gnomad.header.vcf

# Write header for columns needed for external AF annotation workflow
echo -e "#chrom\tstart\tend\tname\tsvtype\tSVTYPE\tSVLEN\tAF\tMALE_AF\tFEMALE_AF\tAFR_AF\tAMI_AF\tAMR_AF\tASJ_AF\tEAS_AF\tFIN_AF\tMID_AF\tNFE_AF\tOTH_AF\tSAS_AF" \
> $WRKDIR/gnomad_v4.1_sv.sites.header.bed

# Use bcftools query to extract fields necessary for external AF annotation
bcftools query \
  -i 'INFO/SVTYPE != "CNV"' \
  -f '%CHROM\t%POS\t%END\t%ID\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/AF_controls_and_biobanks\t%INFO/AF_XY\t%INFO/AF_XX\t%INFO/AF_controls_and_biobanks_afr\t%INFO/AF_controls_and_biobanks_ami\t%INFO/AF_controls_and_biobanks_amr\t%INFO/AF_controls_and_biobanks_asj\t%INFO/AF_controls_and_biobanks_eas\t%INFO/AF_controls_and_biobanks_fin\t%INFO/AF_controls_and_biobanks_mid\t%INFO/AF_controls_and_biobanks_nfe\t%INFO/AF_controls_and_biobanks_rmi\t%INFO/AF_controls_and_biobanks_sas\n' \
  gs://gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.non_neuro_controls.sites.vcf.gz \
> $WRKDIR/gnomad_v4.1_sv.sites.records.bed
bcftools query \
  -i 'INFO/SVTYPE == "CNV"' \
  -f '%CHROM\t%POS\t%END\t%ID\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/CN_NONREF_FREQ_controls_and_biobanks\t%INFO/CN_NONREF_FREQ_XY\t%INFO/CN_NONREF_FREQ_XX\t%INFO/CN_NONREF_FREQ_controls_and_biobanks_afr\t%INFO/CN_NONREF_FREQ_controls_and_biobanks_ami\t%INFO/CN_NONREF_FREQ_controls_and_biobanks_amr\t%INFO/CN_NONREF_FREQ_controls_and_biobanks_asj\t%INFO/CN_NONREF_FREQ_controls_and_biobanks_eas\t%INFO/CN_NONREF_FREQ_controls_and_biobanks_fin\t%INFO/CN_NONREF_FREQ_controls_and_biobanks_mid\t%INFO/CN_NONREF_FREQ_controls_and_biobanks_nfe\t%INFO/CN_NONREF_FREQ_controls_and_biobanks_rmi\t%INFO/CN_NONREF_FREQ_controls_and_biobanks_sas\n' \
  gs://gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.non_neuro_controls.sites.vcf.gz \
>> $WRKDIR/gnomad_v4.1_sv.sites.records.bed
sort -Vk1,1 -k2,2n -k3,3n -k4,4V $WRKDIR/gnomad_v4.1_sv.sites.records.bed \
| cat $WRKDIR/gnomad_v4.1_sv.sites.header.bed - \
| bgzip -c \
> $WRKDIR/gnomad_v4.1_sv.sites.bed.gz
tabix -f $WRKDIR/gnomad_v4.1_sv.sites.bed.gz

# Note: must copy BED to gs:// bucket for annotation in Terra

# Clean up
cd ~
rm -rf $WRKDIR
unset WRKDIR