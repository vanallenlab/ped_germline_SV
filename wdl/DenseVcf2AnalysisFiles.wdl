#########################
#    Germline SVs in    #
#   Pediatric Cancers   #
#########################
#
# DenseVcf2AnalysisFiles.wdl
#
# Convert a dense (genotyped) VCF to slimmer files for downstream analysis
#
# Copyright (c) 2023-Present, the Van Allen laboratory at Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0


version 1.0


import "https://raw.githubusercontent.com/vanallenlab/pancan_germline_wgs/main/wdl/ExcludeSamplesFromVcf.wdl" as exclude


workflow DenseVcf2AnalysisFiles {
  input {
    File vcf
    File vcf_idx
    String prefix
    File? exclude_samples_list
    String pedsv_docker
    String pedsv_r_docker
  }

  if(defined(exclude_samples_list)){
    call exclude.ExcludeSamples as ExcludeSamples {
      input:
        vcf = vcf,
        exclude_samples_list = exclude_samples_list,
        update_info = false,
        docker = pedsv_docker
    }
  }

  File upstream_vcf = select_first([ExcludeSamples.out_vcf, vcf])
  File upstream_vcf_idx = select_first([ExcludeSamples.out_vcf_idx, vcf_idx])

  call MakeSitesVcf {
    input:
      vcf = upstream_vcf,
      docker = pedsv_docker
  }

  call Vcf2Bed {
    input:
      vcf = upstream_vcf,
      docker = pedsv_docker
  }

  call MakeAdMatrix {
    input:
      vcf = upstream_vcf,
      vcf_idx = upstream_vcf_idx,
      docker = pedsv_r_docker
  }

  output {
    File dense_vcf = upstream_vcf
    File dense_vcf_idx = upstream_vcf_idx
    File sites_vcf = MakeSitesVcf.vcf_out
    File sites_vcf_idx = MakeSitesVcf.vcf_out_idx
    File slim_bed = Vcf2Bed.bed
    File slim_bed_idx = Vcf2Bed.bed_idx
    File ad_matrix = MakeAdMatrix.ad_matrix
    File ad_matrix_idx = MakeAdMatrix.ad_matrix_idx
  }
}


task Vcf2Bed {
  input {
    File vcf
    String docker
  }

  String out_filename = basename(vcf, ".vcf.gz") + ".bed.gz"

  command <<<

    set -eu -o pipefail

    svtk vcf2bed \
      --no-samples \
      --info ALL \
      --include-filters \
      ~{vcf} - \
    | bgzip -c > ~{out_filename}

    tabix -p bed -f ~{out_filename}

  >>>

  output {
    File bed = "~{out_filename}"
    File bed_idx = "~{out_filename}.tbi"
  }
  
  runtime {
    cpu: 2
    memory: "7.5 GiB"
    disks: "local-disk " + ceil(3 * size(vcf, "GB")) + 20 + " HDD"
    bootDiskSizeGb: 10
    docker: docker
    preemptible: 3
    maxRetries: 1
  }
}


task MakeSitesVcf {
  input {
    File vcf
    String docker
  }

  String out_filename = basename(vcf, ".vcf.gz") + ".sites.vcf.gz"

  command <<<

    set -eu -o pipefail

    zcat ~{vcf} | cut -f1-8 | bgzip -c > ~{out_filename}

    tabix -p vcf -f ~{out_filename}

  >>>

  output {
    File bed = "~{out_filename}"
    File bed_idx = "~{out_filename}.tbi"
  }
  
  runtime {
    cpu: 2
    memory: "7.5 GiB"
    disks: "local-disk " + ceil(2 * size(vcf, "GB")) + 20 + " HDD"
    bootDiskSizeGb: 10
    docker: docker
    preemptible: 3
    maxRetries: 1
  }
}


task MakeAdMatrix {
  input {
    File vcf
    File vcf_idx
    String docker
  }

  String out_filename = basename(vcf, ".vcf.gz") + ".allele_dosages.bed.gz"

  command <<<

    set -eu -o pipefail

    # Write header
    tabix -H ~{vcf} | fgrep -v "##" | cut -f10- \
    | paste <( echo -e "#chr\tstart\tend\tID" ) - \
    > header.bed

    # Convert AD matrix for biallelic variants
    bcftools query \
      --exclude 'FILTER ~ "MULTIALLELIC" | INFO/SVTYPE == "CNV"' \
      --format '%CHROM\t%POS\t%INFO/END\t%ID\t[%GT\t]\n' \
      ~{vcf} \
    | sed -e 's/\.\/\./NA/g' -e 's/0\/0/0/g' -e 's/0\/1/1/g' -e 's/1\/1/2/g' \
    > biallelic.ad.bed

    # Convert & normalize AD matrix for mCNVs
    bcftools query \
      --include 'FILTER ~ "MULTIALLELIC" | INFO/SVTYPE == "CNV"' \
      --format '%CHROM\t%POS\t%INFO/END\t%ID\t[%RD_CN\t]\n' \
      ~{vcf} \
    | sed -e 's/\t\./\tNA/g' \
    | /opt/ped_germline_SV/gatksv_scripts/normalize_mcnv_ad_values.R \
    > multiallelic.ad.bed

    # Merge, sort, compress, and index
    cat biallelic.ad.bed multiallelic.ad.bed \
    | sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
    | cat header.bed - \
    | bgzip -c \
    > ~{out_filename}
    tabix -p bed -f ~{out_filename}

  >>>

  output {
    File ad_matrix = "~{out_filename}"
    File ad_matrix_idx = "~{out_filename}.tbi"
  }
  
  runtime {
    cpu: 2
    memory: "3.75 GiB"
    disks: "local-disk " + ceil(10 * size(vcf, "GB")) + 20 + " HDD"
    bootDiskSizeGb: 10
    docker: docker
    preemptible: 3
    maxRetries: 1
  }
}

