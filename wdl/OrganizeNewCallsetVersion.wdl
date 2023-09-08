#########################
#    Germline SVs in    #
#   Pediatric Cancers   #
#########################
#
# OrganizeNewCallsetVersion.wdl
#
# Reorganize a new version of the PedSV callset
#
# Copyright (c) 2023-Present, the Van Allen laboratory at Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0


version 1.0


workflow OrganizeNewCallsetVersion {
  input {
    File full_cohort_w_relatives_dense_vcf
    File full_cohort_w_relatives_dense_vcf_idx
    File full_cohort_w_relatives_sites_vcf
    File full_cohort_w_relatives_sites_vcf_idx
    File full_cohort_w_relatives_sites_bed
    File full_cohort_w_relatives_sites_bed_idx
    File full_cohort_w_relatives_ad_bed
    File full_cohort_w_relatives_ad_bed_idx

    File full_cohort_dense_vcf
    File full_cohort_dense_vcf_idx
    File full_cohort_sites_vcf
    File full_cohort_sites_vcf_idx
    File full_cohort_sites_bed
    File full_cohort_sites_bed_idx
    File full_cohort_ad_bed
    File full_cohort_ad_bed_idx

    File case_control_cohort_dense_vcf
    File case_control_cohort_dense_vcf_idx
    File case_control_cohort_sites_vcf
    File case_control_cohort_sites_vcf_idx
    File case_control_cohort_sites_bed
    File case_control_cohort_sites_bed_idx
    File case_control_cohort_ad_bed
    File case_control_cohort_ad_bed_idx

    File trio_cohort_dense_vcf
    File trio_cohort_dense_vcf_idx
    File trio_cohort_sites_vcf
    File trio_cohort_sites_vcf_idx
    File trio_cohort_sites_bed
    File trio_cohort_sites_bed_idx
    File trio_cohort_ad_bed
    File trio_cohort_ad_bed_idx

    String new_version_number
    String old_version_number
    Boolean archive_old_version = true
    String target_bucket = "gs://vanallen-pedsv-analysis"
    String pedsv_docker
  }

  call RelocateFiles as OrganizeDenseVcfs {
    input:
      full_cohort_w_relatives_file = full_cohort_w_relatives_dense_vcf,
      full_cohort_w_relatives_idx = full_cohort_w_relatives_dense_vcf_idx,
      full_cohort_file = full_cohort_dense_vcf,
      full_cohort_idx = full_cohort_dense_vcf_idx,
      case_control_cohort_file = case_control_cohort_dense_vcf,
      case_control_cohort_idx = case_control_cohort_dense_vcf_idx,
      trio_cohort_file = trio_cohort_dense_vcf,
      trio_cohort_idx = trio_cohort_dense_vcf_idx,
      new_version_number = new_version_number,
      old_version_number = old_version_number,
      archive_old_version = archive_old_version,
      target_bucket = target_bucket + "/vcfs",
      file_suffix = "wGTs.vcf.gz",
      docker = pedsv_docker
  }

  call RelocateFiles as OrganizeSitesVcfs {
    input:
      full_cohort_w_relatives_file = full_cohort_w_relatives_sites_vcf,
      full_cohort_w_relatives_idx = full_cohort_w_relatives_sites_vcf_idx,
      full_cohort_file = full_cohort_sites_vcf,
      full_cohort_idx = full_cohort_sites_vcf_idx,
      case_control_cohort_file = case_control_cohort_sites_vcf,
      case_control_cohort_idx = case_control_cohort_sites_vcf_idx,
      trio_cohort_file = trio_cohort_sites_vcf,
      trio_cohort_idx = trio_cohort_sites_vcf_idx,
      new_version_number = new_version_number,
      old_version_number = old_version_number,
      archive_old_version = archive_old_version,
      target_bucket = target_bucket + "/vcfs",
      file_suffix = "sites.vcf.gz",
      docker = pedsv_docker
  }

  call RelocateFiles as OrganizeSitesBeds {
    input:
      full_cohort_w_relatives_file = full_cohort_w_relatives_sites_bed,
      full_cohort_w_relatives_idx = full_cohort_w_relatives_sites_bed_idx,
      full_cohort_file = full_cohort_sites_bed,
      full_cohort_idx = full_cohort_sites_bed_idx,
      case_control_cohort_file = case_control_cohort_sites_bed,
      case_control_cohort_idx = case_control_cohort_sites_bed_idx,
      trio_cohort_file = trio_cohort_sites_bed,
      trio_cohort_idx = trio_cohort_sites_bed_idx,
      new_version_number = new_version_number,
      old_version_number = old_version_number,
      archive_old_version = archive_old_version,
      target_bucket = target_bucket + "/beds",
      file_suffix = "sites.bed.gz",
      docker = pedsv_docker
  }

  call RelocateFiles as OrganizeDosageBeds {
    input:
      full_cohort_w_relatives_file = full_cohort_w_relatives_ad_bed,
      full_cohort_w_relatives_idx = full_cohort_w_relatives_ad_bed_idx,
      full_cohort_file = full_cohort_ad_bed,
      full_cohort_idx = full_cohort_ad_bed_idx,
      case_control_cohort_file = case_control_cohort_ad_bed,
      case_control_cohort_idx = case_control_cohort_ad_bed_idx,
      trio_cohort_file = trio_cohort_ad_bed,
      trio_cohort_idx = trio_cohort_ad_bed_idx,
      new_version_number = new_version_number,
      old_version_number = old_version_number,
      archive_old_version = archive_old_version,
      target_bucket = target_bucket + "/beds",
      file_suffix = "allele_dosages.bed.gz",
      docker = pedsv_docker
  }

  output {}
}


task RelocateFiles {
  input {
    File full_cohort_w_relatives_file
    File full_cohort_w_relatives_idx

    File full_cohort_file
    File full_cohort_idx
    
    File case_control_cohort_file
    File case_control_cohort_idx
    
    File trio_cohort_file
    File trio_cohort_idx

    String new_version_number
    String old_version_number
    Boolean archive_old_version
    String target_bucket
    String file_suffix
    String docker
  }

  Int disk_gb = ceil(2 * size([full_cohort_w_relatives_file, full_cohort_file, case_control_cohort_file, trio_cohort_file], "GB")) + 10

  command <<<

    set -eu -o pipefail

    # Relocate all new files
    gsutil -m cp \
      ~{full_cohort_w_relatives_file} \
      "~{target_bucket}/PedSV.v~{new_version_number}.full_cohort_w_relatives_1000G.~{file_suffix}"
    gsutil -m cp \
      ~{full_cohort_w_relatives_idx} \
      "~{target_bucket}/PedSV.v~{new_version_number}.full_cohort_w_relatives_1000G.~{file_suffix}.tbi"
    gsutil -m cp \
      ~{full_cohort_file} \
      "~{target_bucket}/PedSV.v~{new_version_number}.full_cohort.analysis_samples.~{file_suffix}"
    gsutil -m cp \
      ~{full_cohort_idx} \
      "~{target_bucket}/PedSV.v~{new_version_number}.full_cohort.analysis_samples.~{file_suffix}.tbi"
    gsutil -m cp \
      ~{case_control_cohort_file} \
      "~{target_bucket}/PedSV.v~{new_version_number}.case_control_cohort.analysis_samples.~{file_suffix}"
    gsutil -m cp \
      ~{case_control_cohort_idx} \
      "~{target_bucket}/PedSV.v~{new_version_number}.case_control_cohort.analysis_samples.~{file_suffix}.tbi"
    gsutil -m cp \
      ~{trio_cohort_file} \
      "~{target_bucket}/PedSV.v~{new_version_number}.trio_cohort.analysis_samples.~{file_suffix}"
    gsutil -m cp \
      ~{trio_cohort_idx} \
      "~{target_bucket}/PedSV.v~{new_version_number}.trio_cohort.analysis_samples.~{file_suffix}.tbi"

    if [ ~{archive_old_version} == "true" ]; then
      gsutil -m mv \
        "~{target_bucket}/PedSV.v~{old_version_number}.full_cohort_w_relatives_1000G.~{file_suffix}" \
        "~{target_bucket}/PedSV.v~{old_version_number}.full_cohort_w_relatives_1000G.~{file_suffix}.tbi" \
        "~{target_bucket}/PedSV.v~{old_version_number}.full_cohort.analysis_samples.~{file_suffix}" \
        "~{target_bucket}/PedSV.v~{old_version_number}.full_cohort.analysis_samples.~{file_suffix}.tbi" \
        "~{target_bucket}/PedSV.v~{old_version_number}.case_control_cohort.analysis_samples.~{file_suffix}" \
        "~{target_bucket}/PedSV.v~{old_version_number}.case_control_cohort.analysis_samples.~{file_suffix}.tbi" \
        "~{target_bucket}/PedSV.v~{old_version_number}.trio_cohort.analysis_samples.~{file_suffix}" \
        "~{target_bucket}/PedSV.v~{old_version_number}.trio_cohort.analysis_samples.~{file_suffix}.tbi" \
        ~{target_bucket}/archive/v~{old_version_number}/
    fi

  >>>

  output {}
  
  runtime {
    cpu: 1
    memory: "1.75 GiB"
    disks: "local-disk " + disk_gb + " HDD"
    bootDiskSizeGb: 10
    docker: docker
    preemptible: 3
    maxRetries: 1
  }
}

