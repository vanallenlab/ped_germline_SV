#########################
#    Germline SVs in    #
#   Pediatric Cancers   #
#########################
#
# AnnotateAfs.wdl
#
# Annotate allele frequencies of all variants (with optional sample pruning)
#
# Copyright (c) 2023-Present, the Van Allen laboratory at Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0


version 1.0


workflow PedSvMainAnalysis {
	input {
		File sample_metadata_tsv
    String study_prefix

		File trio_sites_vcf
		File trio_sites_vcf_idx
		File trio_bed
		File trio_bed_idx
		File trio_ad_matrix
		File trio_ad_matrix_idx
		File trio_samples_list

		File validation_sites_vcf
		File validation_sites_vcf_idx
		File validation_bed
		File validation_bed_idx
		File validation_ad_matrix
		File validation_ad_matrix_idx
		File validation_samples_list

    String pedsv_r_docker
    String ubuntu_docker = "ubuntu:latest"
	}

  call ConcatTextFiles as ConcatSampleLists {
    input:
      infiles = [trio_samples_list, validation_samples_list],
      outfile_name = study_prefix + ".both_cohorts.samples.list",
      docker = ubuntu_docker
  }
  File all_samples_list = ConcatSampleLists.outfile

	call StudyWideSummaryPlots {
    input:
      sample_metadata_tsv = sample_metadata_tsv,
      samples_list = all_samples_list,
      prefix = study_prefix,
      docker = pedsv_r_docker
  }
}


# Helper task to concatenate two or more text files
task ConcatTextFiles {
  input {
    Array[File] infiles
    String outfile_name
    String docker
  }

  Int disk_gb = ceil(3 * size(infiles, "GB")) + 10

  command <<<
    set -eu -o pipefail

    cat ~{sep=" " infiles} > ~{outfile_name}
  >>>

  output {
    File outfile = "~{outfile_name}"
  }

  runtime {
    docker: docker
    memory: "1.5 GB"
    cpu: 1
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}


# Basic study-wide summary plots
task StudyWideSummaryPlots {
  input {
    File sample_metadata_tsv
    File samples_list
    String prefix
    String docker
  }

  command <<<
    set -eu -o pipefail

    # Prep output directory
    mkdir StudyWideSummaryPlots
    for subdir in pca; do
      mkdir StudyWideSummaryPlots/$subdir
    done

    # Plot PCs colored by ancestry
    /opt/ped_germline_SV/analysis/cohort_summaries/plot_pcs.R \
      --subset-samples ~{samples_list} \
      --out-prefix StudyWideSummaryPlots/pca/~{study_prefix} \
      ~{sample_metadata_tsv}
  >>>
}