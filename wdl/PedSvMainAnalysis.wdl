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
    File ref_fai
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

  Array[Array[String]] contiglist = read_tsv(ref_fai)

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

  scatter ( contig_info in contiglist ) {
    String contig = contig_info[0]

    call GetSVsPerSample as GetTrioSVsPerSample {
      input:
        sv_bed = trio_bed,
        sv_bed_idx = trio_bed_idx,
        ad_matrix = trio_ad_matrix,
        ad_matrix_idx = trio_ad_matrix_idx,
        contig = contig,
        prefix = study_prefix + ".trio_cohort." + contig,
        docker = pedsv_r_docker
    }

    call GetSVsPerSample as GetValidationSVsPerSample {
      input:
        sv_bed = validation_bed,
        sv_bed_idx = validation_bed_idx,
        ad_matrix = validation_ad_matrix,
        ad_matrix_idx = validation_ad_matrix_idx,
        contig = contig,
        prefix = study_prefix + ".validation_cohort." + contig,
        docker = pedsv_r_docker
    }
  }

  call MergeSVsPerSample as MergeTrioSVsPerSample {
    input:
      tarballs = GetTrioSVsPerSample.per_sample_tarball,
      sample_list = trio_samples_list,
      prefix = study_prefix + ".trio_cohort",
      docker = ubuntu_docker
  }

  call MergeSVsPerSample as MergeValidationSVsPerSample {
    input:
      tarballs = GetValidationSVsPerSample.per_sample_tarball,
      sample_list = validation_samples_list,
      prefix = study_prefix + ".validation_cohort",
      docker = ubuntu_docker
  }

  call CohortSummaryPlots as TrioCohortSummaryPlots {
    input:
      bed = trio_bed,
      bed_idx = trio_bed_idx,
      prefix = study_prefix + ".trio_cohort",
      af_field = "parent_AF",
      ac_field = "parent_AC",
      docker = pedsv_r_docker
  }

  call CohortSummaryPlots as ValidationCohortSummaryPlots {
    input:
      bed = validation_bed,
      bed_idx = validation_bed_idx,
      prefix = study_prefix + ".validation_cohort",
      docker = pedsv_r_docker
  }

  call UnifyOutputs {
    input:
      tarballs = [StudyWideSummaryPlots.plots_tarball,
                  TrioCohortSummaryPlots.plots_tarball,
                  ValidationCohortSummaryPlots.plots_tarball],
      prefix = study_prefix + "_analysis_outputs",
      docker = ubuntu_docker
  }

  output {
    File plots_tarball = UnifyOutputs.merged_tarball
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

  Int disk_gb = ceil(2 * size([sample_metadata_tsv], "GB")) + 10

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
      --out-prefix StudyWideSummaryPlots/pca/~{prefix} \
      ~{sample_metadata_tsv}

    # Compress output
    tar -czvf StudyWideSummaryPlots.tar.gz StudyWideSummaryPlots
  >>>

  output {
    File plots_tarball = "StudyWideSummaryPlots.tar.gz"
  }

  runtime {
    docker: docker
    memory: "15.5 GB"
    cpu: 4
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}


# Get lists of variants per sample for a single chromosome
task GetSVsPerSample {
  input {
    File sv_bed
    File sv_bed_idx
    File ad_matrix
    File ad_matrix_idx
    String contig
    String prefix
    String docker
  }

  Int disk_gb = ceil(3 * size([ad_matrix], "GB")) + 10

  command <<<
    set -eu -o pipefail

    mkdir ~{prefix}_sample_lists/

    tabix -h ~{sv_bed} ~{contig} | bgzip -c > sv_info.~{contig}.bed.gz

    /opt/ped_germline_SV/analysis/utilities/get_SVs_per_sample.R \
      --query ~{contig} \
      --outdir ~{prefix}_sample_lists \
      sv_info.~{contig}.bed.gz \
      ~{ad_matrix}

    find ~{prefix}_sample_lists/ -name "*.list" | xargs -I {} gzip {}

    tar -czvf ~{prefix}_sample_lists.tar.gz ~{prefix}_sample_lists
  >>>

  output {
    File per_sample_tarball = "~{prefix}_sample_lists.tar.gz"
  }

  runtime {
    docker: docker
    memory: "15.5 GB"
    cpu: 4
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}


# Merge outputs of GetSVsPerSample across multiple chromosomes
task MergeSVsPerSample {
  input {
    Array[File] tarballs
    File sample_list
    String prefix
    String docker
  }

  Int disk_gb = ceil(3 * size(tarballs, "GB")) + 10

  command <<<
    set -eu -o pipefail

    mkdir ~{prefix}_SVs_per_sample

    cat ~{write_lines(tarballs)} | xargs -I {} tar -xzvf {}

    while read ID; do
      find ./ -name "$ID.svids.list.gz" \
      | xargs -I {} zcat {} \
      | sort -V \
      | gzip -c \
      > ~{prefix}_SVs_per_sample/$ID.svids.list.gz
    done < sample_list

    tar -czvf ~{prefix}_SVs_per_sample.tar.gz ~{prefix}_SVs_per_sample
  >>>

  output {
    File per_sample_tarball = "~{prefix}_SVs_per_sample.tar.gz"
  }

  runtime {
    docker: docker
    memory: "15.5 GB"
    cpu: 4
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}


# Basic summary plots for a single cohort
task CohortSummaryPlots {
  input {
    File bed
    File bed_idx
    String prefix
    
    String af_field = "AF"
    String ac_field = "AC"

    String docker
  }

  Int disk_gb = ceil(2 * size([bed], "GB")) + 10

  command <<<
    set -eu -o pipefail

    # Prep output directory
    mkdir ~{prefix}

    # Plot SV counts and sizes
    /opt/ped_germline_SV/analysis/landscape/plot_sv_counts_and_sizes.R \
      --af-field ~{af_field} \
      --ac-field ~{ac_field} \
      --out-prefix ~{prefix}/~{prefix} \
      ~{bed}

    # Compress output
    tar -czvf ~{prefix}.tar.gz ~{prefix}
  >>>

  output {
    File plots_tarball = "~{prefix}.tar.gz"
  }

  runtime {
    docker: docker
    memory: "15.5 GB"
    cpu: 4
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}


# Merge multiple tarballs
task UnifyOutputs {
  input {
    Array[File] tarballs
    String prefix
    String docker
  }

  Int disk_gb = ceil(10 * size(tarballs, "GB")) + 10

  command <<<
    set -eu -o pipefail

    # Prep output directory
    mkdir ~{prefix}

    # Untar each tarball into output directory
    while read tarball; do
      tar -xzvf $tarball -C ~{prefix}
    done < <( cat ~{write_lines(tarballs)} )

    # Compress output directory
    tar -czvf ~{prefix}.tar.gz ~{prefix}
  >>>

  output {
    File merged_tarball = "~{prefix}.tar.gz"
  }

  runtime {
    docker: docker
    memory: "1.5 GB"
    cpu: 1
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}

