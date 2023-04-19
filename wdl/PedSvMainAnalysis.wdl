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

  call BurdenTests as StudyWideBurdenTests {
    input:
      beds = [trio_bed, validation_bed],
      bed_idxs = [trio_bed_idx, validation_bed_idx],
      ad_matrixes = [trio_ad_matrix, validation_ad_matrix],
      ad_matrix_idxs = [trio_ad_matrix_idx, validation_ad_matrix_idx],
      sample_metadata_tsv = sample_metadata_tsv,
      samples_list = all_samples_list,
      af_fields = ["parent_AF", "AF"],
      ac_fields = ["parent_AC", "AC"],
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
      sample_metadata_tsv = sample_metadata_tsv,
      sample_list = trio_samples_list,
      per_sample_tarball = MergeTrioSVsPerSample.per_sample_tarball,
      af_field = "parent_AF",
      ac_field = "parent_AC",
      docker = pedsv_r_docker
  }

  call CohortSummaryPlots as ValidationCohortSummaryPlots {
    input:
      bed = validation_bed,
      bed_idx = validation_bed_idx,
      sample_metadata_tsv = sample_metadata_tsv,
      sample_list = validation_samples_list,
      per_sample_tarball = MergeValidationSVsPerSample.per_sample_tarball,
      prefix = study_prefix + ".validation_cohort",
      docker = pedsv_r_docker
  }

  call BurdenTests as TrioCohortBurdenTests {
    input:
      beds = [trio_bed],
      bed_idxs = [trio_bed_idx],
      ad_matrixes = [trio_ad_matrix],
      ad_matrix_idxs = [trio_ad_matrix_idx],
      sample_metadata_tsv = sample_metadata_tsv,
      samples_list = trio_samples_list,
      af_fields = ["parent_AF"],
      ac_fields = ["parent_AC"],
      prefix = study_prefix + ".trio_cohort",
      docker = pedsv_r_docker
  }

  call BurdenTests as ValidationCohortBurdenTests {
    input:
      beds = [validation_bed],
      bed_idxs = [validation_bed_idx],
      ad_matrixes = [validation_ad_matrix],
      ad_matrix_idxs = [validation_ad_matrix_idx],
      sample_metadata_tsv = sample_metadata_tsv,
      samples_list = validation_samples_list,
      af_fields = ["parent_AF"],
      ac_fields = ["parent_AC"],
      prefix = study_prefix + ".validation_cohort",
      docker = pedsv_r_docker
  }

  call UnifyOutputs {
    input:
      tarballs = [StudyWideSummaryPlots.plots_tarball,
                  StudyWideBurdenTests.plots_tarball,
                  TrioCohortSummaryPlots.plots_tarball,
                  ValidationCohortSummaryPlots.plots_tarball,
                  TrioCohortBurdenTests.plots_tarball,
                  ValidationCohortBurdenTests.plots_tarball],
      prefix = study_prefix + "_analysis_outputs",
      docker = ubuntu_docker,
      relocate_stats = true
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
    gzip -f StudyWideSummaryPlots/pca/*.tsv

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


# Conduct standard global burden tests
task BurdenTests {
  input {
    Array[File] beds
    Array[File] bed_idxs
    Array[File] ad_matrixes
    Array[File] ad_matrix_idxs
    File sample_metadata_tsv
    File samples_list
    Array[String] af_fields
    Array[String] ac_fields
    String prefix
    String docker
  }

  Int disk_gb = ceil(2 * (size(flatten([beds, ad_matrixes]), "GB"))) + 10

  command <<<
    set -eu -o pipefail

    # Make output directory
    mkdir ~{prefix}.BurdenTests

    # Run burden tests
    /opt/ped_germline_SV/analysis/landscape/global_burden_tests.R \
      ~{sep=" --bed " beds} \
      ~{sep=" --ad " ad_matrixes} \
      --metadata ~{sample_metadata_tsv} \
      --subset-samples ~{samples_list} \
      ~{sep=" --af-field " af_fields} \
      ~{sep=" --ac-field " ac_fields} \
      --out-prefix ~{prefix}.BurdenTests/~{prefix}
    gzip -f ~{prefix}.BurdenTests/*.tsv

    # Compress output
    tar -czvf ~{prefix}.BurdenTests.tar.gz ~{prefix}.BurdenTests
  >>>

  output {
    File plots_tarball = "~{prefix}.BurdenTests.tar.gz"
  }

  runtime {
    docker: docker
    memory: "31.5 GB"
    cpu: 16
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

  Int disk_gb = ceil(10 * size(tarballs, "GB")) + 20

  command <<<
    set -eu -o pipefail

    mkdir ~{prefix}_SVs_per_sample

    cat ~{write_lines(tarballs)} | xargs -I {} tar -xzvf {}

    while read ID; do
      find ./ -name "$ID.svids.list.gz" \
      | xargs -I {} zcat {} \
      | sort -V \
      | gzip -c \
      > ~{prefix}_SVs_per_sample/$ID.SV_IDs.list.gz
    done < ~{sample_list}

    tar -czvf ~{prefix}_SVs_per_sample.tar.gz ~{prefix}_SVs_per_sample
  >>>

  output {
    File per_sample_tarball = "~{prefix}_SVs_per_sample.tar.gz"
  }

  runtime {
    docker: docker
    memory: "3.5 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}


# Basic summary plots for a single cohort
task CohortSummaryPlots {
  input {
    File bed
    File bed_idx
    File sample_metadata_tsv
    File sample_list
    File per_sample_tarball
    String prefix
    
    String af_field = "AF"
    String ac_field = "AC"

    String docker
  }

  Int disk_gb = ceil(2 * size([bed], "GB")) + 10

  command <<<
    set -eu -o pipefail

    # Prep output directory
    mkdir ~{prefix}.SummaryPlots

    # Plot SV counts and sizes
    /opt/ped_germline_SV/analysis/landscape/plot_sv_counts_and_sizes.R \
      --af-field ~{af_field} \
      --ac-field ~{ac_field} \
      --out-prefix ~{prefix}.SummaryPlots/~{prefix} \
      ~{bed}

    # Format SV counts per sample
    mkdir perSample_data
    tar -xzvf ~{per_sample_tarball} -C perSample_data/
    while read sample; do
      find perSample_data/ -name "$sample.SV_IDs.list.gz" \
      | xargs -I {} zcat {} | wc -l \
      | paste <( echo $sample ) -
    done < ~{sample_list} \
    | cat <( echo -e "sample\tcount" ) - \
    > sv_counts_per_sample.tsv

    # Plot SVs per sample
    /opt/ped_germline_SV/analysis/landscape/plot_svs_per_sample.R \
      --metadata ~{sample_metadata_tsv} \
      --out-prefix ~{prefix}.SummaryPlots/~{prefix} \
      sv_counts_per_sample.tsv

    # Compress output
    find ~{prefix}.SummaryPlots/ -name "*.tsv" | xargs -I {} gzip -f {}
    tar -czvf ~{prefix}.SummaryPlots.tar.gz ~{prefix}.SummaryPlots
  >>>

  output {
    File plots_tarball = "~{prefix}.SummaryPlots.tar.gz"
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

    Boolean relocate_stats = false
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

    # Relocate all text files to stats subdirectory
    if [ ~{relocate_stats} == "true" ]; then
      mkdir ~{prefix}/stats
      find ~{prefix}/ -name "*.tsv.gz" | xargs -I {} mv {} ~{prefix}/stats/
    fi

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

