#########################
#    Germline SVs in    #
#   Pediatric Cancers   #
#########################
#
# PedSvMainAnalysis.wdl
#
# Main analysis WDL for global/landscaping arm of pediatric cancer germline WGS study
#
# Copyright (c) 2023-Present, the Van Allen laboratory at Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0


version 1.0


import "https://raw.githubusercontent.com/vanallenlab/ped_germline_SV/rlc_landscape_v2.2/wdl/SvGenicRvas.wdl" as Rvas


workflow PedSvMainAnalysis {
	input {
		File sample_metadata_tsv
    File ref_fai
    File gtf
    String study_prefix

    File trio_bed
    File trio_bed_idx
		File trio_ad_matrix
		File trio_ad_matrix_idx
		File trio_samples_list
    File? trio_variant_exclusion_list

    File case_control_bed
    File case_control_bed_idx
		File case_control_ad_matrix
		File case_control_ad_matrix_idx
		File case_control_samples_list
    File? case_control_variant_exclusion_list

    # Note: the below inputs are expected to contain unrelated analysis samples
    File full_cohort_bed
    File full_cohort_bed_idx
    File full_cohort_ad_matrix
    File full_cohort_ad_matrix_idx
    File full_cohort_samples_list

    # Note: the below inputs should contain ALL individuals in study, including relatives
    File full_cohort_w_relatives_bed
    File full_cohort_w_relatives_bed_idx
    File full_cohort_w_relatives_ad_matrix
    File full_cohort_w_relatives_ad_matrix_idx
    File full_cohort_w_relatives_samples_list

    Array[File]? gene_lists
    Array[String]? gene_list_names
    File? genomic_disorder_bed
    File? rvas_exclude_regions
    Float? rvas_exclusion_frac_overlap

    String pedsv_docker
    String pedsv_r_docker
    String ubuntu_docker = "ubuntu:latest"
	}

  Array[Array[String]] contiglist = read_tsv(ref_fai)

  call ConcatTextFiles as ConcatSampleLists {
    input:
      infiles = [trio_samples_list, case_control_samples_list],
      outfile_name = study_prefix + ".both_cohorts.samples.list",
      docker = ubuntu_docker
  }
  File all_samples_list = ConcatSampleLists.outfile

  if (defined(trio_variant_exclusion_list) || defined(case_control_variant_exclusion_list)) {
    call ConcatTextFiles as ConcatVariantExclusionLists {
      input:
        infiles = select_all([trio_variant_exclusion_list, 
                              case_control_variant_exclusion_list]),
        outfile_name = study_prefix + ".both_cohorts.variant_exclusion.list",
        docker = ubuntu_docker
    }
  }
  File? sv_exclusion_list = ConcatVariantExclusionLists.outfile

  call StudyWideSummaryPlots {
    input:
      sample_metadata_tsv = sample_metadata_tsv,
      samples_list = full_cohort_samples_list,
      analysis_samples_list = all_samples_list,
      prefix = study_prefix,
      docker = pedsv_r_docker
  }

  call Rvas.SvGenicRvas as StudyWideRvas {
    input:
      gtf = gtf,
      sites_bed = full_cohort_bed,
      sites_bed_idx = full_cohort_bed_idx,
      ad_matrix = full_cohort_ad_matrix,
      ad_matrix_idx = full_cohort_ad_matrix_idx,
      sample_metadata_tsv = sample_metadata_tsv,
      samples_list = all_samples_list,
      ref_fai = ref_fai,
      prefix = study_prefix,
      exclude_regions = rvas_exclude_regions,
      exclusion_frac_overlap = rvas_exclusion_frac_overlap,
      pedsv_docker = pedsv_docker,
      pedsv_r_docker = pedsv_r_docker
  }

  # TODO: plot RVAS meta-analysis

  scatter ( contig_info in contiglist ) {
    String contig = contig_info[0]

    call GetSVsPerSample as GetFullCohortSVsPerSample {
      input:
        sv_bed = full_cohort_w_relatives_bed,
        sv_bed_idx = full_cohort_w_relatives_bed_idx,
        ad_matrix = full_cohort_w_relatives_ad_matrix,
        ad_matrix_idx = full_cohort_w_relatives_ad_matrix_idx,
        contig = contig,
        prefix = study_prefix + ".full_cohort_w_relatives." + contig,
        docker = pedsv_r_docker,
        mem_gb = 31,
        n_cpu = 8
    }

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

    call GetSVsPerSample as GetCaseControlSVsPerSample {
      input:
        sv_bed = case_control_bed,
        sv_bed_idx = case_control_bed_idx,
        ad_matrix = case_control_ad_matrix,
        ad_matrix_idx = case_control_ad_matrix_idx,
        contig = contig,
        prefix = study_prefix + ".case_control_cohort." + contig,
        docker = pedsv_r_docker
    }
  }

  call MergeSVsPerSample as MergeFullCohortSVsPerSample {
    input:
      tarballs = GetFullCohortSVsPerSample.per_sample_tarball,
      sample_list = full_cohort_w_relatives_samples_list,
      prefix = study_prefix + ".full_cohort_w_relatives",
      docker = ubuntu_docker
  }

  call MergeSVsPerSample as MergeTrioSVsPerSample {
    input:
      tarballs = GetTrioSVsPerSample.per_sample_tarball,
      sample_list = trio_samples_list,
      prefix = study_prefix + ".trio_cohort",
      docker = ubuntu_docker
  }

  call MergeSVsPerSample as MergeCaseControlSVsPerSample {
    input:
      tarballs = GetCaseControlSVsPerSample.per_sample_tarball,
      sample_list = case_control_samples_list,
      prefix = study_prefix + ".case_control_cohort",
      docker = ubuntu_docker
  }

  call CohortSummaryPlots as FullCohortSummaryPlots {
    input:
      bed = full_cohort_w_relatives_bed,
      bed_idx = full_cohort_w_relatives_bed_idx,
      prefix = study_prefix + ".full_cohort_w_relatives",
      sample_metadata_tsv = sample_metadata_tsv,
      sample_list = full_cohort_w_relatives_samples_list,
      per_sample_tarball = MergeFullCohortSVsPerSample.per_sample_tarball,
      af_field = "POPMAX_AF",
      ac_field = "AC",
      docker = pedsv_r_docker
  }

  call CohortSummaryPlots as TrioCohortSummaryPlots {
    input:
      bed = trio_bed,
      bed_idx = trio_bed_idx,
      prefix = study_prefix + ".trio_cohort",
      sample_metadata_tsv = sample_metadata_tsv,
      sample_list = trio_samples_list,
      per_sample_tarball = MergeTrioSVsPerSample.per_sample_tarball,
      af_field = "POPMAX_AF",
      ac_field = "AC",
      docker = pedsv_r_docker
  }

  call CohortSummaryPlots as CaseControlCohortSummaryPlots {
    input:
      bed = case_control_bed,
      bed_idx = case_control_bed_idx,
      sample_metadata_tsv = sample_metadata_tsv,
      sample_list = case_control_samples_list,
      per_sample_tarball = MergeCaseControlSVsPerSample.per_sample_tarball,
      prefix = study_prefix + ".case_control_cohort",
      docker = pedsv_r_docker
  }

  call BurdenTests as StudyWideBurdenTests {
    input:
      beds = [trio_bed, case_control_bed],
      bed_idxs = [trio_bed_idx, case_control_bed_idx],
      ad_matrixes = [trio_ad_matrix, case_control_ad_matrix],
      ad_matrix_idxs = [trio_ad_matrix_idx, case_control_ad_matrix_idx],
      sample_metadata_tsv = sample_metadata_tsv,
      samples_list = all_samples_list,
      af_fields = ["POPMAX_AF", "POPMAX_AF"],
      ac_fields = ["AC", "AC"],
      variant_exclusion_list = sv_exclusion_list,
      gene_lists = gene_lists,
      gene_list_names = gene_list_names,
      genomic_disorder_bed = genomic_disorder_bed,
      prefix = study_prefix,
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
      af_fields = ["POPMAX_AF"],
      ac_fields = ["AC"],
      variant_exclusion_list = sv_exclusion_list,
      gene_lists = gene_lists,
      gene_list_names = gene_list_names,
      genomic_disorder_bed = genomic_disorder_bed,
      prefix = study_prefix + ".trio_cohort",
      docker = pedsv_r_docker
  }

  call BurdenTests as CaseControlCohortBurdenTests {
    input:
      beds = [case_control_bed],
      bed_idxs = [case_control_bed_idx],
      ad_matrixes = [case_control_ad_matrix],
      ad_matrix_idxs = [case_control_ad_matrix_idx],
      sample_metadata_tsv = sample_metadata_tsv,
      samples_list = case_control_samples_list,
      af_fields = ["POPMAX_AF"],
      ac_fields = ["AC"],
      variant_exclusion_list = sv_exclusion_list,
      gene_lists = gene_lists,
      gene_list_names = gene_list_names,
      genomic_disorder_bed = genomic_disorder_bed,
      prefix = study_prefix + ".case_control_cohort",
      docker = pedsv_r_docker
  }

  call Rvas.SvGenicRvas as TrioCohortRvas {
    input:
      gtf = gtf,
      sites_bed = trio_bed,
      sites_bed_idx = trio_bed_idx,
      ad_matrix = trio_ad_matrix,
      ad_matrix_idx = trio_ad_matrix_idx,
      sample_metadata_tsv = sample_metadata_tsv,
      samples_list = trio_samples_list,
      ref_fai = ref_fai,
      prefix = study_prefix + ".trio_cohort",
      exclude_regions = rvas_exclude_regions,
      exclusion_frac_overlap = rvas_exclusion_frac_overlap,
      pedsv_docker = pedsv_docker,
      pedsv_r_docker = pedsv_r_docker
  }

  # TODO: plot trio RVAS

  call Rvas.SvGenicRvas as CaseControlCohortRvas {
    input:
      gtf = gtf,
      sites_bed = case_control_bed,
      sites_bed_idx = case_control_bed_idx,
      ad_matrix = case_control_ad_matrix,
      ad_matrix_idx = case_control_ad_matrix_idx,
      sample_metadata_tsv = sample_metadata_tsv,
      samples_list = case_control_samples_list,
      ref_fai = ref_fai,
      prefix = study_prefix + ".case_control_cohort",
      exclude_regions = rvas_exclude_regions,
      exclusion_frac_overlap = rvas_exclusion_frac_overlap,
      pedsv_docker = pedsv_docker,
      pedsv_r_docker = pedsv_r_docker
  }

  # TODO: plot case/control cohort RVAS

  call UnifyOutputs {
    input:
      tarballs = [StudyWideSummaryPlots.plots_tarball,
                  FullCohortSummaryPlots.plots_tarball,
                  StudyWideBurdenTests.plots_tarball,
                  TrioCohortSummaryPlots.plots_tarball,
                  CaseControlCohortSummaryPlots.plots_tarball,
                  TrioCohortBurdenTests.plots_tarball,
                  CaseControlCohortBurdenTests.plots_tarball],
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
    File? variant_exclusion_list
    Array[File]? gene_lists
    Array[String]? gene_list_names
    File? genomic_disorder_bed
    String prefix
    String docker
  }

  Int disk_gb = ceil(2 * (size(flatten([beds, ad_matrixes]), "GB"))) + 10
  Boolean any_gene_lists = defined(gene_lists)

  command <<<
    set -eu -o pipefail

    # Make output directory
    mkdir ~{prefix}.BurdenTests

    # Prep gene list input, if optioned
    if [ ~{any_gene_lists} == "true" ]; then
      paste \
        ~{write_lines(select_first([gene_list_names]))} \
        ~{write_lines(select_first([gene_lists]))} \
      > gene_list.input.tsv
    fi

    # Build burden test command
    cmd="/opt/ped_germline_SV/analysis/landscape/global_burden_tests.R"
    cmd="$cmd $( awk -v OFS=" " -v ORS=" " '{ print "--bed", $1 }' ~{write_lines(beds)} )"
    cmd="$cmd $( awk -v OFS=" " -v ORS=" " '{ print "--ad", $1 }' ~{write_lines(ad_matrixes)} )"
    cmd="$cmd --metadata ~{sample_metadata_tsv} --subset-samples ~{samples_list}"
    cmd="$cmd $( awk -v OFS=" " -v ORS=" " '{ print "--af-field", $1 }' ~{write_lines(af_fields)} )"
    cmd="$cmd $( awk -v OFS=" " -v ORS=" " '{ print "--ac-field", $1 }' ~{write_lines(ac_fields)} )"
    if [ ~{any_gene_lists} == "true" ]; then
      cmd="$cmd --gene-lists gene_list.input.tsv"
    fi
    if [ ~{defined(variant_exclusion_list)} == "true" ]; then
      cmd="$cmd --exclude-variants ~{variant_exclusion_list}"
    fi
    cmd="$cmd --out-prefix ~{prefix}.BurdenTests/~{prefix}"

    # Prep list of SVs to evaluate for genomic disorder analyses, if optioned
    if [ ~{defined(genomic_disorder_bed)} == "true" ]; then
      while read bed; do
        bedtools intersect \
          -r -wa -wb -f 0.8 \
          -a $bed -b ~{select_first([genomic_disorder_bed])} \
        | awk '{ if ($5==$NF) print $4 }'
      done < ~{write_lines(beds)} \
      | sort -V | uniq > gd_hits.vids.list
      cmd="$cmd --genomic-disorder-hits gd_hits.vids.list"
    fi

    # Run burden tests
    echo -e "\n\nNow running burden tests with the following command:\n$cmd\n"
    eval $cmd
    gzip -f ~{prefix}.BurdenTests/*.tsv

    # Generate summary plots based on the stats from burden test
    /opt/ped_germline_SV/analysis/landscape/plot_global_burden_sumstats.R \
      "~{prefix}.BurdenTests/~{prefix}.global_burden_tests.tsv.gz" \
      "~{prefix}.BurdenTests/~{prefix}"

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

    Float mem_gb = 15.5
    Int n_cpu = 4
    Int? disk_gb
  }

  Int default_disk_gb = ceil(3 * size([sv_bed, ad_matrix], "GB")) + 20

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
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
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
    File? per_sample_tarball
    String prefix
    
    String af_field = "POPMAX_AF"
    String ac_field = "AC"

    String docker
  }

  Int disk_gb = ceil(2 * size([bed], "GB")) + 10

  command <<<
    set -eu -o pipefail

    # Prep output directory
    mkdir ~{prefix}.SummaryPlots

    # Plot SV counts and sizes
    /opt/ped_germline_SV/analysis/landscape/plot_sv_site_summary.R \
      --af-field ~{af_field} \
      --ac-field ~{ac_field} \
      --out-prefix ~{prefix}.SummaryPlots/~{prefix} \
      ~{bed}

    # Format autosomal, biallelic SV counts per sample
    mkdir perSample_data
    if [ ~{defined(per_sample_tarball)} == "true" ]; then
      tar -xzvf ~{per_sample_tarball} -C perSample_data/
      while read sample; do
        find perSample_data/ -name "$sample.SV_IDs.list.gz" \
        | xargs -I {} zcat {} \
        | fgrep -v "_CNV_" \
        | fgrep -v "_chrX" \
        | fgrep -v "_chrY" \
        | wc -l \
        | paste <( echo $sample ) -
      done < ~{sample_list} \
      | cat <( echo -e "sample\tcount" ) - \
      > ~{prefix}.SummaryPlots/sv_counts_per_sample.tsv

      # Plot SVs per sample
      /opt/ped_germline_SV/analysis/landscape/plot_svs_per_sample.R \
        --metadata ~{sample_metadata_tsv} \
        --out-prefix ~{prefix}.SummaryPlots/~{prefix} \
        ~{prefix}.SummaryPlots/sv_counts_per_sample.tsv
     fi

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

