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


import "https://raw.githubusercontent.com/vanallenlab/ped_germline_SV/main/wdl/SvGenicRvas.wdl" as Rvas


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
    File trio_dense_vcf
    File trio_dense_vcf_idx
		File trio_samples_list
    File? trio_variant_exclusion_list
    # String trio_gwas_bcftools_query_options = ""

    File case_control_bed
    File case_control_bed_idx
		File case_control_ad_matrix
		File case_control_ad_matrix_idx
    File case_control_dense_vcf
    File case_control_dense_vcf_idx
		File case_control_samples_list
    File? case_control_variant_exclusion_list
    # String case_control_gwas_bcftools_query_options = ""

    # Note: the below inputs are expected to contain unrelated analysis samples
    File full_cohort_bed
    File full_cohort_bed_idx
    File full_cohort_ad_matrix
    File full_cohort_ad_matrix_idx
    File full_cohort_dense_vcf
    File full_cohort_dense_vcf_idx
    File full_cohort_samples_list
    # String full_cohort_gwas_bcftools_query_options = ""

    # Note: the below inputs should contain ALL individuals in study, including relatives
    File full_cohort_w_relatives_vcf
    File full_cohort_w_relatives_vcf_idx
    File full_cohort_w_relatives_bed
    File full_cohort_w_relatives_bed_idx
    File full_cohort_w_relatives_ad_matrix
    File full_cohort_w_relatives_ad_matrix_idx
    File full_cohort_w_relatives_samples_list

    Array[File]? gene_lists
    Array[String]? gene_list_names
    File? genomic_disorder_bed
    Float genomic_disorder_recip_frac = 0.5
    File? locus_assoc_exclude_regions
    Float? rvas_exclusion_frac_overlap
    Float? sliding_window_exclusion_frac_overlap

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
      exclude_regions = locus_assoc_exclude_regions,
      exclusion_frac_overlap = rvas_exclusion_frac_overlap,
      pedsv_docker = pedsv_docker,
      pedsv_r_docker = pedsv_r_docker
  }

  # TODO: plot RVAS meta-analysis

  scatter ( contig_info in contiglist ) {
    String contig = contig_info[0]
    Int contig_len = contig_info[1]

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

    call PrecomputePerSampleBurdens {
      input:
        vcf = full_cohort_w_relatives_vcf,
        vcf_idx = full_cohort_w_relatives_vcf_idx,
        contig = contig,
        prefix = study_prefix + ".full_cohort_w_relatives",
        docker = pedsv_docker
    }

    call SlidingWindowTest as FullCohortSlidingWindowTest {
      input:
        sv_bed = full_cohort_bed,
        sv_bed_idx = full_cohort_bed_idx,
        ad_matrix = full_cohort_ad_matrix,
        ad_matrix_idx = full_cohort_ad_matrix_idx,
        contig = contig,
        contig_len = contig_len,
        exclude_regions = locus_assoc_exclude_regions,
        exclusion_frac_overlap = sliding_window_exclusion_frac_overlap,
        sample_metadata_tsv = sample_metadata_tsv,
        samples_list = all_samples_list,
        prefix = study_prefix + ".full_cohort." + contig,
        docker = pedsv_r_docker
    }

    call SlidingWindowTest as TrioCohortSlidingWindowTest {
      input:
        sv_bed = trio_bed,
        sv_bed_idx = trio_bed_idx,
        ad_matrix = trio_ad_matrix,
        ad_matrix_idx = trio_ad_matrix_idx,
        contig = contig,
        contig_len = contig_len,
        exclude_regions = locus_assoc_exclude_regions,
        exclusion_frac_overlap = sliding_window_exclusion_frac_overlap,
        sample_metadata_tsv = sample_metadata_tsv,
        samples_list = trio_samples_list,
        prefix = study_prefix + ".trio_cohort." + contig,
        docker = pedsv_r_docker
    }

    call SlidingWindowTest as CaseControlCohortSlidingWindowTest {
      input:
        sv_bed = case_control_bed,
        sv_bed_idx = case_control_bed_idx,
        ad_matrix = case_control_ad_matrix,
        ad_matrix_idx = case_control_ad_matrix_idx,
        contig = contig,
        contig_len = contig_len,
        exclude_regions = locus_assoc_exclude_regions,
        exclusion_frac_overlap = sliding_window_exclusion_frac_overlap,
        sample_metadata_tsv = sample_metadata_tsv,
        samples_list = case_control_samples_list,
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

  call CombinePerSampleBurdens {
    input:
      rare_sv_count_tsvs = PrecomputePerSampleBurdens.rare_sv_counts,
      vrare_sv_count_tsvs = PrecomputePerSampleBurdens.vrare_sv_counts,
      singleton_sv_count_tsvs = PrecomputePerSampleBurdens.singleton_sv_counts,
      rare_del_imbalance_tsvs = PrecomputePerSampleBurdens.rare_del_imbalance,
      vrare_del_imbalance_tsvs = PrecomputePerSampleBurdens.vrare_del_imbalance,
      singleton_del_imbalance_tsvs = PrecomputePerSampleBurdens.singleton_del_imbalance,
      rare_dup_imbalance_tsvs = PrecomputePerSampleBurdens.rare_dup_imbalance,
      vrare_dup_imbalance_tsvs = PrecomputePerSampleBurdens.vrare_dup_imbalance,
      singleton_dup_imbalance_tsvs = PrecomputePerSampleBurdens.singleton_dup_imbalance,
      prefix = study_prefix + ".full_cohort_w_relatives",
      docker = pedsv_docker
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

  call CohortSummaryPlots as AnalysisCohortSummaryPlots {
    input:
      bed = full_cohort_bed,
      bed_idx = full_cohort_bed_idx,
      prefix = study_prefix + ".full_cohort",
      sample_metadata_tsv = sample_metadata_tsv,
      sample_list = full_cohort_samples_list,
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
      beds = [full_cohort_bed],
      bed_idxs = [full_cohort_bed_idx],
      ad_matrixes = [full_cohort_ad_matrix],
      ad_matrix_idxs = [full_cohort_ad_matrix_idx],
      sample_metadata_tsv = sample_metadata_tsv,
      samples_list = all_samples_list,
      af_fields = ["POPMAX_AF"],
      ac_fields = ["AC"],
      variant_exclusion_list = sv_exclusion_list,
      gene_lists = gene_lists,
      gene_list_names = gene_list_names,
      genomic_disorder_bed = genomic_disorder_bed,
      genomic_disorder_recip_frac = genomic_disorder_recip_frac,
      precomputed_burden_stats = CombinePerSampleBurdens.precomputed_burden_stats_tsv,
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
      genomic_disorder_recip_frac = genomic_disorder_recip_frac,
      precomputed_burden_stats = CombinePerSampleBurdens.precomputed_burden_stats_tsv,
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
      genomic_disorder_recip_frac = genomic_disorder_recip_frac,
      precomputed_burden_stats = CombinePerSampleBurdens.precomputed_burden_stats_tsv,
      prefix = study_prefix + ".case_control_cohort",
      docker = pedsv_r_docker
  }

  call PlotSlidingWindowTests as PlotFullCohortSlidingWindows {
    input:
      stats = FullCohortSlidingWindowTest.bin_stats,
      prefix = study_prefix,
      docker = pedsv_r_docker
  }

  call PlotSlidingWindowTests as PlotTrioCohortSlidingWindows {
    input:
      stats = TrioCohortSlidingWindowTest.bin_stats,
      prefix = study_prefix + ".trio_cohort",
      docker = pedsv_r_docker
  }

  call PlotSlidingWindowTests as PlotCaseControlCohortSlidingWindows {
    input:
      stats = CaseControlCohortSlidingWindowTest.bin_stats,
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
      exclude_regions = locus_assoc_exclude_regions,
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
      exclude_regions = locus_assoc_exclude_regions,
      exclusion_frac_overlap = rvas_exclusion_frac_overlap,
      pedsv_docker = pedsv_docker,
      pedsv_r_docker = pedsv_r_docker
  }

  # TODO: plot case/control cohort RVAS

  call UnifyOutputs {
    input:
      tarballs = [FullCohortSummaryPlots.plots_tarball,
                  AnalysisCohortSummaryPlots.plots_tarball,
                  StudyWideBurdenTests.plots_tarball,
                  TrioCohortSummaryPlots.plots_tarball,
                  CaseControlCohortSummaryPlots.plots_tarball,
                  TrioCohortBurdenTests.plots_tarball,
                  CaseControlCohortBurdenTests.plots_tarball,
                  PlotFullCohortSlidingWindows.results_tarball,
                  PlotTrioCohortSlidingWindows.results_tarball,
                  PlotCaseControlCohortSlidingWindows.results_tarball],
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
    Float genomic_disorder_recip_frac = 0.5
    File? precomputed_burden_stats
    
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
    if [ ~{defined(precomputed_burden_stats)} == "true" ]; then
      cmd="$cmd --precomputed-burden-stats ~{precomputed_burden_stats}"
    fi
    cmd="$cmd --out-prefix ~{prefix}.BurdenTests/~{prefix}"

    # Prep list of SVs to evaluate for genomic disorder analyses, if optioned
    if [ ~{defined(genomic_disorder_bed)} == "true" ]; then
      while read bed; do
        bedtools intersect \
          -r -wa -wb -f ~{genomic_disorder_recip_frac} \
          -a $bed -b ~{select_first([genomic_disorder_bed])} \
        | awk -v OFS="\t" '{ if ($5==$NF) print $4, $(NF-1) }'
      done < ~{write_lines(beds)} \
      | sort -V | uniq > gd_hits.pairs.tsv
      cut -f1 gd_hits.pairs.tsv | sort -V | uniq > gd_hits.vids.list
      cmd="$cmd --genomic-disorder-hits gd_hits.vids.list"

      # If GD BED provided, also collect rates per disease per GD and run association tests per GD
      gd_cmd="/opt/ped_germline_SV/analysis/association/gd_assoc.R"
      gd_cmd="$gd_cmd $( awk -v OFS=" " -v ORS=" " '{ print "--bed", $1 }' ~{write_lines(beds)} )"
      gd_cmd="$gd_cmd $( awk -v OFS=" " -v ORS=" " '{ print "--ad", $1 }' ~{write_lines(ad_matrixes)} )"
      gd_cmd="$gd_cmd --metadata ~{sample_metadata_tsv} --subset-samples ~{samples_list}"
      gd_cmd="$gd_cmd --genomic-disorder-hits gd_hits.pairs.tsv --outfile ~{prefix}.BurdenTests/~{prefix}.gd_association_stats.tsv"
      echo -e "\n\nNow running genomic disorder association tests with the following command:\n$gd_cmd\n"
      eval $gd_cmd
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


# Run a sliding window CNV burden test for a single chromosome
task SlidingWindowTest {
  input {
    File sv_bed
    File sv_bed_idx
    File ad_matrix
    File ad_matrix_idx
    String contig
    Int contig_len
    File? exclude_regions
    Float exclusion_frac_overlap = 0.5
    File sample_metadata_tsv
    File samples_list
    Int bin_size = 1000000
    Int step_size = 250000
    Int min_sv_size = 50000
    String prefix
    String docker

    Float mem_gb = 3.75
    Int n_cpu = 2
    Int? disk_gb
  }

  Int default_disk_gb = ceil(3 * size([sv_bed, ad_matrix], "GB")) + 20

  command <<<
    set -eu -o pipefail

    # Make genomic bins
    paste \
      <( seq 0 ~{step_size} ~{contig_len} ) \
      <( seq ~{bin_size} ~{step_size} $(( ~{contig_len} + ~{bin_size} )) ) \
    | awk -v contig="~{contig}" -v OFS="\t" '{ print contig, $1, $2 }' \
    | awk -v bin_size=~{bin_size} -v last_end=$(( ~{contig_len} + ~{step_size} )) \
      '{ if ($3-$2 <= bin_size && $3 <= last_end) print }' \
    | bgzip -c \
    > ~{contig}.bins.pre_mask.bed.gz
    if [ ~{defined(exclude_regions)} == "true" ]; then
      bedtools coverage \
        -a ~{contig}.bins.pre_mask.bed.gz \
        -b ~{select_first(exclude_regions)} \
      | awk -v OFS="\t" -v frac=~{exclusion_frac_overlap} \
        '{ if ($NF < frac) print $1, $2, $3 }' \
      | bgzip -c \
      > ~{contig}.bins.bed.gz
    else
      cp ~{contig}.bins.pre_mask.bed.gz ~{contig}.bins.bed.gz
    fi
    tabix -p bed -f ~{contig}.bins.bed.gz

    # Filter SV data to contig of interest & CNVs greater than min_sv_size
    tabix -h ~{sv_bed} ~{contig} | bgzip -c > ~{contig}.all_svs.bed.gz
    tabix -p bed -f ~{contig}.all_svs.bed.gz
    /opt/ped_germline_SV/analysis/landscape/filter_cnvs_for_sliding_window_test.R \
      --bed-in ~{contig}.all_svs.bed.gz \
      --minimum-size ~{min_sv_size} \
      --bed-out ~{contig}.large_rare_cnvs.bed
    while read chrom start end vid cnv intervals; do
      if [ $cnv == "CPX" ]; then
        echo $intervals | sed -e 's/,/\n/g' -e 's/[_:-]/\t/g' \
        | awk -v min_size=~{min_sv_size} -v OFS="\t" -v vid=$vid \
          '{ if ($1 ~ /DEL|DUP/ && $3-$2 >= min_size) print $2, $3, $4, vid, $1 }'
      else
        echo -e "$chrom\t$start\t$end\t$vid\t$cnv"
      fi
    done < <( grep -ve '^#' ~{contig}.large_rare_cnvs.bed ) \
    | sort -Vk1,1 -k2,2n -k3,3n -k4,4V -k5,5V \
    | cat <( head -n1 ~{contig}.large_rare_cnvs.bed ) - \
    | bgzip -c \
    > ~{contig}.large_rare_cnvs.bed.gz
    tabix -p bed -f ~{contig}.large_rare_cnvs.bed.gz

    # Preprocess AD matrix to reduce memory requirement downstream
    tabix ~{ad_matrix} ~{contig} \
    | fgrep -wf <( zcat ~{contig}.large_rare_cnvs.bed.gz | cut -f4 | sed '1d' ) \
    | cat <( tabix -H ~{ad_matrix} ) - \
    | bgzip -c \
    > ~{contig}.large_rare_cnvs.ad.bed.gz
    tabix -p bed -f ~{contig}.large_rare_cnvs.ad.bed.gz

    # Find overlap between SVs and bins
    for CNV in DEL DUP; do
      zcat ~{contig}.large_rare_cnvs.bed.gz \
      | awk -v cnv=$CNV '{ if ($NF==cnv) print }' \
      > ~{contig}.large_rare_cnvs.$CNV.bed
      if [ $( cat ~{contig}.large_rare_cnvs.$CNV.bed | wc -l ) -gt 0 ]; then
        bedtools map \
          -a ~{contig}.bins.bed.gz \
          -b ~{contig}.large_rare_cnvs.$CNV.bed \
          -c 4 \
          -o distinct \
          -f $( echo -e "~{min_sv_size}\t~{bin_size}" \
                | awk '{ printf "%.6f\n", $1 / $2 }' ) \
        | cut -f4 \
        | awk '{ if ($1==".") print "NA"; else print $1 }' \
        > ~{contig}.hits.$CNV.txt
      else
        zcat ~{contig}.bins.bed.gz \
        | awk '{ print "NA" }' \
        > ~{contig}.hits.$CNV.txt
      fi
    done
    zcat ~{contig}.bins.bed.gz \
    | paste \
      - \
      ~{contig}.hits.DEL.txt \
      ~{contig}.hits.DUP.txt \
    | cat \
      <( echo -e "#chrom\tstart\tend\tDELs\tDUPs" ) - \
    | bgzip -c \
    > ~{contig}.counts.bed.gz

    # Compute test statistics for each bin
    /opt/ped_germline_SV/analysis/association/sliding_window_test.R \
      --bed ~{contig}.counts.bed.gz \
      --ad ~{contig}.large_rare_cnvs.ad.bed.gz \
      --metadata ~{sample_metadata_tsv} \
      --subset-samples ~{samples_list} \
      --outfile ~{prefix}.sliding_window_stats.bed
    bgzip -f ~{prefix}.sliding_window_stats.bed
  >>>

  output {
    File bin_stats = "~{prefix}.sliding_window_stats.bed.gz"
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


# Precompute per-sample genome-wide burden statistics that are too 
# memory intensive to do efficiently in R
task PrecomputePerSampleBurdens {
  input {
    File vcf
    File vcf_idx
    String contig
    String prefix
    String docker

    Float mem_gb = 7.5
    Int? disk_gb
  }

  Int default_disk_gb = ceil(3 * size([vcf], "GB")) + 20

  command <<<
    set -eu -o pipefail

    # Subset to chromosome of interest
    tabix -h ~{vcf} ~{contig} | bgzip -c > ~{contig}.vcf.gz

    # Get dummy list of samples in VCF with zero counts for filling missing samples
    bcftools query -l ~{contig}.vcf.gz | awk -v OFS="\t" '{ print $1, "0" }' \
    > zeroes.tsv

    # Filter VCF to isolate rare SVs
    bcftools view \
      --include 'INFO/gnomad_v3.1_sv_POPMAX_AF < 0.01 & AF < 0.01 & FILTER == "PASS"' \
      -Oz -o rare.vcf.gz \
      ~{contig}.vcf.gz
    tabix -p vcf -f rare.vcf.gz

    # Count rare SVs per sample
    bcftools query \
      --include 'GT="alt"' \
      --format '[%SAMPLE\t1\n]' \
      rare.vcf.gz \
    | cat - zeroes.tsv \
    | /opt/ped_germline_SV/analysis/utilities/sum_values_per_sample.py \
    | cat <( echo -e "#sample\trare_sv_count" ) - \
    | gzip -c \
    > ~{prefix}.rare_sv_per_sample.~{contig}.tsv.gz

    # Count very rare SVs per sample
    bcftools query \
      --include 'GT="alt" & AF < 0.001' \
      --format '[%SAMPLE\t1\n]' \
      rare.vcf.gz \
    | cat - zeroes.tsv \
    | /opt/ped_germline_SV/analysis/utilities/sum_values_per_sample.py \
    | cat <( echo -e "#sample\tvrare_sv_count" ) - \
    | gzip -c \
    > ~{prefix}.vrare_sv_per_sample.~{contig}.tsv.gz

    # Count singleton SVs per sample
    bcftools query \
      --include 'INFO/AC == 1 & GT="alt"' \
      --format '[%SAMPLE\t1\n]' \
      rare.vcf.gz \
    | cat - zeroes.tsv \
    | /opt/ped_germline_SV/analysis/utilities/sum_values_per_sample.py \
    | cat <( echo -e "#sample\tsingleton_sv_count" ) - \
    | gzip -c \
    > ~{prefix}.singleton_sv_per_sample.~{contig}.tsv.gz

    # Prepare for dosage imbalance-based calculations
    /opt/ped_germline_SV/analysis/landscape/prep_dosage_imbalance_vcf.py \
      --vcf-out rare.imbalance.vcf.gz \
      rare.vcf.gz

    # Collect CNV type-specific genomic imbalance per frequency tranche
    for CNV in DEL DUP; do
      case $CNV in
        DEL)
          info=DEL_IMBALANCE
          ;;
        DUP)
          info=DUP_IMBALANCE
          ;;
      esac

      # Sum rare imbalance per sample
      bcftools query \
        --include "GT=\"alt\" & INFO/$info > 0" \
        --format "[%SAMPLE\t%INFO/$info\n]" \
        rare.imbalance.vcf.gz \
      | cat - zeroes.tsv \
      | /opt/ped_germline_SV/analysis/utilities/sum_values_per_sample.py \
      | cat <( echo -e "#sample\trare_imbalance_$CNV" ) - \
      | gzip -c \
      > ~{prefix}.rare_imbalance_per_sample.$CNV.~{contig}.tsv.gz

      # Sum very rare imbalance per sample
      bcftools query \
        --include "GT=\"alt\" & INFO/AF < 0.001 & INFO/$info > 0" \
        --format "[%SAMPLE\t%INFO/$info\n]" \
        rare.imbalance.vcf.gz \
      | cat - zeroes.tsv \
      | /opt/ped_germline_SV/analysis/utilities/sum_values_per_sample.py \
      | cat <( echo -e "#sample\tvrare_imbalance_$CNV" ) - \
      | gzip -c \
      > ~{prefix}.vrare_imbalance_per_sample.$CNV.~{contig}.tsv.gz

      # Sum singleton imbalance per sample
      bcftools query \
        --include "GT=\"alt\" & INFO/AC == 1 & INFO/$info > 0" \
        --format "[%SAMPLE\t%INFO/$info\n]" \
        rare.imbalance.vcf.gz \
      | cat - zeroes.tsv \
      | /opt/ped_germline_SV/analysis/utilities/sum_values_per_sample.py \
      | cat <( echo -e "#sample\tsingleton_imbalance_$CNV" ) - \
      | gzip -c \
      > ~{prefix}.singleton_imbalance_per_sample.$CNV.~{contig}.tsv.gz

    done

  >>>

  output {
    File rare_sv_counts = "~{prefix}.rare_sv_per_sample.~{contig}.tsv.gz"
    File vrare_sv_counts = "~{prefix}.vrare_sv_per_sample.~{contig}.tsv.gz"
    File singleton_sv_counts = "~{prefix}.singleton_sv_per_sample.~{contig}.tsv.gz"
    File rare_del_imbalance = "~{prefix}.rare_imbalance_per_sample.DEL.~{contig}.tsv.gz"
    File vrare_del_imbalance = "~{prefix}.vrare_imbalance_per_sample.DEL.~{contig}.tsv.gz"
    File singleton_del_imbalance = "~{prefix}.singleton_imbalance_per_sample.DEL.~{contig}.tsv.gz"
    File rare_dup_imbalance = "~{prefix}.rare_imbalance_per_sample.DUP.~{contig}.tsv.gz"
    File vrare_dup_imbalance = "~{prefix}.vrare_imbalance_per_sample.DUP.~{contig}.tsv.gz"
    File singleton_dup_imbalance = "~{prefix}.singleton_imbalance_per_sample.DUP.~{contig}.tsv.gz"
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    preemptible: 3
  }
}


# Combine outputs generated by PrecomputePerSampleBurdens
task CombinePerSampleBurdens {
  input {
    Array[File] rare_sv_count_tsvs
    Array[File] vrare_sv_count_tsvs
    Array[File] singleton_sv_count_tsvs
    Array[File] rare_del_imbalance_tsvs
    Array[File] vrare_del_imbalance_tsvs
    Array[File] singleton_del_imbalance_tsvs
    Array[File] rare_dup_imbalance_tsvs
    Array[File] vrare_dup_imbalance_tsvs
    Array[File] singleton_dup_imbalance_tsvs
    String prefix
    String docker

    Float mem_gb = 3.5
    Int n_cpu = 2
    Int? disk_gb
  }

  Int default_disk_gb = ceil(2 * size(flatten([rare_sv_count_tsvs, singleton_sv_count_tsvs]), "GB")) + 20

  command <<<
    set -eu -o pipefail

    # Merge rare SV count data
    zcat ~{sep=" " rare_sv_count_tsvs} \
    | fgrep -v "#" \
    | /opt/ped_germline_SV/analysis/utilities/sum_values_per_sample.py \
    | cat <( echo -e "#sample\trare_sv_count" ) - \
    | gzip -c \
    > ~{prefix}.rare_sv_per_sample.tsv.gz

    # Merge vrare SV count data
    zcat ~{sep=" " vrare_sv_count_tsvs} \
    | fgrep -v "#" \
    | /opt/ped_germline_SV/analysis/utilities/sum_values_per_sample.py \
    | cat <( echo -e "#sample\tvrare_sv_count" ) - \
    | gzip -c \
    > ~{prefix}.vrare_sv_per_sample.tsv.gz

    # Merge singleton SV count data
    zcat ~{sep=" " singleton_sv_count_tsvs} \
    | fgrep -v "#" \
    | /opt/ped_germline_SV/analysis/utilities/sum_values_per_sample.py \
    | cat <( echo -e "#sample\tsingleton_sv_count" ) - \
    | gzip -c \
    > ~{prefix}.singleton_sv_per_sample.tsv.gz

    # Merge rare DEL imbalance data
    zcat ~{sep=" " rare_del_imbalance_tsvs} \
    | fgrep -v "#" \
    | /opt/ped_germline_SV/analysis/utilities/sum_values_per_sample.py \
    | cat <( echo -e "#sample\trare_imbalance_DEL" ) - \
    | gzip -c \
    > ~{prefix}.rare_del_imbalance_per_sample.tsv.gz

    # Merge vrare DEL imbalance data
    zcat ~{sep=" " vrare_del_imbalance_tsvs} \
    | fgrep -v "#" \
    | /opt/ped_germline_SV/analysis/utilities/sum_values_per_sample.py \
    | cat <( echo -e "#sample\tvrare_imbalance_DEL" ) - \
    | gzip -c \
    > ~{prefix}.vrare_del_imbalance_per_sample.tsv.gz

    # Merge singleton DEL imbalance data
    zcat ~{sep=" " singleton_del_imbalance_tsvs} \
    | fgrep -v "#" \
    | /opt/ped_germline_SV/analysis/utilities/sum_values_per_sample.py \
    | cat <( echo -e "#sample\tsingleton_imbalance_DEL" ) - \
    | gzip -c \
    > ~{prefix}.singleton_del_imbalance_per_sample.tsv.gz

    # Merge rare DUP imbalance data
    zcat ~{sep=" " rare_dup_imbalance_tsvs} \
    | fgrep -v "#" \
    | /opt/ped_germline_SV/analysis/utilities/sum_values_per_sample.py \
    | cat <( echo -e "#sample\trare_imbalance_DUP" ) - \
    | gzip -c \
    > ~{prefix}.rare_dup_imbalance_per_sample.tsv.gz

    # Merge vrare DUP imbalance data
    zcat ~{sep=" " vrare_dup_imbalance_tsvs} \
    | fgrep -v "#" \
    | /opt/ped_germline_SV/analysis/utilities/sum_values_per_sample.py \
    | cat <( echo -e "#sample\tvrare_imbalance_DUP" ) - \
    | gzip -c \
    > ~{prefix}.vrare_dup_imbalance_per_sample.tsv.gz

    # Merge singleton DUP imbalance data
    zcat ~{sep=" " singleton_dup_imbalance_tsvs} \
    | fgrep -v "#" \
    | /opt/ped_germline_SV/analysis/utilities/sum_values_per_sample.py \
    | cat <( echo -e "#sample\tsingleton_imbalance_DUP" ) - \
    | gzip -c \
    > ~{prefix}.singleton_dup_imbalance_per_sample.tsv.gz

    # Sum imbalance data across DEL and DUP within a frequency tranche
    for freq in rare vrare singleton; do
      zcat \
        ~{prefix}.${freq}_del_imbalance_per_sample.tsv.gz \
        ~{prefix}.${freq}_dup_imbalance_per_sample.tsv.gz \
      | fgrep -v "#" \
      | /opt/ped_germline_SV/analysis/utilities/sum_values_per_sample.py \
      | cat <( echo -e "#sample\t${freq}_imbalance_CNV" ) - \
      | gzip -c \
      > ~{prefix}.${freq}_cnv_imbalance_per_sample.tsv.gz
    done

    # Column-wise merge of all files above
    /opt/ped_germline_SV/analysis/utilities/join_on_first_column.py \
      --join outer \
      --tsv-out ~{prefix}.precomputed_burden_stats.tsv.gz \
      ~{prefix}.rare_sv_per_sample.tsv.gz \
      ~{prefix}.vrare_sv_per_sample.tsv.gz \
      ~{prefix}.singleton_sv_per_sample.tsv.gz \
      ~{prefix}.rare_del_imbalance_per_sample.tsv.gz \
      ~{prefix}.vrare_del_imbalance_per_sample.tsv.gz \
      ~{prefix}.singleton_del_imbalance_per_sample.tsv.gz \
      ~{prefix}.rare_dup_imbalance_per_sample.tsv.gz \
      ~{prefix}.vrare_dup_imbalance_per_sample.tsv.gz \
      ~{prefix}.singleton_dup_imbalance_per_sample.tsv.gz \
      ~{prefix}.rare_cnv_imbalance_per_sample.tsv.gz \
      ~{prefix}.vrare_cnv_imbalance_per_sample.tsv.gz \
      ~{prefix}.singleton_cnv_imbalance_per_sample.tsv.gz
  >>>

  output {
    File precomputed_burden_stats_tsv = "~{prefix}.precomputed_burden_stats.tsv.gz"
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
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

  Int disk_gb = ceil(2 * size([bed, sample_metadata_tsv], "GB")) + 10

  command <<<
    set -eu -o pipefail

    # Prep output directory
    mkdir ~{prefix}.SummaryPlots

    # Plot PCs colored by ancestry
    /opt/ped_germline_SV/analysis/cohort_summaries/plot_pcs.R \
      --subset-samples ~{sample_list} \
      --out-prefix ~{prefix}.SummaryPlots/~{prefix} \
      ~{sample_metadata_tsv}

    # Plot sex ploidy colored by sex assignment
    /opt/ped_germline_SV/analysis/cohort_summaries/plot_sex_ploidy.R \
      --subset-samples ~{sample_list} \
      --out-prefix ~{prefix}.SummaryPlots/~{prefix} \
      ~{sample_metadata_tsv}

    # Plot SV counts and sizes
    /opt/ped_germline_SV/analysis/landscape/plot_sv_site_summary.R \
      --af-field ~{af_field} \
      --ac-field ~{ac_field} \
      --out-prefix ~{prefix}.SummaryPlots/~{prefix} \
      ~{bed}

    # Plot AF correlations vs. gnomAD v4 SVs
    /opt/ped_germline_SV/analysis/landscape/gnomad_af_comparison.R \
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


# Merge & plot sliding window test statistics
task PlotSlidingWindowTests {
  input {
    Array[File] stats
    String prefix
    String docker
  }

  Int disk_gb = ceil(2 * size(stats, "GB")) + 10

  command <<<
    set -eu -o pipefail

    mkdir ~{prefix}.SlidingWindowResults

    # Merge stats
    zcat ~{sep=" " stats} | grep -ve '^#' \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | cat <( zcat ~{stats[0]} | head -n1 ) - \
    | bgzip -c \
    > ~{prefix}.sliding_window_stats.bed.gz
    tabix -p bed -f ~{prefix}.sliding_window_stats.bed.gz

    # Visualize
    # TODO: implement this

    # Compress output
    tar -czvf ~{prefix}.SlidingWindowResults.tar.gz ~{prefix}.SlidingWindowResults
  >>>

  output {
    File results_tarball = "~{prefix}.SlidingWindowResults.tar.gz"
  }

  runtime {
    docker: docker
    memory: "3.5 GB"
    cpu: 2
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
      find ~{prefix}/ -name "*.bed.gz" | xargs -I {} mv {} ~{prefix}/stats/
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

