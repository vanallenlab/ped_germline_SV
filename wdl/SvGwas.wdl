#########################
#    Germline SVs in    #
#   Pediatric Cancers   #
#########################
#
# SvGwas.wdl
#
# Perform genome-wide case:control association testing for individual SVs
#
# Copyright (c) 2023-Present, the Van Allen laboratory at Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0


version 1.0


workflow SvGwas {
  input {
    File dense_vcf
    File dense_vcf_idx
    File ad_matrix
    File ad_matrix_idx
    File sample_metadata_tsv
    File samples_list
    File ref_fai
    String prefix

    File? custom_gwas_script
    String bcftools_query_options = ""

    Float gwas_mem_gb = 15.5

    String pedsv_docker
    String pedsv_r_docker
  }

  Array[Array[String]] contiglist = read_tsv(ref_fai)

  # Process all chromosomes in parallel
  scatter ( contig_info in contiglist ) {
    String contig = contig_info[0]

    call PrepVariants {
      input:
        dense_vcf = dense_vcf,
        dense_vcf_idx = dense_vcf_idx,
        contig = contig,
        samples_list = samples_list,
        bcftools_query_options = bcftools_query_options,
        prefix = prefix + "." + contig,
        docker = pedsv_docker
    }

    call ContigGwas {
      input:
        ad_matrix = ad_matrix,
        ad_matrix_idx = ad_matrix_idx,
        sample_metadata_tsv = sample_metadata_tsv,
        samples_list = samples_list,
        variant_ids_to_test = PrepVariants.vids_list,
        prefix = prefix + "." + contig,
        custom_gwas_script = custom_gwas_script,
        mem_gb = gwas_mem_gb,
        docker = pedsv_r_docker
    }
  }

  # Concatenate sumstats across all chromosomes
  call ConcatSumstats {
    input:
      tsvs = ContigGwas.sumstats,
      prefix = prefix,
      docker = pedsv_docker
  }

  output {
    File gwas_sumstats = ConcatSumstats.combined_sumstats
  }
}


task PrepVariants {
  input {
    File dense_vcf
    File dense_vcf_idx
    String contig
    File samples_list
    String bcftools_query_options
    String prefix
    String docker
  }

  Int disk_gb = ceil(2 * size(dense_vcf, "GB"))

  command <<<
    set -eu -o pipefail

    # Fill missing annotations, filter, and extract variant IDs to test
    bcftools view \
      --regions "~{contig}" \
      --samples-file ~{samples_list} \
      --force-samples \
      ~{dense_vcf} \
    | bcftools +fill-tags - \
      -- -t HWE,F_MISSING \
    | bcftools query - \
      -f "%ID\n" \
      ~{bcftools_query_options} \
    > "~{prefix}.vids.list"
  >>>

  output {
    File vids_list = "~{prefix}.vids.list"
  }

  runtime {
    docker: docker
    memory: "3.75 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}


task ContigGwas {
  input {
    File ad_matrix
    File ad_matrix_idx
    File sample_metadata_tsv
    File samples_list
    File variant_ids_to_test
    String prefix
    String docker

    File? custom_gwas_script

    Float mem_gb
  }

  Int disk_gb = ceil(2 * size(ad_matrix, "GB"))
  String gwas_script = select_first([custom_gwas_script, "/opt/ped_germline_SV/analysis/association/sv_gwas.R"])

  command <<<
    set -eu -o pipefail

    # Run association test for all qualifying variants
    ~{gwas_script} \
      --ad ~{ad_matrix} \
      --vids ~{variant_ids_to_test} \
      --metadata ~{sample_metadata_tsv} \
      --subset-samples ~{samples_list} \
      --out-tsv ~{prefix}.sv_gwas_sumstats.tsv
    gzip -f ~{prefix}.sv_gwas_sumstats.tsv

  >>>

  output {
    File sumstats = "~{prefix}.sv_gwas_sumstats.tsv.gz"
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}


task ConcatSumstats {
  input {
    Array[File] tsvs
    String prefix
    String docker
  }

  Int disk_gb = ceil(2.5 * size(tsvs, "GB"))

  command <<<
    set -eu -o pipefail

    zcat ~{tsvs[0]} | head -n1 > ~{prefix}.sv_gwas.sumstats.tsv
    zcat ~{sep=" " tsvs} | grep -ve '^#' >> ~{prefix}.sv_gwas.sumstats.tsv
    gzip -f ~{prefix}.sv_gwas.sumstats.tsv
  >>>

  output {
    File combined_sumstats = "~{prefix}.sv_gwas.sumstats.tsv.gz"
  }

  runtime {
    docker: docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}

