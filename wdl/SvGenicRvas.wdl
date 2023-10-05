#########################
#    Germline SVs in    #
#   Pediatric Cancers   #
#########################
#
# SvRvas.wdl
#
# Perform exome-wide gene-based case:control association testing for rare SVs
#
# Copyright (c) 2023-Present, the Van Allen laboratory at Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0


version 1.0


workflow SvGenicRvas {
  input {
    File gtf
    File sites_bed
    File sites_bed_idx
    File ad_matrix
    File ad_matrix_idx
    File sample_metadata_tsv
    File samples_list
    File ref_fai
    String prefix

    File? custom_rvas_script
    File? exclude_regions
    Float? exclusion_frac_overlap

    String pedsv_docker
    String pedsv_r_docker
  }

  Array[Array[String]] contiglist = read_tsv(ref_fai)

  # Process all chromosomes in parallel
  scatter ( contig_info in contiglist ) {
    String contig = contig_info[0]

    call PreprocessGtf {
      input:
        gtf = gtf,
        contig = contig,
        exclude_regions = exclude_regions,
        exclusion_frac_overlap = exclusion_frac_overlap,
        docker = pedsv_docker
    }

    call ContigRvas {
      input:
        sites_bed = sites_bed,
        sites_bed_idx = sites_bed_idx,
        ad_matrix = ad_matrix,
        ad_matrix_idx = ad_matrix_idx,
        sample_metadata_tsv = sample_metadata_tsv,
        samples_list = samples_list,
        eligible_genes_bed = PreprocessGtf.eligible_genes_bed,
        prefix = prefix + "." + contig,
        custom_rvas_script = custom_rvas_script,
        docker = pedsv_r_docker
    }
  }

  # Concatenate sumstats across all chromosomes
  call ConcatSumstats {
    input:
      beds = ContigRvas.sumstats,
      prefix = prefix,
      docker = pedsv_docker
  }

  output {
    File rvas_sumstats = ConcatSumstats.combined_sumstats
  }
}


task PreprocessGtf {
  input {
    File gtf
    String contig
    File? exclude_regions
    Float exclusion_frac_overlap = 0.5
    String docker
  }

  command <<<
    set -eu -o pipefail

    if [ ~{defined(exclude_regions)} == "true" ]; then
      exclude_option="--exclusion-bed ~{exclude_regions} --exclusion-frac ~{exclusion_frac_overlap}"
    fi

    /opt/ped_germline_SV/analysis/utilities/preprocess_gtf_for_rvas.py \
      ~{gtf} \
      --chromosome ~{contig} \
      $exclude_option \
      --outfile "eligible_genes.~{contig}.bed"
    bgzip -f "eligible_genes.~{contig}.bed"
    tabix -f "eligible_genes.~{contig}.bed.gz"
  >>>

  output {
    File eligible_genes_bed = "eligible_genes.~{contig}.bed.gz"
  }

  runtime {
    docker: docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk 20 HDD"
    preemptible: 3
  }
}


task ContigRvas {
  input {
    File sites_bed
    File sites_bed_idx
    File ad_matrix
    File ad_matrix_idx
    File sample_metadata_tsv
    File samples_list
    File eligible_genes_bed
    String prefix
    String docker

    File? custom_rvas_script

    Float mem_gb = 15.5
    Int n_cpu = 8
  }

  Int disk_gb = ceil(1.5 * size([sites_bed, ad_matrix], "GB"))
  String rvas_script = if defined(custom_rvas_script) then basename(select_first([custom_rvas_script])) else "/opt/ped_germline_SV/analysis/association/sv_genic_rvas.R"

  command <<<
    set -eu -o pipefail

    # Relocate custom script to execution directory if necessary
    if [ ~{defined(custom_rvas_script)} == "true" ]; then
      mv ~{select_first([custom_rvas_script])} ./
      chmod a+x ~{rvas_script}
    fi

    # Run association test for all eligible genes
    if [ $( cat ~{eligible_genes_bed} | wc -l ) -eq 0 ]; then
      touch ~{prefix}.sv_rvas_sumstats.bed
    else
      ~{rvas_script} \
        --bed ~{sites_bed} \
        --ad ~{ad_matrix} \
        --genes ~{eligible_genes_bed} \
        --metadata ~{sample_metadata_tsv} \
        --subset-samples ~{samples_list} \
        --out-bed ~{prefix}.sv_rvas_sumstats.bed
    fi
    bgzip -f ~{prefix}.sv_rvas_sumstats.bed
  >>>

  output {
    File sumstats = "~{prefix}.sv_rvas_sumstats.bed.gz"
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}


task ConcatSumstats {
  input {
    Array[File] beds
    String prefix
    String docker
  }

  Int disk_gb = ceil(10 * size(beds, "GB")) + 25

  command <<<
    set -eu -o pipefail

    zcat ~{sep=" " beds} \
    | grep -ve '^#' | grep -vP '\tstart\tend\tgene\t' \
    | sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
    | cat <( zcat ~{beds[0]} | head -n1 ) - \
    | bgzip -c \
    > ~{prefix}.sv_rvas.sumstats.bed.gz
    tabix -f ~{prefix}.sv_rvas.sumstats.bed.gz
    ls -ltrh
  >>>

  output {
    File combined_sumstats = "~{prefix}.sv_rvas.sumstats.bed.gz"
  }

  runtime {
    docker: docker
    memory: "3.75 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    max_retries: 3
  }
}

