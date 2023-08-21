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


workflow SvRvas {
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

    String pedsv_docker
    String pedsv_r_docker
  }

  Array[Array[String]] contiglist = read_tsv(ref_fai)

  # Process all chromosomes in parallel
  scatter ( contig_info in contiglist ) {
    String contig = contig_info[0]

    call ContigRvas {
      input:
        gtf = gtf,
        sites_bed = sites_bed,
        sites_bed_idx = sites_bed_idx,
        ad_matrix = ad_matrix,
        ad_matrix_idx = ad_matrix_idx,
        sample_metadata_tsv = sample_metadata_tsv,
        samples_list = samples_list,
        contig = contig,
        prefix = prefix + "." + contig,
        docker = pedsv_r_docker
    }
  }

  # Concatenate sumstats across all chromosomes
  call ConcatSumstats {
    input:
      tsvs = ContigRvas.sumstats,
      prefix = prefix,
      docker = pedsv_docker
  }

  output {
    File rvas_sumstats = ConcatSumstats.combined_sumstats
  }
}


task ContigRvas {
  input {
    File gtf
    File sites_bed
    File sites_bed_idx
    File ad_matrix
    File ad_matrix_idx
    File sample_metadata_tsv
    File samples_list
    String contig
    String prefix
    String docker

    Float mem_gb = 15.5
    Int n_cpu = 8
  }

  Int disk_gb = ceil(1.5 * size([gtf, sites_bed, ad_matrix], "GB"))

  command <<<
    set -eu -o pipefail

    # Extract list of genes to analyze
    if [ $( file ~{gtf} | fgrep zip | wc -l ) -gt 0 ]; then
      zcat ~{gtf} 
    else
      cat ~{gtf}
    fi \
    | grep -e '^~{contig}' \
    | cut -f9 \
    | sed 's/\; /\n/g' \
    | fgrep -w gene_name \
    | sed 's/gene_name\ //g' \
    | tr -d '";' \
    | sort -V \
    | uniq \
    > eligible_genes.list

    # Run association test for all eligible genes
    if [ $( cat eligible_genes.list | wc -l ) -eq 0 ]; then
      touch ~{prefix}.sv_rvas_sumstats.tsv
    else
      /opt/ped_germline_SV/analysis/association/sv_genic_rvas.R \
        --bed ~{sites_bed} \
        --ad ~{ad_matrix} \
        --genes eligible_genes.list \
        --metadata ~{sample_metadata_tsv} \
        --subset-samples ~{samples_list} \
        --out-tsv ~{prefix}.sv_rvas_sumstats.tsv
      gzip -f ~{prefix}.sv_rvas_sumstats.tsv
    fi

  >>>

  output {
    File sumstats = "~{prefix}.sv_rvas_sumstats.tsv.gz"
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
    Array[File] tsvs
    String prefix
    String docker
  }

  Int disk_gb = ceil(2.5 * size(tsvs, "GB"))

  command <<<
    set -eu -o pipefail

    zcat ~{tsvs[0]} | head -n1 > ~{prefix}.sv_rvas.sumstats.tsv
    zcat ~{sep=" " tsvs} | grep -ve '^#' >> ~{prefix}.sv_rvas.sumstats.tsv
    gzip -f ~{prefix}.sv_rvas.sumstats.tsv
  >>>

  output {
    File combined_sumstats = "~{prefix}.sv_rvas.sumstats.tsv.gz"
  }

  runtime {
    docker: docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}

