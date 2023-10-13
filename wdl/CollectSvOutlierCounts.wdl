#########################
#    Germline SVs in    #
#   Pediatric Cancers   #
#########################
#
# CollectSvOutlierCounts.wdl
#
# Helper WDL to collect counts of various subsets of SVs for defining outliers
#
# Copyright (c) 2023-Present, the Van Allen laboratory at Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0


version 1.0


import "https://raw.githubusercontent.com/vanallenlab/pancan_germline_wgs/main/wdl/Utilities.wdl" as Utils


workflow CollectSvOutlierCounts {
  input {
    File vcf
    File vcf_idx
    File tss_bed
    File exons_bed
    File autosomes_fai
    String prefix

    String sv_pipeline_docker
  }

  call Utils.GetContigsFromFai as GetContigs {
    input:
      ref_fai = autosomes_fai,
      docker = sv_pipeline_docker
  }

  call CollectCounts as CountAll {
    input:
      vcf = vcf,
      vcf_idx = vcf_idx,
      contigs = GetContigs.contigs,
      bcftools_options = "--include \"FILTER == 'PASS'\"",
      output_prefix = "all",
      docker = sv_pipeline_docker
  }

  call CollectCounts as CountRare {
    input:
      vcf = vcf,
      vcf_idx = vcf_idx,
      contigs = GetContigs.contigs,
      bcftools_options = "--include \"FILTER == 'PASS' & INFO/AF < 0.01\"",
      output_prefix = "rare",
      docker = sv_pipeline_docker
  }

  call CollectCounts as CountSingletons {
    input:
      vcf = vcf,
      vcf_idx = vcf_idx,
      contigs = GetContigs.contigs,
      bcftools_options = "--include \"FILTER == 'PASS' & INFO/AC == 1\"",
      output_prefix = "singleton",
      docker = sv_pipeline_docker
  }

  call CollectCounts as CountArtifactDels {
    input:
      vcf = vcf,
      vcf_idx = vcf_idx,
      contigs = GetContigs.contigs,
      bcftools_options = "--include \"FILTER == 'PASS' & INFO/SVTYPE == 'DEL' & INFO/SVLEN > 300 & INFO/SVLEN < 1000 & INFO/AF < 0.01\"",
      output_prefix = "artifact",
      docker = sv_pipeline_docker
  }

  call CollectCounts as CountRdOnly {
    input:
      vcf = vcf,
      vcf_idx = vcf_idx,
      contigs = GetContigs.contigs,
      bcftools_options = "--include \"FILTER == 'PASS' & INFO/SVTYPE == 'DEL,DUP' & (INFO/ALGORITHMS == 'depth' | INFO/EVIDENCE == 'RD')\"",
      output_prefix = "depth_only",
      docker = sv_pipeline_docker
  }

  call CollectCounts as CountTssCnvs {
    input:
      vcf = vcf,
      vcf_idx = vcf_idx,
      contigs = GetContigs.contigs,
      bcftools_options = "--include \"FILTER == 'PASS' & INFO/SVTYPE == 'DEL,DUP'\"",
      intersect_file = tss_bed,
      output_prefix = "tss",
      docker = sv_pipeline_docker
  }

  call CollectCounts as CountExonic {
    input:
      vcf = vcf,
      vcf_idx = vcf_idx,
      contigs = GetContigs.contigs,
      bcftools_options = "--include \"FILTER == 'PASS' & INFO/SVTYPE == 'DEL,DUP,INS'\"",
      intersect_file = exons_bed,
      output_prefix = "exonic",
      docker = sv_pipeline_docker
  }

  Array[File] counts = [CountAll.counts, CountRare.counts, CountSingletons.counts,
                        CountArtifactDels.counts, CountRdOnly.counts, 
                        CountTssCnvs.counts, CountExonic.counts]

  call Utils.ConcatTextFiles as MergeCounts {
    input:
      shards = counts,
      sort_command = "sort -Vk1,1 -k2,2V",
      compression_command = "gzip -c",
      input_has_header = false,
      output_filename = prefix + ".all_outlier_counts.tsv.gz",
      docker = sv_pipeline_docker
  }

  output {
    File merged_counts = MergeCounts.merged_file
  }
}


task CollectCounts {
  input {
    File vcf
    File vcf_idx
    Array[String] contigs
    String bcftools_options = ""
    File? intersect_file
    String output_prefix
    String docker
  }

  String bedtools_cmd = if defined(intersect_file) then "| bedtools intersect -u -wa -header -a - -b " + basename(select_first([intersect_file, ""])) else ""
  String outfile = output_prefix + ".counts.tsv"
  Int disk_gb = ceil(3 * size(vcf, "GB")) + 25

  command <<<
    set -eu -o pipefail

    if [ ~{defined(intersect_file)} == "true" ]; then
      mv ~{default="" intersect_file} ./
    fi

    bcftools view \
      ~{bcftools_options} \
      --regions "~{sep=',' contigs}" \
      ~{vcf} \
    ~{bedtools_cmd} \
    | bgzip -c > filtered.vcf.gz
    tabix -p vcf -f filtered.vcf.gz

    svtk count-svtypes \
      --no-header \
      filtered.vcf.gz \
    | awk -v OFS="\t" -v prefix=~{output_prefix} \
      '{ print $1, prefix"_"$2, $3 }' \
    > ~{outfile}
  >>>

  output {
    File counts = "~{outfile}"
  }

  runtime {
    docker: docker
    memory: "3.5 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}