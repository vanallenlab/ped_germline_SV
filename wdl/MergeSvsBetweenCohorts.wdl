#########################
#    Germline SVs in    #
#   Pediatric Cancers   #
#########################
#
# MergeSvsBetweenCohorts.wdl
#
# Generate a map of SV IDs between two independent GATK-SV callsets
#
# Copyright (c) 2023-Present, the Van Allen laboratory at Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0


version 1.0


workflow MergeSvsBetweenCohorts {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    String prefix
    String sv_pipeline_docker
  }

  # 1. Create a sites VCF for each input
  scatter ( vcf_info in zip(vcfs, vcf_idxs) ) {
    call MakeSites {
      input:
        vcf = vcf_info.left,
        vcf_idx = vcf_info.right,
        out_filename = basename(vcf_info.left, ".vcf.gz") + ".sites.vcf.gz",
        docker = sv_pipeline_docker
    }
  }

  # 2. Map variant IDs between vcfs using svtk vcfcluster
  call ClusterSites {
    input:
      vcfs = MakeSites.vcf_out,
      vcf_idxs = MakeSites.vcf_idx_out,
      prefix = prefix,
      docker = sv_pipeline_docker
  }

  # 3. Rename all records to prep for ID-based merge
  scatter ( vcf_info in zip(vcfs, vcf_idxs) ) {
    call RenameRecords {
      input:
        vcf = vcf_info.left,
        vcf_idx = vcf_info.right,
        vid_map = ClusterSites.vid_map,
        vid_list = ClusterSites.vid_list,
        out_filename = basename(vcf_info.left, ".vcf.gz") + ".renamed.vcf.gz",
        docker = sv_pipeline_docker
    }
  }

  # 4. Merge VCF genotypes for records that have one (and only one) corresponding match in each cohort
  call MergeVcfs {
    input:
      vcfs = RenameRecords.renamed_vcf,
      vcf_idxs = RenameRecords.renamed_vcf_idx,
      prefix = prefix,
      docker = sv_pipeline_docker
  }

  output {
    File merged_vcf = MergeVcfs.merged_vcf
    File merged_vcf_idx = MergeVcfs.merged_vcf_idx
  }
}


task MakeSites {
  input {
    File vcf
    File vcf_idx
    String out_filename
    String docker
  }

  command <<<

    set -euo pipefail

    echo "##source=$( basename ~{vcf} | cut -f1 -d\. )" > header.txt

    bcftools annotate \
      --header-lines header.txt \
      ~{vcf} \
    | cut -f1-8 \
    | awk -v FS="\t" -v OFS="\t" \
      '{ if ($1 ~ "^##") print $0; \
         else if ($1 ~ "#") print $0, "FORMAT", "dummy"; \
         else print $0, "GT", "1/1" }' | bgzip -c \
    > ~{out_filename}

    tabix -p vcf -f ~{out_filename}

  >>>

  output {
    File vcf_out = "~{out_filename}"
    File vcf_idx_out = "~{out_filename}.tbi"
  }
  
  runtime {
    cpu: 1
    memory: "3.5 GiB"
    disks: "local-disk " + ceil(3 * size([vcf], "GB")) + 20 + " HDD"
    bootDiskSizeGb: 10
    docker: docker
    preemptible: 3
    maxRetries: 1
  }
}


task ClusterSites {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    String prefix
    String docker
    Int cluster_dist = 250
    Float recip_overlap = 0.8
  }

  command <<<

    set -euo pipefail

    # Cluster records
    svtk vcfcluster \
      --dist ~{cluster_dist} \
      --frac ~{recip_overlap} \
      --prefix ~{prefix} \
      --svtypes DEL,DUP,INS,INV,CPX,CTX,BND \
      --sample-overlap 0 \
      --preserve-ids \
      ~{write_lines(vcfs)} \
      /dev/stdout \
    | cut -f1-8 \
    | bcftools sort -O z -o ~{prefix}.clustered.sites.vcf.gz

    tabix -p vcf -f ~{prefix}.clustered.sites.vcf.gz

    # Extract members comprising clusters of exactly two records
    bcftools query \
      --format '%INFO/MEMBERS\n' \
      ~{prefix}.clustered.sites.vcf.gz \
    | sed 's/,/\t/g' | awk '{ if (NF==2) print $0 }' \
    > ~{prefix}.VID_pairs.prefilter.tsv

    # Filter candidate variant pairs s/t first and second members of pair have
    # different prefixes
    paste \
      <( cut -f1 -d_ ~{prefix}.VID_pairs.prefilter.tsv ) \
      <( cut -f2 ~{prefix}.VID_pairs.prefilter.tsv | cut -f1 -d_ ) \
      ~{prefix}.VID_pairs.prefilter.tsv \
    | awk -v FS="\t" -v OFS="\t" '{ if ($1 != $2) print $3, $4 }' \
    > ~{prefix}.VID_pairs.tsv

    # Make variant ID remapping table
    awk -v FS="\t" -v OFS="\t" '{ print $0, "merged_sv_"NR }' ~{prefix}.VID_pairs.tsv \
    | awk -v FS="\t" -v OFS="\t" '{ print $1, $3"\n"$2, $3 }' \
    > ~{prefix}.vid_remap.tsv

    # Make list of all variant IDs present in any qualifying pair
    sed 's/\t/\n/g' ~{prefix}.VID_pairs.tsv > ~{prefix}.all_VIDs.list

  >>>

  output {
    File clustered_vcf = "~{prefix}.clustered.sites.vcf.gz"
    File clustered_vcf_idx = "~{prefix}.clustered.sites.vcf.gz.tbi"
    File vid_pairs = "~{prefix}.VID_pairs.tsv"
    File vid_map = "~{prefix}.vid_remap.tsv"
    File vid_list = "~{prefix}.all_VIDs.list"
  }
  
  runtime {
    cpu: 1
    memory: "7.5 GiB"
    disks: "local-disk " + ceil(2 * size(vcfs, "GB")) + 20 + " HDD"
    bootDiskSizeGb: 10
    docker: docker
    preemptible: 3
    maxRetries: 1
  }
}


task RenameRecords {
  input {
    File vcf
    File vcf_idx
    File vid_map
    File vid_list
    String out_filename
    String docker
  }

  command <<<

    set -eu -o pipefail

    /opt/ped_germline_SV/gatksv_scripts/rename_records.py \
      --vcf-in ~{vcf} \
      --vid-map ~{vid_map} \
      --vcf-out ~{out_filename}

    tabix -f ~{out_filename}

  >>>

  output {
    File renamed_vcf = "~{out_filename}"
    File renamed_vcf_idx = "~{out_filename}.tbi"
  }
  
  runtime {
    cpu: 2
    memory: "7.5 GiB"
    disks: "local-disk " + ceil(4 * size([vcf], "GB")) + 20 + " HDD"
    bootDiskSizeGb: 15
    docker: docker
    preemptible: 3
    maxRetries: 1
  }
}


task MergeVcfs {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    String prefix
    String docker
  }
  String out_filename = prefix + ".merged.vcf.gz"

  command <<<

    set -eu -o pipefail

    bcftools merge -m id --force-samples ~{sep=" " vcfs} \
    | bcftools sort -O z -o ~{out_filename} -

    tabix -f ~{out_filename}

  >>>

  output {
    File merged_vcf = "~{out_filename}"
    File merged_vcf_idx = "~{out_filename}.tbi"
  }
  
  runtime {
    cpu: 1
    memory: "7.5 GiB"
    disks: "local-disk " + ceil(4 * size(vcfs, "GB")) + 20 + " HDD"
    bootDiskSizeGb: 10
    docker: docker
    preemptible: 3
    maxRetries: 1
  }  
}

