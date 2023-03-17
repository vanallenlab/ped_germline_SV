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
    String pedsv_docker
  }

  Array[Pair[File, File]] vcf_infos = zip(vcfs, vcf_idxs)

  # 1. Find overlapping samples in each VCF
  scatter ( vcf_info in vcf_infos ) {
    call GetSamples {
      input:
        vcf = vcf_info.left,
        vcf_idx = vcf_info.right,
        docker = pedsv_docker
    }
  }
  call OverlapSamples {
    input:
      sample_lists = GetSamples.samples_list,
      docker = pedsv_docker
  }

  # 2. Subset VCF to overlapping samples prior to clustering
  scatter ( vcf_info in vcf_infos ) {
    call PrepVcfForClustering {
      input:
        vcf = vcf_info.left,
        vcf_idx = vcf_info.right,
        samples_list = OverlapSamples.overlapping_samples,
        out_filename = basename(vcf_info.left, ".vcf.gz") + ".precluster.vcf.gz",
        docker = pedsv_docker
    }
  }

  # 3. Map variant IDs between vcfs using svtk vcfcluster
  call ClusterSites {
    input:
      vcfs = PrepVcfForClustering.vcf_out,
      vcf_idxs = PrepVcfForClustering.vcf_idx_out,
      prefix = prefix,
      docker = pedsv_docker
  }

  # 4. Rename all records to prep for ID-based merge
  scatter ( vcf_info in vcf_infos ) {
    call RenameRecords {
      input:
        vcf = vcf_info.left,
        vcf_idx = vcf_info.right,
        vid_map = ClusterSites.vid_map,
        vid_list = ClusterSites.vid_list,
        out_filename = basename(vcf_info.left, ".vcf.gz") + ".renamed.vcf.gz",
        docker = pedsv_docker
    }
  }

  # 5. Merge VCF genotypes for records that have one (and only one) corresponding match in each cohort
  call MergeVcfs {
    input:
      vcfs = RenameRecords.renamed_vcf,
      vcf_idxs = RenameRecords.renamed_vcf_idx,
      prefix = prefix,
      docker = pedsv_docker
  }

  output {
    File merged_vcf = MergeVcfs.merged_vcf
    File merged_vcf_idx = MergeVcfs.merged_vcf_idx
  }
}


task GetSamples {
  input {
    File vcf
    File vcf_idx
    String docker
  }

  String out_filename = basename(vcf, ".vcf.gz") + ".samples.list"

  command <<<

    set -euo pipefail

    bcftools query -l ~{vcf} > ~{out_filename}

  >>>

  output {
    File samples_list = "~{out_filename}"
  }
  
  runtime {
    cpu: 1
    memory: "2.0 GiB"
    disks: "local-disk " + ceil(size([vcf], "GB")) + 20 + " HDD"
    bootDiskSizeGb: 10
    docker: docker
    preemptible: 3
    maxRetries: 1
  }
}


task OverlapSamples {
  input {
    Array[File] sample_lists
    String docker
  }

  command <<<

    set -euo pipefail

    cp ~{sample_lists[0]} overlapping_samples.list

    while read file; do
      fgrep -wf $file overlapping_samples.list \
      > overlapping_samples.list2
      mv overlapping_samples.list2 overlapping_samples.list
    done < <( sed '1d' ~{write_lines(sample_lists)} )

  >>>

  output {
    File overlapping_samples = "overlapping_samples.list"
  }
  
  runtime {
    cpu: 1
    memory: "2.0 GiB"
    disks: "local-disk  20 HDD"
    bootDiskSizeGb: 10
    docker: docker
    preemptible: 3
    maxRetries: 1
  }
}


task PrepVcfForClustering {
  input {
    File vcf
    File vcf_idx
    File samples_list
    String out_filename
    String docker
  }

  command <<<

    set -euo pipefail

    echo "##source=$( basename ~{vcf} | cut -f1 -d\. )" > header.txt

    bcftools annotate \
      --header-lines header.txt \
      ~{vcf} \
    | bcftools view \
      --force-samples \
      --samples-file ~{samples_list} \
      -O z -o ~{out_filename}

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
    Int cluster_dist = 500
    Float recip_overlap = 0.5
    Float sample_overlap = 0.8
  }

  command <<<

    set -euo pipefail

    # Cluster records
    svtk vcfcluster \
      --dist ~{cluster_dist} \
      --frac ~{recip_overlap} \
      --prefix ~{prefix} \
      --svtypes DEL,DUP,INS,INV,CPX,CTX,BND \
      --sample-overlap ~{sample_overlap} \
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

