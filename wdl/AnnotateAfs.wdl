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


workflow AnnotateAfs {
  input {
    File vcf
    File vcf_idx
    Int records_per_shard
    File? keep_samples_list
    File? ancestry_labels
    File? disease_labels
    File? family_labels
    String? prefix
    Boolean make_bed = false
    Boolean make_ad_matrix = false
    String pedsv_docker
  }

  String prefix = select_first([prefix, basename(vcf, ".vcf.gz") + ".wAFs"])

  call ShardVcf {
    input:
      vcf = vcf,
      vcf_idx = vcf_idx,
      records_per_shard = records_per_shard,
      prefix = basename(vcf, ".vcf.gz") + ".shard",
      docker = pedsv_docker
  }

  scatter ( shard in ShardVcf.vcf_shards ) {
    call PruneAndAddAfs {
      input:
        vcf = shard,
        keep_samples_list = keep_samples_list,
        ancestry_labels = ancestry_labels,
        disease_labels = disease_labels,
        family_labels = family_labels,
        prefix = basename(shard, ".vcf.gz") + ".wAFs",
        docker = pedsv_docker
    }
  }

  call MergeShards {
    input:
      vcfs = PruneAndAddAfs.vcf_out,
      prefix = select_first([prefix, basename(vcf, ".vcf.gz") + ".wAFs"]),
      docker = pedsv_docker
  }

  if ( make_bed ) {
    call Vcf2Bed {
      input:
        vcf = MergeShards.merged_vcf,
        docker = pedsv_docker
    }
  }

  if ( make_ad_matrix ) {
    call MakeAdMatrix {
      input:
        vcf = MergeShards.merged_vcf,
        vcf_idx = MergeShards.merged_vcf_idx,
        docker = pedsv_docker
    }
  }

  output {
    File annotated_vcf = MergeShards.merged_vcf
    File annotated_vcf_idx = MergeShards.merged_vcf_idx
    File? annotated_bed = Vcf2Bed.bed
    File? annotated_bed_idx = Vcf2Bed.bed_idx
    File? ad_matrix = MakeAdMatrix.ad_matrix
    File? ad_matrix_idx = MakeAdMatrix.ad_matrix_idx
  }
}


task ShardVcf {
  input {
    File vcf
    File vcf_idx
    Int records_per_shard
    String prefix
    String docker
  }

  command <<<

    set -eu -o pipefail

    # Estimate number of total shards required
    n_splits=$( bcftools query --format '%ID\n' ~{vcf} | wc -l \
                | awk -v shard_size=~{records_per_shard} '{ print $1 / shard_size }' \
                | cut -f1 -d\. )
    n_digits=${#n_splits}

    # Extract VCF header
    tabix -H ~{vcf} > header.vcf

    # Shard VCF using GNU split
    bcftools view --no-header ~{vcf} > records.vcf
    split -d -n l/$n_splits -a $n_digits \
      --numeric-suffixes=$( printf "%0${n_digits}d" 1 ) \
      records.vcf ~{prefix}

    # Add header to each split and compress
    for records in $( find ./ -name "~{prefix}*" ); do
      cat header.vcf $records | bgzip -c \
      > ${records}.vcf.gz
    done

  >>>

  output {
    Array[File] vcf_shards = glob("~{prefix}*.vcf.gz")
  }

  runtime {
    cpu: 1
    memory: "3.5 GiB"
    disks: "local-disk " + ceil(30 * size(vcf, "GB")) + 20 + " HDD"
    bootDiskSizeGb: 10
    docker: docker
    preemptible: 3
    maxRetries: 1
  }
}


task PruneAndAddAfs {
  input {
    File vcf
    File? keep_samples_list
    File? ancestry_labels
    File? disease_labels
    File? family_labels
    String prefix
    String docker
  }

  String out_filename = prefix + ".vcf.gz"

  command <<<

    set -eu -o pipefail

    # Build command
    cmd="/opt/ped_germline_SV/gatksv_scripts/annotate_afs.py"
    if [ ~{defined(keep_samples_list)} == "true" ]; then
      cmd="$cmd --keep-samples ~{keep_samples_list}"
    fi
    if [ ~{defined(ancestry_labels)} == "true" ]; then
      cmd="$cmd --ancestry-labels ~{ancestry_labels}"
    fi
    if [ ~{defined(disease_labels)} == "true" ]; then
      cmd="$cmd --disease-labels ~{disease_labels}"
    fi
    if [ ~{defined(family_labels)} == "true" ]; then
      cmd="$cmd --family-labels ~{family_labels}"
    fi
    cmd="$cmd --pace ~{vcf} ~{out_filename}"
    echo -e "Annotating VCF with the following command:\n$cmd"
    eval "$cmd"

  >>>

  output {
    File vcf_out = "~{out_filename}"
  }

  runtime {
    cpu: 2
    memory: "7.5 GiB"
    disks: "local-disk " + ceil(4 * size(vcf, "GB")) + 20 + " HDD"
    bootDiskSizeGb: 10
    docker: docker
    preemptible: 3
    maxRetries: 1
  }
}


task MergeShards {
  input {
    Array[File] vcfs
    String prefix
    String docker
  }
  
  String out_filename = prefix + ".vcf.gz"

  command <<<

    set -eu -o pipefail

    bcftools concat \
      --no-version \
      --naive \
      --file-list ~{write_lines(vcfs)} \
    | bcftools view \
      --include 'INFO/SVTYPE == "CNV" | FILTER ~ "MULTIALLELIC" | INFO/AC > 0' \
      -Oz -o ~{out_filename}
    tabix -f ~{out_filename}

  >>>

  output {
    File merged_vcf = "~{out_filename}"
    File merged_vcf_idx = "~{out_filename}.tbi"
  }
  
  runtime {
    cpu: 2
    memory: "3.5 GiB"
    disks: "local-disk " + ceil(4 * size(vcfs, "GB")) + 20 + " HDD"
    bootDiskSizeGb: 10
    docker: docker
    preemptible: 3
    maxRetries: 1
  }
}


task Vcf2Bed {
  input {
    File vcf
    String docker
  }

  String out_filename = basename(vcf, ".vcf.gz") + ".bed.gz"

  command <<<

    set -eu -o pipefail

    svtk vcf2bed \
      --no-samples \
      --info ALL \
      --include-filters \
      ~{vcf} - \
    | bgzip -c > ~{out_filename}

    tabix -p bed -f ~{out_filename}

  >>>

  output {
    File bed = "~{out_filename}"
    File bed_idx = "~{out_filename}.tbi"
  }
  
  runtime {
    cpu: 2
    memory: "7.5 GiB"
    disks: "local-disk " + ceil(3 * size(vcf, "GB")) + 20 + " HDD"
    bootDiskSizeGb: 10
    docker: docker
    preemptible: 3
    maxRetries: 1
  }
}


task MakeAdMatrix {
  input {
    File vcf
    File vcf_idx
    String docker
  }

  String out_filename = basename(vcf, ".vcf.gz") + "allele_dosages.bed.gz"

  command <<<

    set -eu -o pipefail

    # Write header
    tabix -H ~{vcf} | fgrep -v "##" | cut -f10- \
    | paste <( echo -e "#chr\tstart\tend\tID" ) - \
    > header.bed

    # Convert AD matrix for biallelic variants
    bcftools query \
      --exclude 'FILTER ~ "MULTIALLELIC" | INFO/SVTYPE == "CNV"' \
      --format '%CHROM\t%POS\t%INFO/END\t%ID\t[%GT\t]\n' \
      ~{vcf} \
    | sed -e 's/\.\/\./NA/g' -e 's/0\/0/0/g' -e 's/0\/1/1/g' -e 's/1\/1/2/g' \
    > biallelic.ad.bed

    # Convert AD matrix for mCNVs
    bcftools query \
      --include 'FILTER ~ "MULTIALLELIC" | INFO/SVTYPE == "CNV"' \
      --format '%CHROM\t%POS\t%INFO/END\t%ID\t[%RD_CN\t]\n' \
      ~{vcf} \
    | sed -e 's/\t\./\tNA/g' \
    > multiallelic.ad.bed

    # Merge, sort, compress, and index
    cat biallelic.ad.bed multiallelic.ad.bed \
    | sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
    | cat header.bed - \
    | bgzip -c \
    > ~{out_filename}
    tabix -p bed -f ~{out_filename}

  >>>

  output {
    File ad_matrix = "~{out_filename}"
    File ad_matrix_idx = "~{out_filename}.tbi"
  }
  
  runtime {
    cpu: 2
    memory: "3.75 GiB"
    disks: "local-disk " + ceil(10 * size(vcf, "GB")) + 20 + " HDD"
    bootDiskSizeGb: 10
    docker: docker
    preemptible: 3
    maxRetries: 1
  }
}

