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
        prefix = basename(shard, ".vcf.gz"),
        docker = pedsv_docker
    }
  }

  call MergeShards {
    input:
      vcfs = PruneAndAddAfs.vcf_out,
      prefix = prefix,
      docker = pedsv_docker
  }

  output {
    File annotated_vcf = MergeShards.merged_vcf
    File annotated_vcf_idx = MergeShards.merged_vcf_idx
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
    bcftools view --no-header ~{vcf} \
    | split -n l/$n_splits -a $n_digits \
      --numeric-suffixes $(printf "%0${N_DIGITS}d" 1) \
      - ~{prefix}

    # Add header to each split and compress
    for records in $( find ./ -name "~{prefix}" ); do
      cat header.vcf $records | bgzip -c \
      > ${prefix}.vcf.gz
    done

  >>>

  output {
    Array[File] vcf_shards = glob("~{prefix}*.vcf.gz")
  }

  runtime {
    cpu: 1
    memory: "3.5 GiB"
    disks: "local-disk " + ceil(20 * size(vcf, "GB")) + 20 + " HDD"
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

  command <<<

    set -eu -o pipefail

    # Build command
    cmd="/opt/ped_germline_SV/gatksv_scripts/annotate_afs.py"
    if [ ~{defined(keep_samples_list)} ]; then
      cmd="$cmd --keep-samples ~{keep_samples_list}"
    fi
    if [ ~{defined(ancestry_labels)} ]; then
      cmd="$cmd --ancestry-labels ~{ancestry_labels}"
    fi
    if [ ~{defined(disease_labels)} ]; then
      cmd="$cmd --disease-labels ~{disease_labels}"
    fi
    if [ ~{defined(family_labels)} ]; then
      cmd="$cmd --family-labels ~{family_labels}"
    fi
    $CODEDIR/gatksv_scripts/annotate_afs.py \
       $WRKDIR/data/ancestry_and_relatedness/PedSV.v1.trio_cohort_final_samples.list \
      --ancestry-labels $WRKDIR/data/ancestry_and_relatedness/merged/PedSV.merged.ancestry_labels.tsv \
      --disease-labels $WRKDIR/data/PedSV.all_samples.phenotype_labels.tsv \
      --family-labels $WRKDIR/data/ancestry_and_relatedness/PedSV.v1.trio_membership_labels.from_RG.tsv \
      ~/scratch/test.vcf.gz \
      ~/scratch/test.wAFs.vcf.gz

  >>>
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
      -Oz -o ~{output_filename} \
      --file-list ~{write_lines(vcfs)}
    tabix -f ~{out_filename}

  >>>

  output {
    File merged_vcf = "~{out_filename}"
    File merged_vcf_idx = "~{out_filename}.tbi"
  }
  
  runtime {
    cpu: 1
    memory: "3.5 GiB"
    disks: "local-disk " + ceil(4 * size(vcfs, "GB")) + 20 + " HDD"
    bootDiskSizeGb: 10
    docker: docker
    preemptible: 3
    maxRetries: 1
  }  
}

