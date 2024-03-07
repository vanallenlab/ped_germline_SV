#########################
#    Germline SVs in    #
#   Pediatric Cancers   #
#########################
#
# PedSvMainAnalysis.wdl
#
# Generate manuscript-calibre panels of read depth visualization for selected SVs
#
# Copyright (c) 2024-Present, the Van Allen laboratory at Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0


version 1.0


workflow PlotRdPanels {
  input {
    File variant_info_tsv
    String prefix

    Array[File] bincov_matrixes
    Array[File] bincov_matrix_indexes
    Array[File] median_coverage_tsvs

    String gatk_sv_pipeline_docker
    String linux_docker
    String pedsv_r_docker
  }

  call PrepMedianCov {
    input:
      median_cov_tsvs = median_coverage_tsvs,
      prefix = prefix,
      docker = linux_docker
  }

  Array[Array[String]] variant_info_array = read_tsv(variant_info_tsv)

  scatter ( variant_info in variant_info_array ) {
    call PrepBincov {
      input:
        chrom = variant_info[0],
        sv_start = variant_info[1],
        sv_end = variant_info[2],
        svid = variant_info[3],
        bincov_matrixes = bincov_matrixes,
        bincov_matrix_indexes = bincov_matrix_indexes,
        docker = gatk_sv_pipeline_docker
    }

    # call PlotPanel {
    #   chrom = variant_info[0],
    #   sv_start = variant_info[1],
    #   sv_end = variant_info[2],
    #   svid = variant_info[3],
    #   bincov = PrepBincov.matrix,
    #   bincov_idx = PrepBincov.index,
    #   docker = pedsv_r_docker
    # }
  }
}


task PrepMedianCov {
  input {
    Array[File] median_cov_tsvs
    String prefix
    String docker
  }

  Int disk_gb = ceil(size(median_cov_tsvs, "GB")) + 20

  command <<<
    set -eu -o pipefail


    while read file; do
      awk \
        '{ 
              for (i=1; i<=NF; i++)  {
                  a[NR,i] = $i
              }
          }
          NF>p { p = NF }
          END {    
              for(j=1; j<=p; j++) {
                  str=a[1,j]
                  for(i=2; i<=NR; i++){
                      str=str"\t"a[i,j];
                  }
                  print str
              }
          }' \
        $file
    done < ~{write_lines(median_cov_tsvs)} \
    | sort -Vk1,1 -k2,2n \
    | cat <( echo -e "sample_id\tmedian_cov" ) - \
    | gzip -c \
    > ~{prefix}.median_coverage.tsv.gz
  >>>

  output {
    File median_cov = "~{prefix}.median_coverage.tsv.gz"
  }

  runtime {
    docker: docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}


# Slice, merge, and compress bincov matrixes to a query interval of interest
task PrepBincov {
  input {
    String chrom
    Int sv_start
    Int sv_end
    String svid

    Array[String] bincov_matrixes
    Array[File] bincov_matrix_indexes

    Int min_query_pad = 10000
    Int max_bins = 150
    Int disk_gb = 275

    String docker
  }

  Int svlen = ceil(sv_end - sv_start)
  Int interval_length = if ( svlen >= min_query_pad ) then ceil(2 * svlen / 3) else min_query_pad
  Int query_start = if ( sv_start - interval_length < 0 ) then 0 else sv_start - interval_length
  Int query_end = sv_end + interval_length
  String query = chrom + ":" + floor(query_start) + "-" + ceil(query_end)
  Int n_bincovs = length(bincov_matrixes)

  command <<<
    # Relocate all indexes to working directory for tabix read access
    echo "Relocating all bincov tabix indexes to current working directory"
    while read tbi; do
      mv $tbi ./
    done < ~{write_lines(bincov_matrix_indexes)}

    # Export variable required for remote streaming
    export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

    # Loop over all bincov matrixes and slice each in serial
    i=0
    while read bincov; do

      # Slice first bincov matrix to be parent (i.e., keep coordinates)
      if [ $i -lt 1 ]; then
        echo -e "Slicing first bincov matrix ($bincov) to ~{query}"
        tabix -h $bincov ~{query} > bincov_parent.tsv

      # Slice remaining bincov matrixes and only keep coverage values
      else
        echo -e "Slicing additional bincov matrix $i to ~{query}"
        tabix -h $bincov ~{query} | cut -f4- > bincov_slice.$i.tsv

      fi
      echo -e "Slicing success"
      i=$((i+1))
    done < ~{write_lines(bincov_matrixes)}

    # Merge bincov slices
    echo -e "Merging all bincov slices"
    paste bincov_parent.tsv bincov_slice.*.tsv \
    > "~{svid}.~{chrom}_~{query_start}_~{query_end}.rd.raw.bed"
    rm bincov_slice.*.tsv

    # Compress merged bincov if necessary
    n_bins=$( sed '1d' bincov_parent.tsv | wc -l )
    if [ $n_bins -gt ~{max_bins} ]; then
      (( ratio=(n_bins+~{max_bins}-1)/~{max_bins} ))
      echo -e "Compressing merged bincov matrix using ratio $ratio"
      bash /opt/WGD/bin/compressCov.sh \
        -z -o "~{svid}.~{chrom}_~{query_start}_~{query_end}.rd.bed.gz" \
        "~{svid}.~{chrom}_~{query_start}_~{query_end}.rd.raw.bed" \
        $ratio
    else
      cat "~{svid}.~{chrom}_~{query_start}_~{query_end}.rd.raw.bed" \
      | bgzip -c \
      > "~{svid}.~{chrom}_~{query_start}_~{query_end}.rd.bed.gz"
      tabix -p bed -f "~{svid}.~{chrom}_~{query_start}_~{query_end}.rd.bed.gz"
    fi
  >>>

  output {
    File matrix = "~{svid}.~{chrom}_~{query_start}_~{query_end}.rd.bed.gz"
    File index = "~{svid}.~{chrom}_~{query_start}_~{query_end}.rd.bed.gz.tbi"
  }

  runtime {
    docker: docker
    memory: "7.5 GB"
    cpu: 4
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}

