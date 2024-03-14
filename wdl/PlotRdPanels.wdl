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
    File sample_metadata
    File all_eligible_samples_list
    String prefix
    File variant_info_tsv

    Array[File] bincov_matrixes
    Array[File] bincov_matrix_indexes
    Array[File] median_coverage_tsvs

    File gtf
    File gtf_index

    Boolean plot_genes = true
    Boolean add_idiogram = true

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

    call PlotPanel {
      input:
        chrom = variant_info[0],
        sv_start = variant_info[1],
        sv_end = variant_info[2],
        svid = variant_info[3],
        sv_label = variant_info[4],
        sample_metadata = sample_metadata,
        samples_list = all_eligible_samples_list,
        plot_sample_ids = variant_info[5],
        bincov = PrepBincov.matrix,
        bincov_idx = PrepBincov.index,
        cov_medians = PrepMedianCov.median_cov,
        window_start = PrepBincov.window_start,
        window_end = PrepBincov.window_end,
        gtf = gtf,
        gtf_index = gtf_index,
        gene = variant_info[6],
        add_idiogram = add_idiogram,
        prefix = prefix,
        docker = pedsv_r_docker
    }
  }

  call GatherPlots {
    input:
      plots = PlotPanel.plot,
      prefix = prefix,
      docker = linux_docker
  }

  output {
    File plots_tarball = GatherPlots.tarball
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

    Int min_query_pad = 5000
    Int max_query_pad = 1000000
    Int min_notsmall_compression = 8
    Int max_bins = 70
    Int disk_gb = 275

    String docker
  }

  Int svlen = ceil(sv_end - sv_start)
  Int interval_length_pre = if ( svlen >= min_query_pad ) then ceil(2 * svlen / 3) else min_query_pad
  Int interval_length = if ( interval_length_pre > max_query_pad ) then max_query_pad else interval_length_pre
  Int query_start = if ( sv_start - interval_length < 0 ) then 0 else sv_start - interval_length
  Int query_end = sv_end + interval_length
  String query = chrom + ":" + floor(query_start) + "-" + ceil(query_end)
  Int n_bincovs = length(bincov_matrixes)
  Int min_ratio = if ( svlen < 5000 ) then 1 else min_notsmall_compression

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
      if [ $ratio -lt ~{min_ratio} ]; then
        ratio=~{min_ratio}
      fi
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
    Int window_start = query_start
    Int window_end = query_end
  }

  runtime {
    docker: docker
    memory: "7.5 GB"
    cpu: 4
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}


# Plot a single panel
task PlotPanel {
  input {
    String chrom
    Int sv_start
    Int sv_end
    String svid
    String? sv_label

    File sample_metadata
    File samples_list
    String plot_sample_ids      # note: pipe-delimited
    
    File bincov
    File bincov_idx
    File cov_medians
    Int window_start
    Int window_end

    File gtf
    File gtf_index
    String gene                  # set to "NA" if no gene highlight desired

    Boolean plot_genes = true
    Boolean add_idiogram = true

    String prefix
    String docker
  }

  String window_coords = chrom + ":" + window_start + "-" + window_end
  String sv_coords = chrom + ":" + sv_start + "-" + sv_end
  String outfile = prefix + "." + svid + "." + plot_sample_ids + ".rd_viz.pdf"
  Int disk_gb = ceil(3 * size([bincov], "GB")) + 20

  command <<<
    set -eu -o pipefail

    # Prep gene features
    if [ ~{plot_genes} == "true" ]; then
      tabix ~{gtf} ~{window_coords} | bgzip -c > gfeats.gtf.gz
      if [ ~{gene} != "NA" ]; then
        zcat gfeats.gtf.gz | fgrep -w ~{select_first([gene])} | bgzip -c > gfeats.gtf.gz2
        mv gfeats.gtf.gz2 gfeats.gtf.gz
      fi
      zcat gfeats.gtf.gz | awk -v OFS="\t" '{ print $1, $4, $5, $3 }' | bgzip -c > gfeats.bed.gz
    fi

    # Run plot command
    cmd="/opt/ped_germline_SV/analysis/utilities/plot_rd_single_locus.R"
    cmd="$cmd --bincov ~{bincov} --cov-medians ~{cov_medians}"
    cmd="$cmd --metadata ~{sample_metadata} --subset-samples ~{samples_list}"
    n_samps=0
    while read sid; do
      n_samps=$(( n_samps + 1 ))
      cmd="$cmd --sample-id \"$sid\""
    done < <( echo "~{plot_sample_ids}" | sed 's/|/\n/g' )
    if [ $n_samps -gt 1 ]; then
      cmd="$cmd --no-parents"
    fi
    cmd="$cmd --sv-interval ~{sv_coords}"
    if [ ~{defined(sv_label)} == "true" ]; then
      cmd="$cmd --sv-label '~{select_first([sv_label])}'"
    fi
    if [ ~{plot_genes} == "true" ]; then
      if [ ~{gene} != "NA" ]; then
        cmd="$cmd --highlight-gene-label '~{select_first([gene])}'"
        cmd="$cmd --highlight-gene-features gfeats.bed.gz"
      else
        cmd="$cmd --background-gene-features gfeats.bed.gz"
      fi
    fi
    if [ ~{add_idiogram} == "false" ]; then
      cmd="$cmd --no-idiogram"
    fi
    cmd="$cmd --outfile \"~{outfile}\""
    echo -e "Now plotting with the following command:\n$cmd"
    eval $cmd
  >>>

  output {
    File plot = outfile
  }

  runtime {
    docker: docker
    memory: "7.5 GB"
    cpu: 4
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}


task GatherPlots {
  input {
    Array[File] plots
    String prefix
    String docker
  }

  String out_prefix = prefix + ".polished_rd_viz"
  Int disk_gb = ceil(3 * size(plots, "GB")) + 20

  command <<<
    set -eu -o pipefail

    mkdir ~{out_prefix}

    while read file; do
      mv "$file" ~{out_prefix}/
    done < ~{write_lines(plots)}

    tar -czvf ~{out_prefix}.tar.gz ~{out_prefix}
  >>>

  output {
    File tarball = out_prefix + ".tar.gz"
  }

  runtime {
    docker: docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}

