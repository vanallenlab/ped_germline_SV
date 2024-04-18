#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2024-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Filter an SV BED to prepare for sliding window large CNV analysis

# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(PedSV, quietly=TRUE)

# Parse command line arguments and options
parser <- ArgumentParser(description="Filter CNVs for sliding window analysis")
parser$add_argument("--bed-in", metavar=".tsv", type="character", required=TRUE,
                    help="SV sites .bed to be filtered. Can be bgzipped.")
parser$add_argument("--minimum-size", metavar="float", type="numeric",
                    help="minimum dosage imbalance [default: 50kb]",
                    default=50000)
parser$add_argument("--bed-out", metavar=".tsv", type="character", required=TRUE,
                    help="Path to filtered .bed. Will be written as uncompressed text.")
args <- parser$parse_args()


# # Dev arguments
# args <- list("bed_in" = "~/scratch/chr1.all_svs.bed.gz",
#              "minimum_size" = 50000,
#              "bed_out" = "~/scratch/chr1.large_rare_cnvs.bed")


# Read SV BED
bed <- load.sv.bed(args$bed_in)

# Filter on frequency
bed <- filter.bed(bed, query="rare")

# Filter on genomic imbalance
bed <- bed[which(calc.genomic.imbalance(bed) >= args$minimum_size), ]

# Sort by coordinate
bed <- bed[with(bed, order(chrom, start, end)), ]

# Write simplified output BED
bed.out <- data.frame(bed[, 1:3], rownames(bed), bed$SVTYPE, bed$CPX_INTERVALS)
colnames(bed.out) <- c("#chrom", "start", "end", "vid", "cnv", "cpx_intervals")
write.table(bed.out, args$bed_out, col.names=T, row.names=F, quote=F, sep="\t")
