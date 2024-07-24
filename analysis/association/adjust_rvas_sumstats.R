#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2024-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Perform saddlepoint approximation on all gene-based SV RVAS summary statistics


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(PedSV, quietly=TRUE)
require(EQL, quietly=TRUE)
PedSV::load.constants("all")


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Correct RVAS summary statistics")
parser$add_argument("--stats", metavar=".bed", type="character",
                    help="RVAS association stats .bed", required=TRUE)
parser$add_argument("--outfile", metavar=".bed", type="character", required=TRUE,
                    help="path to uncompressed output BED")
args <- parser$parse_args()
#
# # DEV:
# args <- list("stats" = "~/scratch/PedSV.v2.5.3.sv_rvas.sumstats.bed.gz",
#              "outfile" = "~/scratch/PedSV.v2.5.3.sv_rvas.sumstats.adj.bed")

# # YL DEV:
# args <- list("stats" = "~/scratch/YL.SV.v1.1.RvasResults/YL.SV.v1.1.sv_rvas.sumstats.bed.gz",
#              "outfile" = "~/scratch/YL.SV.v1.1.RvasResults/YL.SV.v1.1.sv_rvas.sumstats.adj.bed")

# Load summary statistics
stats <- read.table(args$stats, header=T, sep="\t", comment.char="", check.names=F)

# Iterate over cancers
for(cancer in names(cancer.colors)){

  # Only process cancers present in header of --stats
  beta.colname <- paste(cancer, "beta", sep=".")
  if(!beta.colname %in% colnames(stats)){
    next
  }

  # Force all test statistics to be Z-score equivalents (even for FLIC)
  se.colname <- paste(cancer, "beta_se", sep=".")
  zscore.colname <- paste(cancer, "zscore", sep=".")
  stats[, zscore.colname] <- stats[, beta.colname] / stats[, se.colname]

  # Adjust each consequence separately
  for(csq in unique(stats$consequence)){
    csq.idx <- which(stats$consequence == csq)
    adj.vals <- saddlepoint.adj(stats[csq.idx, zscore.colname])
    stats[csq.idx, zscore.colname] <- adj.vals$zscores
    stats[csq.idx, paste(cancer, "neglog10_p", sep=".")] <- -log10(adj.vals$pvalues)
    stats[csq.idx, beta.colname] <- adj.vals$zscores * stats[csq.idx, se.colname]
  }
}

# Write adjusted stats to file
write.table(stats, args$outfile, col.names=T, row.names=F, sep="\t", quote=F)

