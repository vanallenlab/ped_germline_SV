#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Collect lists of SV IDs carried by each sample in an allele dosage matrix


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(PedSV, quietly=TRUE)
require(argparse, quietly=TRUE)
require(bedr, quietly=TRUE)


##################
# Data functions #
##################
# Write all SV IDs from a single sample to a file
write.sample.svids.list <- function(ad, sample, outdir){
  svids <- rownames(ad)[which(ad[, sample] != 0)]
  outpath <- paste(gsub("/$", "", outdir), "/", sample, ".svids.list", sep="")
  write.table(svids, outpath, quote=F, sep="\t", col.names=F, row.names=F)
}

###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Get lists of SV IDs per sample")
parser$add_argument("bed", metavar=".tsv", type="character",
                    help="SV sites .bed")
parser$add_argument("ad_matrix", metavar="BED", type="character",
                    help="Allele dosage matrix (.bed/.tsv)")
parser$add_argument("--query", metavar="interval", type="character",
                    help="Tabix-compliant query region.")
parser$add_argument("--outdir", metavar="directory", type="character",
                    default="./", help="output directory for all files")
args <- parser$parse_args()

# # DEV:
# args <- list("ad_matrix" = "~/scratch/PedSV.v1.1.trio_cohort.analysis_samples.wAFs.allele_dosages.bed.gz",
#              "outdir" = "~/scratch/ad_query_test",
#              "query" = "chr22")

# Load BED and extract passing SV IDs
keeper.ids <- rownames(load.sv.bed(args$bed, pass.only=TRUE))

# Format query
query <- as.vector(unlist(sapply(unlist(strsplit(args$query, split=":")), strsplit, split="-")))

# Load AD query into memory
ad <- query.ad.matrix(args$ad_matrix, query.regions=list(c(query)),
                      query.ids=keeper.ids, action="verbose")

# Apply over IDs and write one file of SV IDs per sample
sapply(colnames(ad), write.sample.svids.list, ad=ad, outdir=args$outdir)
