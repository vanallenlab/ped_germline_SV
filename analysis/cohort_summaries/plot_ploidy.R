#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Plot autosomal and sex ploidy


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(PedSV, quietly=TRUE)
require(argparse, quietly=TRUE)
PedSV::load.constants("all")


##################
# Data functions #
##################
# TBD


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Plot sample ploidy")
parser$add_argument("metadata", metavar=".tsv", type="character",
                    help="sample metadata .tsv")
parser$add_argument("--subset-samples", metavar=".tsv", type="character",
                    help="list of samples to subset [default: use all samples]")
parser$add_argument("--out-prefix", metavar="path", type="character", required=TRUE,
                    help="path/prefix for all output files")
args <- parser$parse_args()

# # DEV:
# args <- list("metadata" = "~/scratch/PedSV.v2.1.cohort_metadata.w_control_assignments.tsv.gz",
#              "subset_samples" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_v2_callset_generation/PedSV.v2.1.final_analysis_cohort.samples.list",
#              "out_prefix" = "~/scratch/PedSV.v2.1.dev")

# Load ploidy
keepers <- NULL
if(!is.null(args$subset_samples)){
  keepers <- read.table(args$subset_samples, header=F)[, 1]
}
meta <- PedSV::load.sample.metadata(args$metadata, keep.samples=keepers)

# Suffle order of samples to avoid control overplotting, and put controls first
set.seed(2023)
meta <- meta[sample(1:nrow(meta), nrow(meta), replace=F), ]
meta <- meta[c(which(meta$disease == "control"),
               which(meta$disease != "control")), ]

# Assign strict aneuploidies
# TODO: implement this

# Write table of sample IDs with aneuploidies
# TODO: implement this

# Run case/control burden test of all aneuploidies
# TODO: implement this

# Plot autosomal ploidy swarms
# TODO: implement this

# Plot sex chromosome 2d scatterplot
# TODO: implement this
