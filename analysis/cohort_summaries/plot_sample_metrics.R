#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2024-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Plot various quality control metrics and sample metadata


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(PedSV, quietly=TRUE)
require(argparse, quietly=TRUE)
PedSV::load.constants("all")


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Plot sample QC metrics")
parser$add_argument("metadata", metavar=".tsv", type="character",
                    help="sample metadata .tsv")
parser$add_argument("--subset-samples", metavar=".tsv", type="character",
                    help="list of samples to subset [default: use all samples]")
parser$add_argument("--out-prefix", metavar="path", type="character", required=TRUE,
                    help="path/prefix for all output files")
args <- parser$parse_args()

# # DEV:
# args <- list("metadata" = "~/scratch/PedSV.v2.5.3.cohort_metadata.w_control_assignments.tsv.gz",
#              "subset_samples" = "~/scratch/PedSV.v2.5.3.final_analysis_cohort.samples.list",
#              "out_prefix" = "~/scratch/PedSV.v2.5.3.dev")
# args <- list("metadata" = "/Users/ryan/Desktop/Collins/VanAllen/jackie_younglung/younglung_metadata/YL.SV.v1.1.analysis_metadata.tsv.gz",
#              "subset_samples" = "/Users/ryan/Desktop/Collins/VanAllen/jackie_younglung/YL_analysis/YL.analysis_samples.list",
#              "out_prefix" = "~/scratch/YL.SV.v1.1")

# Load metadata
keepers <- NULL
if(!is.null(args$subset_samples)){
  keepers <- read.table(args$subset_samples, header=F)[, 1]
}
meta <- PedSV::load.sample.metadata(args$metadata, keep.samples=keepers)

# Set plot dimensions
swarmplot.height <- 3.75
swarmplot.width <- 2.2 + (length(unique(meta$disease)) / 2)
swarmplot.height.small <- 2.75
swarmplot.width.small <- 1.25 + (length(unique(meta$disease)) / 2)

# Plot coverage by cancer type
v.cov <- as.numeric(meta$median_coverage)
names(v.cov) <- rownames(meta)
pdf(paste(args$out_prefix, "wgs_coverage_by_cancer.pdf", sep="."),
    height=swarmplot.height, width=swarmplot.width)
swarmplot.by.phenotype(v.cov, meta, median.labels=T, y.axis.title="WGS coverage",
                       shorten.cancer.names=T, add.sample.size=F,
                       sample.size.mirror.buffer=0.15)
dev.off()

# Smaller coverage plot (better viewing on slides)
pdf(paste(args$out_prefix, "wgs_coverage_by_cancer.small.pdf", sep="."),
    height=swarmplot.height.small, width=swarmplot.width.small)
swarmplot.by.phenotype(v.cov, meta, median.labels=T, y.axis.title="WGS coverage",
                       shorten.cancer.names=T, add.sample.size=F,
                       sample.size.mirror.buffer=0.15)
dev.off()
