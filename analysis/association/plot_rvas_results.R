#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2024-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Plot Q-Qs for all RVAS summary statistics


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(PedSV, quietly=TRUE)
PedSV::load.constants("all")


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Plot RVAS QQs")
parser$add_argument("--stats", metavar=".bed", type="character",
                    help="RVAS association stats .bed", required=TRUE)
parser$add_argument("--out-prefix", metavar="path", type="character", required=TRUE,
                    help="path/prefix for all output files")
args <- parser$parse_args()

# # DEV:
# args <- list("stats" = "~/scratch/PedSV.v2.5.3.sv_rvas.sumstats.bed.gz",
#              "out_prefix" = "~/scratch/PedSV.v2.5.3.sv_rvas")

# Load summary statistics
stats <- read.table(args$stats, header=T, sep="\t", comment.char="")

# Infer multiple test correction
p.cutoff <- 0.05 / length(unique(stats$gene))

# Iterate over cancers
sapply(names(cancer.colors), function(cancer){

  # Only plot cancers present in header of --stats
  p.colname <- paste(cancer, "neglog10_p", sep=".")
  if(!p.colname %in% colnames(stats)){
    return()
  }

  # One QQ for each consequence
  sapply(unique(stats$consequence), function(csq){
    png(paste(args$out_prefix, cancer, csq, "rvas", "qq", "png", sep="."),
        height=2.25*400, width=2.25*400, res=400)
    par(family="Arial")
    plot.qq(pvals=10^-stats[which(stats$consequence == csq), p.colname],
            pt.color=cancer.colors[cancer], fdr.color=cancer.palettes[[cancer]]["dark2"],
            cutoff=p.cutoff, do.fdr=FALSE)
    dev.off()
  })
})

