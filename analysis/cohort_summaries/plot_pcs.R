#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Plot principal components for manuscript


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
parser <- ArgumentParser(description="Assign sample ancestry")
parser$add_argument("metadata", metavar=".tsv", type="character",
                    help="sample metadata .tsv")
parser$add_argument("--subset-samples", metavar=".tsv", type="character",
                    help="list of samples to subset [default: use all samples]")
parser$add_argument("--out-prefix", metavar="path", type="character", required=TRUE,
                    help="path/prefix for all output files")
args <- parser$parse_args()

# # DEV:
# args <- list("metadata" = "~/scratch/gatk_sv_pediatric_cancers_combined_cohort_metadata_3_31_23.txt",
#              "subset_samples" = "/Users/collins/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/data/ancestry_and_relatedness/PedSV.v1.both_cohorts_final_samples.list",
#              "out_prefix" = "~/scratch/PedSV.dev")

# Load PCs
keepers <- NULL
if(!is.null(args$subset_samples)){
  keepers <- read.table(args$subset_samples, header=F)[, 1]
}
meta <- PedSV::load.sample.metadata(args$metadata, keep.samples=keepers)

# Plot PCs
pdf.dim <- 3
sapply(list(c(1, 2), c(3, 4)), function(pc.idxs){

  # Colored by ancestry
  pdf(paste(args$out_prefix, ".pc", pc.idxs[1], "_vs_pc",
            pc.idxs[2], ".by_ancestry.pdf", sep=""),
      height=pdf.dim, width=pdf.dim)
  pc.scatterplot(meta,
                 paste("PC", pc.idxs[1], sep=""),
                 paste("PC", pc.idxs[2], sep=""),
                 colors=pop.colors[meta$inferred_ancestry],
                 x.title=paste("SV Principal Component", pc.idxs[1]),
                 x.title.line=0.65,
                 y.title=paste("SV Principal Component", pc.idxs[2]),
                 y.title.line=0.65,
                 legend.vals=pop.colors,
                 legend.labels=pop.names.short[names(pop.colors)],
                 parmar=c(2.5, 2.5, 0.5, 0.5))
  dev.off()

  # Colored by cancer type
  pdf(paste(args$out_prefix, ".pc", pc.idxs[1], "_vs_pc",
            pc.idxs[2], ".by_disease.pdf", sep=""),
      height=pdf.dim, width=pdf.dim)
  set.seed(2023)
  pc.scatterplot(meta[sample(1:nrow(meta), nrow(meta), replace=F), ],
                 paste("PC", pc.idxs[1], sep=""),
                 paste("PC", pc.idxs[2], sep=""),
                 colors=cancer.colors[metadata.cancer.label.map[meta$disease]],
                 x.title=paste("SV Principal Component", pc.idxs[1]),
                 x.title.line=0.65,
                 y.title=paste("SV Principal Component", pc.idxs[2]),
                 y.title.line=0.65,
                 legend.vals=cancer.colors[which(names(cancer.colors) != "pancan")],
                 legend.labels=cancer.names.short[setdiff(names(cancer.colors), "pancan")],
                 parmar=c(2.5, 2.5, 0.5, 0.5))
  dev.off()
})

