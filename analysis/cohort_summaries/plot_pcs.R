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

# Suffle order of samples to avoid control overplotting, and put controls first
set.seed(2023)
meta <- meta[sample(1:nrow(meta), nrow(meta), replace=F), ]
meta <- meta[c(which(meta$disease == "control"),
               which(meta$disease != "control")), ]

# Plot PCs
pdf.dim <- 3.5
axis.title.line <- 0.65
parmar <- c(2.65, 2.65, 0.5, 0.5)
pop.legend.labels <- sapply(names(pop.colors), function(pop){
  paste(pop.names.short[pop], " (",
        round(100 * length(which(meta$inferred_ancestry == pop)) / nrow(meta), 0),
        "%)", sep="")
})
disease.legend.labels <- sapply(names(cancer.colors), function(disease){
  paste(cancer.names.short[disease], " (",
        round(100 * length(which(metadata.cancer.label.map[meta$disease] == disease)) / nrow(meta), 0),
        "%)", sep="")
})
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
                 x.title.line=axis.title.line,
                 y.title=paste("SV Principal Component", pc.idxs[2]),
                 y.title.line=axis.title.line,
                 legend.vals=pop.colors[which(names(pop.colors) != "OTH")],
                 legend.labels=pop.legend.labels[which(names(pop.colors) != "OTH")],
                 parmar=parmar)
  dev.off()

  # Colored by cancer type
  pdf(paste(args$out_prefix, ".pc", pc.idxs[1], "_vs_pc",
            pc.idxs[2], ".by_disease.pdf", sep=""),
      height=pdf.dim, width=pdf.dim)
  pc.scatterplot(meta,
                 paste("PC", pc.idxs[1], sep=""),
                 paste("PC", pc.idxs[2], sep=""),
                 colors=cancer.colors[metadata.cancer.label.map[meta$disease]],
                 x.title=paste("SV Principal Component", pc.idxs[1]),
                 x.title.line=axis.title.line,
                 y.title=paste("SV Principal Component", pc.idxs[2]),
                 y.title.line=axis.title.line,
                 legend.vals=cancer.colors[which(names(cancer.colors) != "pancan")],
                 legend.labels=disease.legend.labels[setdiff(names(cancer.colors), "pancan")],
                 parmar=parmar)
  dev.off()
})
