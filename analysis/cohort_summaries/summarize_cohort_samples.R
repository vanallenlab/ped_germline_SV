#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Summarize samples by family membership & disease cohort


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(PedSV, quietly=TRUE)
require(argparse, quietly=TRUE)
PedSV::load.constants("names")


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Assign sample ancestry")
parser$add_argument("metadata", metavar=".tsv", type="character",
                    help="sample metadata .tsv")
parser$add_argument("--subset-samples", metavar=".tsv", type="character",
                    help="list of samples to subset [default: use all samples]")
args <- parser$parse_args()

# # DEV:
# args <- list("metadata" = "~/scratch/gatk_sv_pediatric_cancers_combined_cohort_metadata_3_31_23.txt",
#              "subset_samples" = "/Users/collins/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/data/ancestry_and_relatedness/PedSV.v1.trio_cohort_final_samples.list")

# Load metadata
keepers <- NULL
if(!is.null(args$subset_samples)){
  keepers <- read.table(args$subset_samples, header=F)[, 1]
}
meta <- PedSV::load.sample.metadata(args$metadata, keep.samples=keepers)

# Summarize
for(dis in sort(unique(meta$disease))){
  n.pro <- length(which(meta$proband & meta$disease == dis))
  n.fa <- length(which(!meta$proband & meta$disease == dis & meta$chrX_CopyNumber < 1.5))
  n.mo <- length(which(!meta$proband & meta$disease == dis & meta$chrX_CopyNumber >= 1.5))
  n.case <- length(which(is.na(meta$proband) & meta$disease == dis & is.na(meta$family_id)))
  cat(paste(dis, ":\n  - ", prettyNum(n.fa, big.mark=","),
            " fathers\n  - ", prettyNum(n.mo, big.mark=","),
            " mothers\n  - ", prettyNum(n.pro, big.mark=","),
            " probands\n  - ", prettyNum(n.case, big.mark=","),
            " isolated individuals\n\n", sep=""))
}
