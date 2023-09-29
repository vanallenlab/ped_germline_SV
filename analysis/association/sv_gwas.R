#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Conduct association tests for individual SVs in all cancer types


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(PedSV, quietly=TRUE)
require(argparse, quietly=TRUE)
PedSV::load.constants("all")
cancers <- setdiff(names(cancer.colors), "control")


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Individual SV association tests")
parser$add_argument("--ad", metavar=".tsv", type="character", required=TRUE,
                    help="Allele dosage .bed. Required.")
parser$add_argument("--vids", metavar=".txt", type="character",
                    help=paste("List of variant IDs to be tested",
                               "[default: test all variants]"))
parser$add_argument("--metadata", metavar=".tsv", type="character",
                    help="sample metadata .tsv", required=TRUE)
parser$add_argument("--subset-samples", metavar=".tsv", type="character",
                    help="list of samples to subset [default: use all samples]")
parser$add_argument("--out-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="Path to output .tsv of summary statistics. Required.")
args <- parser$parse_args()

# # DEV:
# args <- list("vids" = "~/scratch/PedSV.v2.4.case_control_cohort.chr21.vids.list",
#              "ad" = "~/scratch/PedSV.v2.4.case_control_cohort.analysis_samples.allele_dosages.bed.gz",
#              "metadata" = "~/scratch/PedSV.v2.4.1.cohort_metadata.w_control_assignments.tsv.gz",
#              "subset_samples" = "~/scratch/PedSV.v2.4.case_control_analysis_cohort.samples.list",
#              "out_tsv" = "~/scratch/PedSV.sv_gwas.test.tsv")

# Load metadata
keepers <- NULL
if(!is.null(args$subset_samples)){
  keepers <- read.table(args$subset_samples, header=F)[, 1]
}
meta <- load.sample.metadata(args$metadata, keep.samples=keepers,
                             reassign.parents=FALSE)

# Load allele dosage matrix
if(is.null(args$vids)){
  query.ids <- NULL
}else{
  query.ids <- unique(read.table(args$vids, header=F)[, 1])
}
ad <- query.ad.matrix(args$ad, query.ids=query.ids, keep.samples=rownames(meta))

# Apply over variants
gwas.res <- as.data.frame(do.call("rbind", lapply(rownames(ad), function(vid){

  # Apply over diseases
  c(vid, as.vector(sapply(cancers, function(cancer){

    # Run logit glm to test association
    esamps <- get.eligible.samples(meta, cancer)
    keepers <- colnames(ad)[which(!is.na(ad[vid, ]))]
    esamps$cases <- intersect(esamps$cases, keepers)
    esamps$controls <- intersect(esamps$controls, keepers)
    Y <- get.phenotype.vector(esamps$cases, esamps$controls)
    X <- unlist(ad[vid, c(esamps$cases, esamps$controls)])
    ct <- table(data.frame(X, Y))
    if(nrow(ct) > 1 & ncol(ct) > 1){
      glm.res <- pedsv.glm(meta, X=X, Y=Y, family=binomial())
      glm.res[4] <- -log10(as.numeric(glm.res[4]))
    }else{
      glm.res <- c(NA, NA, NA, NA, "glm")
    }

    # Compute summary metrics
    n.case <- length(esamps$cases)
    case.nonref <- length(which(X[esamps$cases] != 0))
    case.mean <- mean(X[esamps$cases], na.rm=T)
    n.control <- length(esamps$controls)
    control.nonref <- length(which(X[esamps$controls] != 0))
    control.mean <- mean(X[esamps$controls], na.rm=T)

    c(n.case, case.nonref, case.mean, n.control, control.nonref, control.mean, glm.res)
  })))
})))
cancer.cols <- as.character(sapply(cancers, function(cancer){
  paste(cancer, c("n_case", "case_nonref", "case_mean",
                  "n_control", "control_nonref", "control_mean",
                  "beta", "beta_se", "test_stat", "neglog10_p", "model"), sep=".")
}))
colnames(gwas.res) <- c("#variant_id", cancer.cols)

# Write results to output tsv
write.table(gwas.res, args$out_tsv, col.names=T, row.names=F, sep="\t", quote=F)

