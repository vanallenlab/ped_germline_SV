#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2024-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Conduct disease association analyses for each of a prespecified list of
# known genomic disorder CNV loci


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(PedSV, quietly=TRUE)
require(Hmisc, quietly=TRUE)
require(argparse, quietly=TRUE)
PedSV::load.constants("all")
cancers <- setdiff(names(cancer.colors), "control")


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Case:control SV burden tests")
parser$add_argument("--bed", metavar=".tsv", type="character", action="append",
                    help=paste("SV sites .bed. Can be supplied one or more",
                               "times. Must be supplied at least once. Order",
                               "must match --ad."), required=TRUE)
parser$add_argument("--ad", metavar=".tsv", type="character", action="append",
                    help=paste("Allele dosage .bed. Can be supplied one or more",
                               "times. Must be supplied at least once. Order must",
                               "match --bed."), required=TRUE)
parser$add_argument("--metadata", metavar=".tsv", type="character",
                    help="sample metadata .tsv", required=TRUE)
parser$add_argument("--genomic-disorder-hits", metavar=".tsv", type="character",
                    required=TRUE, help=paste(".tsv mapping variant IDs to ",
                                              "genomic disorder names"))
parser$add_argument("--subset-samples", metavar=".tsv", type="character",
                    help="list of samples to subset [default: use all samples]")
parser$add_argument("--outfile", metavar="path", type="character", required=TRUE,
                    help="path for output .tsv of association stats")
args <- parser$parse_args()

# # DEV:
# args <- list("bed" = c("~/scratch/PedSV.v2.5.3.full_cohort.analysis_samples.sites.bed.gz"),
#              "ad" = c("~/scratch/PedSV.v2.5.3.full_cohort.analysis_samples.allele_dosages.bed.gz"),
#              "metadata" = "~/scratch/PedSV.v2.5.3.cohort_metadata.w_control_assignments.tsv.gz",
#              "genomic_disorder_hits" = "~/scratch/gd_hits.pairs.tsv",
#              "subset_samples" = "~/scratch/PedSV.v2.5.3.final_analysis_cohort.samples.list",
#              "outfile" = "~/scratch/PedSV.v2.5.3.full_cohort.dev.gd_assoc_results.tsv")

# Load list of IDs to count for each GD
gd.hits <- read.table(args$genomic_disorder_hits, header=F, sep="\t")
colnames(gd.hits) <- c("VID", "GD")
gd.ids <- sort(unique(gd.hits$GD))

# Load BEDs, restrict to variants matching GDs, and pair AD paths with each BED
data <- lapply(1:length(args$bed), function(i){
  bed <- load.sv.bed(args$bed[i])
  bed <- bed[which(rownames(bed) %in% gd.hits$VID), ]
  list("bed" = bed,
       "ad" = query.ad.from.sv.bed(args$ad[i], bed))
})

# Load metadata
keepers <- NULL
if(!is.null(args$subset_samples)){
  keepers <- read.table(args$subset_samples, header=F)[, 1]
}
meta <- load.sample.metadata(args$metadata, keep.samples=keepers, reassign.parents=FALSE)

# Run association test for each GD with at least one hit
gd.res <- as.data.frame(do.call("rbind", lapply(gd.ids, function(gd.id){
  # Run association test for each cancer type
  subres <- unlist(lapply(cancers, function(cancer){

    # Run logit glm to test association
    esamps <- get.eligible.samples(meta, cancer)
    Y <- get.phenotype.vector(esamps$cases, esamps$controls)
    vids <- gd.hits[which(gd.hits$GD == gd.id), "VID"]
    X <- unlist(lapply(data, function(dl){compress.ad.matrix(dl$ad[vids, ], action="any")}))
    keepers <- intersect(names(X), names(Y))
    X <- X[keepers]
    Y <- Y[keepers]
    ct <- table(data.frame(X, Y))
    if(nrow(ct) > 1 & ncol(ct) > 1){
      glm.res <- pedsv.glm(meta, X=X, Y=Y, family=binomial(),
                           extra.terms="cohort", firth.fallback=TRUE)
      glm.res[4] <- -log10(as.numeric(glm.res[4]))
    }else{
      glm.res <- c(NA, NA, NA, NA, "glm")
    }

    # Compute summary metrics
    n.sv <- length(vids)
    n.case <- length(esamps$cases)
    case.nonref <- length(which(X[esamps$cases] != 0))
    case.mean <- mean(X[esamps$cases], na.rm=T)
    n.control <- length(esamps$controls)
    control.nonref <- length(which(X[esamps$controls] != 0))
    control.mean <- mean(X[esamps$controls], na.rm=T)

    c(n.sv, n.case, case.nonref, case.mean, n.control, control.nonref, control.mean, glm.res)
  }))
  c(gd.id, subres)
})))
cancer.cols <- as.character(sapply(cancers, function(cancer){
  paste(cancer, c("n_sv", "n_case", "case_nonref", "case_mean",
                  "n_control", "control_nonref", "control_mean",
                  "beta", "beta_se", "test_stat", "neglog10_p", "model"), sep=".")
}))
colnames(gd.res) <- c("gd.id", cancer.cols)

# Sort results by most significant P-value for any cancer and write to --outfile
sig.order <- order(as.numeric(apply(gd.res[, grep("neglog10_p",colnames(gd.res))],
                                    1, max, na.rm=T)), decreasing=T)
write.table(gd.res[sig.order, ], args$outfile,
            col.names=T, row.names=F, sep="\t", quote=F)

