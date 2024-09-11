#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Summarize samples lost at various stages of callset creation
# Used to populate REMARK-style supplementary figure


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(PedSV, quietly=TRUE)
require(argparse, quietly=TRUE)
PedSV::load.constants("all")
remark.order <- c("ewing", "neuroblastoma", "osteosarcoma", "parent", "control", "1000G")
onc.carriers <- c("TPMCCDG8038", "PT_3RVRTCYB", "PT_WVTAZ6FJ", "PT_N4YDXR75")


##################
# Data functions #
##################
# Load sample metadata and reformat to minimal information needed for REMARK diagram
load.meta.tsv <- function(tsv.in){
  # Read data
  meta <- read.table(tsv.in, header=T, sep="\t", comment.char="")

  # Set rownames = sample names
  rownames(meta) <- meta$entity.sample_id
  meta$entity.sample_id <- NULL

  # Clean sample grouping for REMARK format
  meta$group <- meta$disease
  meta$group[which(meta$study == "1000G")] <- "1000G"
  parent.idxs <- (meta$disease != "control" & meta$proband != "Yes")
  meta$group[parent.idxs] <- "parent"

  # Only return relevant columns
  keep.cols <- c("global_qc_pass", "batch_qc_pass", "group",
                 colnames(meta)[grep("_(case|control)$", colnames(meta))])
  meta[, keep.cols]
}

# Gather & format information for an individual REMARK cell
remark.cell <- function(meta, prev.meta=NULL){
  k <- sapply(remark.order, function(q){length(which(meta$group == q))})
  k <- c(sum(k), k)
  if(is.null(prev.meta)){
    k.prev <- k
  }else{
    k.prev <- sapply(remark.order, function(q){length(which(prev.meta$group == q))})
    k.prev <- c(sum(k.prev), k.prev)
  }
  delta.v <- 100 * k / k.prev
  delta.v[which(is.na(delta.v))] <- 0
  delta <- sapply(delta.v, function(v){round(v, if(v >= 99.9){0}else if(v >= 1){1}else{2})})
  cat(paste(paste("n=", prettyNum(k, big.mark=","), " (", delta, "%)\n", sep=""), collapse=""))
}

# Simple function to read a list of sample IDs
load.sample.ids <- function(txt.in){
  sort(unique(read.table(txt.in, header=F)[, 1]))
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Gather REMARK diagram data")
parser$add_argument("--metadata", metavar=".tsv", type="character",
                    help="sample metadata .tsv", required=TRUE)
parser$add_argument("--cleanvcf-samples", metavar="list", type="character",
                    help="list of sample IDs present in VCF after CleanVCF",
                    required=TRUE)
parser$add_argument("--final-vcf-samples", metavar="list", type="character",
                    help="list of sample IDs present in VCF after post hoc refinement",
                    required=TRUE)
parser$add_argument("--unrelated-samples", metavar="list", type="character",
                    help="list of unrelated sample IDs to retain for association",
                    required=TRUE)
args <- parser$parse_args()

# # DEV:
# args <- list("metadata" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_v2_callset_generation/v2.5.4/PedSV.v2.5.4.cohort_metadata.w_control_assignments.tsv.gz",
#              "cleanvcf_samples" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/data/remark_misc/PedSV.v2.5.4.samples_in_cleanvcf.list",
#              "final_vcf_samples" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_v2_callset_generation/v2.5.4/PedSV.v2.5.4.full_cohort_w_relatives_1000G.samples.list",
#              "unrelated_samples" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/data/remark_misc/PedSV.v2.5.4.all_samples_no_relatives.list")

# Load input metadata
meta.0 <- load.meta.tsv(args$metadata)

# 0. Summarize initial sample counts
cat("\nGermline WGS harmonized & processed:\n")
remark.cell(meta.0, prev.meta=NULL)

# 1. Global QC
meta.1 <- meta.0[which(meta.0$global_qc_pass == "True"), ]
meta.1l <- meta.0[which(meta.0$global_qc_pass == "False"), ]
cat("\nFailed global quality control:\n")
remark.cell(meta.1l, meta.0)
cat("\nPassed global quality control:\n")
remark.cell(meta.1, meta.0)

# 2. Batch-specific QC
meta.2 <- meta.1[which(meta.1$batch_qc_pass == "True"), ]
meta.2l <- meta.1[which(meta.1$batch_qc_pass == "False"), ]
cat("\nFailed batch-specific quality control:\n")
remark.cell(meta.2l, meta.1)
cat("\nPassed batch-specific quality control:\n")
remark.cell(meta.2, meta.1)

# 3. Post-RF exclusion
cleanvcf.sids <- load.sample.ids(args$cleanvcf_samples)
meta.3 <- meta.2[cleanvcf.sids, ]
meta.3l <- meta.2[setdiff(rownames(meta.2), cleanvcf.sids), ]
cat("\nOutliers after GATK-SV random forest:\n")
remark.cell(meta.3l, meta.2)
cat("\nNot outliers:\n")
remark.cell(meta.3, meta.2)

# 4. Post hoc outlier exclusion
posthoc.survivors <- union(load.sample.ids(args$final_vcf_samples), onc.carriers)
meta.4 <- meta.3[posthoc.survivors, ]
meta.4l <- meta.3[setdiff(rownames(meta.3), posthoc.survivors), ]
cat("\nOutliers after post hoc refinement:\n")
remark.cell(meta.4l, meta.3)
cat("\nNot outliers:\n")
remark.cell(meta.4, meta.3)

# 5. Suspected T-N swap
meta.5 <- meta.4[setdiff(rownames(meta.4), onc.carriers), ]
meta.5l <- meta.4[onc.carriers, ]
cat("\nSuspected somatic contamination:\n")
remark.cell(meta.5l, meta.4)
cat("\nNo evidence of somatic contamination:\n")
remark.cell(meta.5, meta.4)

# 6. Prune relatives + 1000G
unrelateds <- load.sample.ids(args$unrelated_samples)
meta.6 <- meta.5[unrelateds, ]
meta.6l <- meta.5[setdiff(rownames(meta.5), unrelateds), ]
cat("\nRelatives or 1000G:\n")
remark.cell(meta.6l, meta.5)
cat("\nGenetically unrelated:\n")
remark.cell(meta.6, meta.5)

# 7. Ancestry rebalancing
cc.cols <- colnames(meta.6)[grep("_(case|control)$", colnames(meta.6))]
rebal.idxs <- which(apply(meta.6[, cc.cols], 1, function(v){any(v == "True")}))
meta.7 <- meta.6[rebal.idxs, ]
meta.7l <- meta.6[-rebal.idxs, ]
cat("\nExcess from unbalanced ancestries:\n")
remark.cell(meta.7l, meta.6)
cat("\nCohort rebalanced for ancestries:\n")
remark.cell(meta.7, meta.6)
