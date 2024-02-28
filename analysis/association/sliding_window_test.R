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
cnvs <- c("DEL", "DUP", "CNV")


##################
# Data functions #
##################
# Load windows to be tested
load.window.counts <- function(tsv.in){
  df <- read.table(tsv.in, header=T, sep="\t", check.names=F, comment.char="")
  df$DELs <- sapply(df$DELs, strsplit, split=",", fixed=T)
  df$DUPs <- sapply(df$DUPs, strsplit, split=",", fixed=T)
  return(df)
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Case:control SV burden tests")
parser$add_argument("--bed", metavar=".tsv", type="character",
                    help="Sliding window .bed.", required=TRUE)
parser$add_argument("--ad", metavar=".tsv", type="character",
                    help="Allele dosage .bed.", required=TRUE)
parser$add_argument("--metadata", metavar=".tsv", type="character",
                    help="sample metadata .tsv", required=TRUE)
parser$add_argument("--subset-samples", metavar=".tsv", type="character",
                    help="list of samples to subset [default: use all samples]")
parser$add_argument("--outfile", metavar="path", type="character", required=TRUE,
                    help=paste("path for output .bed of association stats.",
                               "Will write as uncompressed text."))
args <- parser$parse_args()

# # DEV:
# args <- list("bed" = "~/scratch/chr20.counts.bed.gz",
#              "ad" = "~/scratch/chr20.large_rare_cnvs.ad.bed.gz",
#              "metadata" = "~/scratch/PedSV.v2.5.3.cohort_metadata.w_control_assignments.tsv.gz",
#              "subset_samples" = "~/scratch/PedSV.v2.5.3.final_analysis_cohort.samples.list",
#              "outfile" = "~/scratch/PedSV.v2.5.3.full_cohort.dev.chr20.sliding_window_stats.bed")

# Load BED of windows to test
windows <- load.window.counts(args$bed)

# Load metadata
keepers <- NULL
if(!is.null(args$subset_samples)){
  keepers <- read.table(args$subset_samples, header=F)[, 1]
}
meta <- load.sample.metadata(args$metadata, keep.samples=keepers, reassign.parents=FALSE)

# Load entire AD matrix into memory
# (Pre-filtering recommended; this is already implemented in PedSvMainAnalysis.wdl)
ad <- query.ad.matrix(args$ad, keep.samples=keepers)

# Run association test for each GD with at least one hit
res <- as.data.frame(do.call("rbind", lapply(1:nrow(windows), function(widx){
  # Run association test for DEL, DUP, and CNV separately
  unlist(lapply(cnvs, function(cnv){
    # Run association test for each cancer type
    subres2 <- unlist(lapply(cancers, function(cancer){

      # Run logit glm to test association
      esamps <- get.eligible.samples(meta, cancer)
      Y <- get.phenotype.vector(esamps$cases, esamps$controls)
      del.vids <- unlist(windows$DELs[widx])
      dup.vids <- unlist(windows$DUPs[widx])
      vids <- if(cnv == "DEL"){del.vids}else if(cnv == "DUP"){dup.vids}else{unique(c(del.vids, dup.vids))}
      vids <- vids[which(!is.na(vids))]
      n.sv <- length(vids)
      if(n.sv > 0){
        X <- compress.ad.matrix(ad[vids, ], action="any")
        keepers <- intersect(names(X), names(Y))
        X <- X[keepers]
        Y <- Y[keepers]
        ct <- table(data.frame(X, Y))
        if(nrow(ct) > 1 & ncol(ct) > 1){
          glm.res <- pedsv.glm(meta, X=X, Y=Y, family=binomial(),
                               extra.terms="cohort")
          glm.res[4] <- -log10(as.numeric(glm.res[4]))
        }else{
          glm.res <- c(NA, NA, NA, NA, "glm")
        }
      }else{
        glm.res <- c(NA, NA, NA, NA, "glm")
      }

      # Compute summary metrics
      n.case <- length(esamps$cases)
      n.control <- length(esamps$controls)
      case.nonref <- if(n.sv > 0){length(which(X[esamps$cases] != 0))}else{0}
      case.mean <- if(n.sv > 0){mean(X[esamps$cases], na.rm=T)}else{0}
      control.nonref <- if(n.sv > 0){length(which(X[esamps$controls] != 0))}else{0}
      control.mean <- if(n.sv > 0){mean(X[esamps$controls], na.rm=T)}else{0}

      # Return summary row
      c(n.sv, n.case, case.nonref, case.mean, n.control, control.nonref, control.mean, glm.res)
    }))
  }))
})))
colnames(res) <- as.character(sapply(cnvs, function(cnv){sapply(cancers, function(cancer){
  paste(cancer, cnv,
        c("n_sv", "n_case", "case_nonref", "case_mean",
          "n_control", "control_nonref", "control_mean",
          "beta", "beta_se", "test_stat", "neglog10_p", "model"), sep=".")
})}))

# Merge results with coordinate column and write to --outfile
write.table(as.data.frame(cbind(windows[, 1:3], res)), args$outfile,
            col.names=T, row.names=F, sep="\t", quote=F)

