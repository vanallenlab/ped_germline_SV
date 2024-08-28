#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2024-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Plot ploidy sex assignment for manuscript


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
parser <- ArgumentParser(description="Plot sample PCs")
parser$add_argument("metadata", metavar=".tsv", type="character",
                    help="sample metadata .tsv")
parser$add_argument("--subset-samples", metavar=".tsv", type="character",
                    help="list of samples to subset [default: use all samples]")
parser$add_argument("--control-label", metavar="string", type="character", default="Adult controls",
                    help="How should the control-only plot be labeled? [default: Adult controls]")
parser$add_argument("--out-prefix", metavar="path", type="character", required=TRUE,
                    help="path/prefix for all output files")
args <- parser$parse_args()

# # DEV:
# args <- list("metadata" = "~/scratch/PedSV.v2.5.3.cohort_metadata.w_control_assignments.tsv.gz",
#              "subset_samples" = "~/scratch/PedSV.v2.5.3.final_analysis_cohort.samples.list",
#              "control_label" = "Adult controls",
#              "out_prefix" = "~/scratch/PedSV.v2.5.3.dev")
# args <- list("metadata" = "/Users/ryan/Desktop/Collins/VanAllen/jackie_younglung/younglung_metadata/YL.SV.v1.1.analysis_metadata.tsv.gz",
#              "subset_samples" = "/Users/ryan/Desktop/Collins/VanAllen/jackie_younglung/YL_analysis/YL.analysis_samples.list",
#              "control_label" = "Controls",
#              "out_prefix" = "~/scratch/YL.SV.v1.1")

# Load sex ploidy
keepers <- NULL
if(!is.null(args$subset_samples)){
  keepers <- read.table(args$subset_samples, header=F)[, 1]
}
meta <- PedSV::load.sample.metadata(args$metadata, keep.samples=keepers)
meta$inferred_sex[which(!meta$inferred_sex %in% c("MALE", "FEMALE") & meta$sex_aneuploidy)] <- "OTHER"
meta$inferred_sex[which(meta$inferred_sex == "TURNER" & !meta$sex_aneuploidy & meta$chrX_CopyNumber < 1.25)] <- "MALE"
meta$inferred_sex[which(meta$inferred_sex == "TURNER" & !meta$sex_aneuploidy & meta$chrX_CopyNumber > 1.25)] <- "FEMALE"

# # Restrict to samples used for formal association testing
# elig.samples <- unlist(get.eligible.samples(meta, "pancan"))
# meta <- meta[intersect(rownames(meta), elig.samples), ]

# Plot X & Y ploidy colored by sex assignment
pdf.dim <- 2.5
axis.title.line <- 0.45
axis.label.line <- -0.65
x.line.adj <- -0.2
parmar <- c(2.45, 2.45, 0.25, 0.25)
sex.legend.labels <- c(
  paste("XX (", round(100 * length(which(meta$inferred_sex == "FEMALE")) / nrow(meta), 1), "%)", sep=""),
  paste("XY (", round(100 * length(which(meta$inferred_sex == "MALE")) / nrow(meta), 1), "%)", sep=""),
  paste("Other (", round(100 * length(which(meta$inferred_sex == "OTHER")) / nrow(meta), 1), "%)", sep="")
)
pdf(paste(args$out_prefix, ".ploidy_sex_assignment.pdf", sep=""),
    height=pdf.dim, width=pdf.dim)
pc.scatterplot(meta,
               "chrX_CopyNumber",
               "chrY_CopyNumber",
               colors=sex.colors[meta$inferred_sex],
               x.title="chrX ploidy",
               x.title.line=axis.title.line + x.line.adj,
               x.label.line=axis.label.line + x.line.adj,
               y.title="chrY ploidy",
               y.title.line=axis.title.line,
               y.label.line=axis.label.line,
               legend.vals=sex.colors[c("FEMALE", "MALE", "OTHER")],
               legend.labels=sex.legend.labels,
               cex=0.4,
               parmar=parmar)
dev.off()

# Plot once separately for cases
meta.sub <- meta[which(meta$disease != "control"), ]
sex.legend.labels <- c(
  paste("XX (", round(100 * length(which(meta.sub$inferred_sex == "FEMALE")) / nrow(meta.sub), 1), "%)", sep=""),
  paste("XY (", round(100 * length(which(meta.sub$inferred_sex == "MALE")) / nrow(meta.sub), 1), "%)", sep=""),
  paste("Other (", round(100 * length(which(meta.sub$inferred_sex == "OTHER")) / nrow(meta.sub), 1), "%)", sep="")
)
pdf(paste(args$out_prefix, ".ploidy_sex_assignment.cases_only.pdf", sep=""),
    height=pdf.dim, width=pdf.dim)
pc.scatterplot(meta.sub,
               "chrX_CopyNumber",
               "chrY_CopyNumber",
               colors=sex.colors[meta.sub$inferred_sex],
               title=paste("Patients (N=",
                           prettyNum(nrow(meta.sub), big.mark=","),
                           ")", sep=""),
               x.title="chrX ploidy",
               x.title.line=axis.title.line + x.line.adj,
               x.label.line=axis.label.line + x.line.adj,
               y.title="chrY ploidy",
               y.title.line=axis.title.line,
               y.label.line=axis.label.line,
               legend.vals=sex.colors[c("FEMALE", "MALE", "OTHER")],
               legend.labels=sex.legend.labels,
               cex=0.4,
               parmar=parmar + c(0, 0, 1, 1))
dev.off()

# Plot once separately for TOPMed controls
meta.sub <- meta[which(meta$disease == "control" & meta$study %in% c("Topmed_BIOME", "Topmed_MESA", "BioMe")), ]
sex.legend.labels <- c(
  paste("XX (", round(100 * length(which(meta.sub$inferred_sex == "FEMALE")) / nrow(meta.sub), 1), "%)", sep=""),
  paste("XY (", round(100 * length(which(meta.sub$inferred_sex == "MALE")) / nrow(meta.sub), 1), "%)", sep=""),
  paste("Other (", round(100 * length(which(meta.sub$inferred_sex == "OTHER")) / nrow(meta.sub), 1), "%)", sep="")
)
pdf(paste(args$out_prefix, ".ploidy_sex_assignment.topmed_controls_only.pdf", sep=""),
    height=pdf.dim, width=pdf.dim)
pc.scatterplot(meta.sub,
               "chrX_CopyNumber",
               "chrY_CopyNumber",
               colors=sex.colors[meta.sub$inferred_sex],
               title=paste(args$control_label, " (N=",
                           prettyNum(nrow(meta.sub), big.mark=","),
                           ")", sep=""),
               x.title="chrX ploidy",
               x.title.line=axis.title.line + x.line.adj,
               x.label.line=axis.label.line + x.line.adj,
               y.title="chrY ploidy",
               y.title.line=axis.title.line,
               y.label.line=axis.label.line,
               legend.vals=sex.colors[c("FEMALE", "MALE", "OTHER")],
               legend.labels=sex.legend.labels,
               cex=0.4,
               parmar=parmar + c(0, 0, 1, 1))
dev.off()

