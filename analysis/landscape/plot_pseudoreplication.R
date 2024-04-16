#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Compose internal pseudo-replication figure panels for selected hypotheses


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(Hmisc, quietly=TRUE)
require(PedSV, quietly=TRUE)
PedSV::load.constants("all")


##################
# Data functions #
##################
# Load ncCWAS sumstats
load.nccwas.ss <- function(tsv.in, cancer){
  ss <- read.table(tsv.in, header=F, sep="\t", quote="")
  colnames(ss) <- c("coefficient", "std.err", "test.statistic", "P.value",
                    "case_count", "case_max",  "n.case.ref", "n.case",
                    "control_count", "control_max",  "n.control.ref", "n.control",
                    "n.SVs", "hypothesis")
  ss$n.case.carrier <- ss$n.case - ss$n.case.ref
  ss$n.control.carrier <- ss$n.control - ss$n.control.ref
  ss$test <- "glm"
  ss$model <- "logit"
  ss$case.mean <- ss$case_count / ss$n.case
  ss$case.stdev <- NA
  ss$control.mean <- ss$control_count / ss$n.control
  ss$control.stdev <- NA
  ss$disease <- cancer
  ss[, c("hypothesis", "disease", "model", "n.case", "n.case.carrier",
         "case.mean", "case.stdev", "n.control", "n.control.carrier",
         "control.mean", "control.stdev", "coefficient", "std.err",
         "test.statistic", "P.value", "test")]
}


# Infer cohort pairings based on study stage & cancer type
get.cohort.names <- function(meta, cancer, stage){
  stage.meta <- meta[which(meta$study_phase == stage), ]
  case.count <- table(stage.meta$study[which(cancer == metadata.cancer.label.map[stage.meta$disease])])
  case.cohort <- cohort.names[tail(names(sort(case.count)), 1)]
  control.count <- table(stage.meta$study[which(stage.meta$disease == "control")])
  control.cohort <- cohort.names[tail(names(sort(control.count)), 1)]
  return(list("case.cohort" = case.cohort, "control.cohort" = control.cohort))
}

# Format data for a single sub-panel
get.subpanel.data <- function(cc.ss, trio.ss, hypothesis, cancer, meta, ci.mode="binomial"){
  # Merge data across study arms
  plot.df <- as.data.frame(do.call("rbind", lapply(list(cc.ss, trio.ss), function(df){
    df[which(df$hypothesis == hypothesis & df$disease == cancer), ]
  })))
  n.case <- plot.df$n.case
  n.control <- plot.df$n.control
  plot.df <- stats2barplotdf(plot.df, ci.mode=ci.mode, order.by.cancer=FALSE)

  # Add descriptive labels for plotting
  rownames(plot.df) <- c("cc", "trio")
  cc.cohort.names <- get.cohort.names(meta, cancer, "case_control")
  trio.cohort.names <- get.cohort.names(meta, cancer, "trio")
  plot.df$label <- c(paste(cc.cohort.names$case.cohort, " (N=",
                           prettyNum(n.case[1], big.mark=","),
                           ")\n", cc.cohort.names$control.cohort, " (N=",
                           prettyNum(n.control[1], big.mark=","),
                           ")", sep=""),
                     paste(trio.cohort.names$case.cohort, " (N=",
                           prettyNum(n.case[2], big.mark=","),
                           ")\n", trio.cohort.names$control.cohort, " (N=",
                           prettyNum(n.control[2], big.mark=","),
                           ")", sep=""))
  plot.df[, c("case.value", "case.ci.lower", "case.ci.upper",
              "control.value", "control.ci.lower", "control.ci.upper",
              "p", "label")]
}


######################
# Plotting functions #
######################
# Plot two sets of two horizontal bars
pseudorep.bars <- function(top.data, top.cancer.types, bottom.data, bottom.cancer.types, panel.spacer=0.5, title=NULL, top.axis.units=NULL, outer.y.axis.labels=NULL, outer.y.axis.line=5.5, parmar=c(0, 8.5, 2, 4.25)){
  # Get plot dimensions
  xlims <- c(0, (4/3) * max(as.numeric(unlist(rbind(top.data, bottom.data)[, c(1, 4)]))))
  ylims <- c(nrow(rbind(top.data, bottom.data)) + panel.spacer, 0)

  # Get bar mids
  top.bar.mids <- (1:nrow(top.data))-0.5
  bottom.bar.mids <- (1:nrow(bottom.data))+panel.spacer+nrow(top.data)-0.5

  # Prep plot area
  prep.plot.area(xlims, ylims, yaxs="r", parmar=parmar)

  # Add bars
  add.pheno.bars(top.data, bar.mids=top.bar.mids, bar.hex=0.45, control.sep=0.175,
                 cancer.types.override=top.cancer.types, add.pvals=TRUE)
  add.pheno.bars(bottom.data, bar.mids=bottom.bar.mids, bar.hex=0.45, control.sep=0.175,
                 cancer.types.override=bottom.cancer.types, add.pvals=TRUE)

  # Add top X axis
  clean.axis(3, label.units=top.axis.units, title=title, infinite.positive=TRUE,
             label.line=-0.8, title.line=0.1, max.ticks=5)

  # Add inner Y labels
  axis(2, at=top.bar.mids, tick=F, las=2, cex.axis=4.5/6, line=-0.8, labels=top.data$label)
  axis(2, at=bottom.bar.mids, tick=F, las=2, cex.axis=4.5/6, line=-0.8, labels=bottom.data$label)

  # Add outer Y labels, if optioned
  if(!is.null(outer.y.axis.labels)){
    axis(2, tck=0, at=range(top.bar.mids)+c(-0.5, 0.5), labels=NA, line=outer.y.axis.line)
    axis(2, at=mean(top.bar.mids), tick=F, line=outer.y.axis.line-0.9, las=2, labels=outer.y.axis.labels[1])
    axis(2, tck=0, at=range(bottom.bar.mids)+c(-0.5, 0.5), labels=NA, line=outer.y.axis.line)
    axis(2, at=mean(bottom.bar.mids), tick=F, line=outer.y.axis.line-0.9, las=2, labels=outer.y.axis.labels[2])
  }
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Generate pseudoreplication plots")
parser$add_argument("--case-control-stats", metavar=".tsv", type="character",
                    help="Summary stats from case-control cohort.", required=TRUE)
parser$add_argument("--trio-stats", metavar=".tsv", type="character",
                    help="Summary stats from trio cohort.", required=TRUE)
parser$add_argument("--metadata", metavar=".tsv", type="character",
                    help="sample metadata .tsv", required=TRUE)
parser$add_argument("--case-control-nbl-noncoding-cwas", required=TRUE,
                    help=paste("Neuroblastoma noncoding CWAS summary statistics",
                               "for case-control cohort"))
parser$add_argument("--trio-nbl-noncoding-cwas", required=TRUE,
                    help=paste("Neuroblastoma noncoding CWAS summary statistics",
                               "for trio cohort"))
parser$add_argument("--subset-samples", metavar=".tsv", type="character",
                    help="list of samples to subset [default: use all samples]")
parser$add_argument("--out-prefix", metavar="path", type="character", required=TRUE,
                    help="path/prefix for all output files")
args <- parser$parse_args()

# # DEV:
# args <- list("case_control_stats" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_MS/PedSV_figures/PedSV.v2.5.3_analysis_outputs/stats/PedSV.v2.5.3.case_control_cohort.global_burden_tests.tsv.gz",
#              "trio_stats" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_MS/PedSV_figures/PedSV.v2.5.3_analysis_outputs/stats/PedSV.v2.5.3.trio_cohort.global_burden_tests.tsv.gz",
#              "metadata" = "~/scratch/PedSV.v2.5.3.cohort_metadata.w_control_assignments.tsv.gz",
#              "case_control_nbl_noncoding_cwas" = "/Users/ryan/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/data/ncCWAS_sumstats/adjusted_rlc_04_12_24/neuroblastoma_discovery_noncoding_cwas_concatenated_glm_results_4_11_24.adjusted.txt",
#              "trio_nbl_noncoding_cwas" = "/Users/ryan/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/data/ncCWAS_sumstats/adjusted_rlc_04_12_24/neuroblastoma_trio_noncoding_cwas_concatenated_glm_results_4_11_24.adjusted.txt",
#              "subset_samples" = "~/scratch/PedSV.v2.5.3.final_analysis_cohort.samples.list",
#              "out_prefix" = "~/scratch/PedSV.v2.5.3.dev")

# Load metadata
keepers <- NULL
if(!is.null(args$subset_samples)){
  keepers <- read.table(args$subset_samples, header=F)[, 1]
}
meta <- load.sample.metadata(args$metadata, keep.samples=keepers)

# Load summary statistics for each cohort
cc.ss <- load.burden.sumstats(args$case_control_stats)
trio.ss <- load.burden.sumstats(args$trio_stats)
do.nccwas <- FALSE
if(!is.null(args$case_control_nbl_noncoding_cwas)
   & !is.null(args$trio_nbl_noncoding_cwas)){
  do.nccwas <- TRUE
  cc.nc.ss <- load.nccwas.ss(args$case_control_nbl_noncoding_cwas, cancer="NBL")
  trio.nc.ss <- load.nccwas.ss(args$trio_nbl_noncoding_cwas, cancer="NBL")
}


# Set plot device parameters
pdf.height <- 2.15
pdf.width <- 3.75


# Rate of autosomal, large, rare, unbalanced SVs in NBL males vs. females
pdf(paste(args$out_prefix, "large_rare_unbalanced.male_vs_female.NBL.pseudoreplication.pdf", sep="."),
    height=pdf.height, width=pdf.width)
pseudorep.bars(top.data=get.subpanel.data(cc.ss, trio.ss, "large.rare.unbalanced.autosomal_only.MALE_only", "NBL", meta),
               top.cancer.types=c("NBL", "NBL"),
               bottom.data=get.subpanel.data(cc.ss, trio.ss, "large.rare.unbalanced.autosomal_only.FEMALE_only", "NBL", meta),
               bottom.cancer.types=c("NBL", "NBL"),
               outer.y.axis.labels=c("Male\n(XY)", "Female\n(XX)"),
               title="Samples w/ rare unbalanced SV >1Mb", top.axis.units="percent")
dev.off()


# Rate of singleton gene-disruptive SVs by cancer type
pdf(paste(args$out_prefix, "singleton_gene_disruptive.pseudoreplication.pdf", sep="."),
    height=pdf.height, width=pdf.width)
pseudorep.bars(top.data=get.subpanel.data(cc.ss, trio.ss, "singleton.gene_disruptive", "NBL", meta, ci.mode="normal"),
               top.cancer.types=c("NBL", "NBL"),
               bottom.data=get.subpanel.data(cc.ss, trio.ss, "singleton.gene_disruptive", "EWS", meta, ci.mode="normal"),
               bottom.cancer.types=c("EWS", "EWS"),
               outer.y.axis.labels=c("Neuro.", "Ewing"),
               title="Singleton gene-disruptive SVs/sample")
dev.off()


# NBL noncoding TAD boundary categories
if(do.nccwas){
  pdf(paste(args$out_prefix, "singleton_tad_boundary.ncCWAS.pseudoreplication.pdf", sep="."),
      height=pdf.height, width=pdf.width)
  pseudorep.bars(top.data=get.subpanel.data(cc.nc.ss, trio.nc.ss, "ANY.SINGLETON.ANY.neuroblastoma_tad_boundary.ANY.ANY.ANY.protein_coding", "NBL", meta),
                 top.cancer.types=c("NBL", "NBL"),
                 bottom.data=get.subpanel.data(cc.nc.ss, trio.nc.ss, "DEL.SINGLETON.ANY.neuroblastoma_tad_boundary.ANY.ANY.ANY.protein_coding", "NBL", meta),
                 bottom.cancer.types=c("NBL", "NBL"),
                 outer.y.axis.labels=c("All\nnoncoding\nSVs", "Noncoding\ndeletions"),
                 title="Singletons disrupting adrenal TADs",
                 parmar=c(0, 9.65, 2, 4.25))
  dev.off()
}



