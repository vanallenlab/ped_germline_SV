#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2024-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Generate Manhattan plots for sliding window association tests


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(PedSV, quietly=TRUE)
PedSV::load.constants("all")


##################
# Data functions #
##################
load.sumstats <- function(tsv.in){
  ss <- read.table(tsv.in, header=T, sep="\t", comment.char="", check.names=F)
  colnames(ss)[1] <- "chrom"
  ss$pos <- sapply(1:nrow(ss), function(i){
    mean(as.numeric(ss[i, c("start", "end")])) + cumsum(c(0, contig.lengths))[which(names(contig.lengths) == ss[i, "chrom"])]
  })
  ss[, c("chrom", "pos", colnames(ss)[grep("neglog10_p", colnames(ss))])]
}


######################
# Plotting functions #
######################
manhattan <- function(stats, p.column, chrom.colors, gw.sig, pt.cex=0.3, title=NULL){
  # Prep plot area
  xlims <- range(stats$pos)
  ylims <- c(0, max(c(stats[, p.column], gw.sig), na.rm=T) + 0.1)
  prep.plot.area(xlims, ylims, parmar=c(2, 2.25, 1, 0.3))
  abline(h=gw.sig, lty=5)

  # Add points
  set.seed(2024)
  stats <- stats[sample(1:nrow(stats), size=nrow(stats)), ]
  points(stats[, c("pos", p.column)], pch=19, cex=pt.cex, col=chrom.colors[stats$chrom])

  # Add axes
  clean.axis(1, at=c(0, cumsum(contig.lengths)), tck=0.025, infinite=TRUE,
             labels=NA, title="Genomic coordinate", title.line=0)
  sapply(1:length(contig.lengths), function(x){
    axis(1, at=((c(0, cumsum(contig.lengths[-24])) + cumsum(contig.lengths))/2)[x],
         tick=F, cex.axis=4/6, labels=gsub("^chr", "", names(contig.lengths)[x]),
         col.axis=chrom.colors[x],
         line=c(-0.8, -1.3)[as.integer((x %% 2) == 1) + 1])
  })
  clean.axis(2, title=bquote(-log[10](italic(P))), infinite.positive=T, title.line=0.2)
  mtext(3, text=title, line=0)
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Plot sliding window Manhattans")
parser$add_argument("--stats", metavar=".bed", type="character",
                    help="sliding window association stats .bed", required=TRUE)
parser$add_argument("--gw-sig", metavar="float", type="numeric", default=10e-8,
                    help="P-value threshold for genome-wide significance")
parser$add_argument("--out-prefix", metavar="path", type="character", required=TRUE,
                    help="path/prefix for all output files")
args <- parser$parse_args()

# # DEV:
# args <- list("stats" = "/Users/ryan/Downloads/PedSV.v2.5.3.case_control_cohort.SlidingWindowResults/PedSV.v2.5.3.case_control_cohort.sliding_window_stats.bed.gz",
#              "gw_sig" = 10e-6,
#              "out_prefix" = "~/scratch/PedSV.v2.5.3.cc.sliding_windows")

# Load summary statistics
stats <- load.sumstats(args$stats)

# Iterate over cancers
sapply(names(cancer.colors), function(cancer){

  # Set cancer palettes
  chrom.colors <- rep(cancer.palettes[[cancer]][c("dark1", "light1")], 12)
  names(chrom.colors) <- names(contig.lengths)

  # Iterate over CNV types
  sapply(c("CNV", "DEL", "DUP"), function(cnv){
    if(cnv == "CNV"){
      title <- paste(cancer.names.long[cancer], "CNVs")
    }else{
      title <- paste(cancer.names.long[cancer], " ", tolower(sv.names[cnv]), "s", sep="")
    }

    # Only plot cancers present in header of --stats
    p.colname <- paste(cancer, cnv, "neglog10_p", sep=".")
    if(!p.colname %in% colnames(stats)){
      return()
    }

    # Plot Manhattan
    png(paste(args$out_prefix, cancer, cnv, "sliding_window", "manhattan", "png", sep="."),
        height=2.25*400, width=3.75*400, res=400)
    manhattan(stats, p.colname, chrom.colors, -log10(args$gw_sig), title=title)
    dev.off()
  })
})

