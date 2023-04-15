#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Plot genome-wide SV count and size distributions for a single cohort


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(PedSV, quietly=TRUE)
require(argparse, quietly=TRUE)
PedSV::load.constants("all")


##################
# Data functions #
##################
# Load sample counts as named vector
load.counts <- function(tsv.in){
  df <- read.table(tsv.in, header=T, sep="\t")
  vals <- as.numeric(df$count)
  names(vals) <- df$sample
  return(vals)
}



######################
# Plotting functions #
######################
# Plot a waterfall of SV counts grouped by ancestry and colored/ordered by disease
plot.waterfall <- function(counts, meta, pop.spacer=0.05){
  # Get plot dimensions
  n.samples <- length(counts)
  n.pops <- length(unique(meta$inferred_ancestry))
  pop.spacer <- ceiling(pop.spacer * n.samples)
  xlims <- c(0, n.samples + (pop.spacer * (n.pops-1)))
  ylims <- c(0, max(counts))

  # Prep plot area
  prep.plot.area(xlims, ylims, parmar=c(2, 2.5, 0.15, 0.15), xaxs="r", yaxs="i")
  clean.axis(2, at=axTicks(2), labels=paste(axTicks(2) / 1000, "k", sep=""),
             title="SVs / genome", infinite=TRUE)

  # Plot ancestries in alphabetical order
  pops <- sort(unique(meta$inferred_ancestry))
  n.plotted <- 0
  for(i in 1:length(pops)){
    pop <- pops[i]
    n.start <- n.plotted
    for(pheno in intersect(names(cancer.colors),
                           unique(metadata.cancer.label.map[meta$disease]))){
      sids <- rownames(meta)[which(meta$inferred_ancestry == pop
                                   & metadata.cancer.label.map[meta$disease] == pheno)]
      vals <- sort(counts[sids], decreasing=T)
      n.vals <- length(vals)
      rect(xleft=n.plotted+(1:n.vals)-1, xright=n.plotted+(1:n.vals),
           ybottom=0, ytop=vals, col=cancer.colors[pheno], border=cancer.colors[pheno])
      segments(x0=n.plotted, x1=n.plotted+n.vals, y0=median(vals), y1=median(vals),
               col=cancer.palettes[[pheno]]["dark1"], lend="butt")
      n.plotted <- n.plotted + n.vals
    }
    axis(1, at=c(n.start, n.plotted), tck=0, labels=NA)
    text(x=mean(c(n.start, n.plotted))+(0.05*diff(par("usr")[1:2])),
         y=par("usr")[3]-(0.075*diff(par("usr")[3:4])),
         pos=2, labels=pop.abbreviations[pop], xpd=T, srt=30)
    n.plotted <- n.plotted + pop.spacer
  }
}



###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Plot SV counts per sample")
parser$add_argument("counts", metavar=".tsv", type="character",
                    help=".tsv of SV counts per sample")
parser$add_argument("--metadata", metavar=".tsv", type="character", required=TRUE,
                    help="sample metadata .tsv")
parser$add_argument("--out-prefix", metavar="path", type="character", required=TRUE,
                    help="path/prefix for all output files")
args <- parser$parse_args()

# # DEV:
# args <- list("metatdata" = "~/scratch/gatk_sv_pediatric_cancers_combined_cohort_metadata_3_31_23.txt",
#              "counts" = "~/scratch/dummy_counts_per_sample.tsv",
#              "out_prefix" = "~/scratch/PedSV.dev")

# Load counts
counts <- load.counts(args$counts)

# Load metadata
meta <- load.sample.metadata(args$metadata, keep.samples=names(counts),
                             reassign.parents=TRUE)

# Plot waterfall
pdf(paste(args$out_prefix, "svs_per_sample.pdf", sep="."),
    height=1.6, width=4)
plot.waterfall(counts, meta)
dev.off()
