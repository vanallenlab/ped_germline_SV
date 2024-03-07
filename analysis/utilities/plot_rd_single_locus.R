#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2024-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Generate a publication-quality read depth visualization of a locus


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
# Load coverage data and normalize by sample median
load.coverage <- function(cov_in, medians_in=NULL, keep.samples=NULL){
  # Read matrix
  cov <- read.table(cov_in, sep="\t", comment.char="", header=T)
  cov[, 2:ncol(cov)] <- apply(cov[, 2:ncol(cov)], 2, as.numeric)
  colnames(cov)[1:3] <- c("chrom", "start", "end")

  # Median normalization & absolute copy number estimation
  if(!is.null(medians_in)){
    # Read & reformat coverage medians
    cov.meds <- read.table(medians_in, header=T, comment.char="")
    cov.med.samples <- cov.meds[, 1]
    cov.meds <- as.numeric(as.vector(cov.meds[, 2]))
    names(cov.meds) <- cov.med.samples

    # Restrict to samples present in both coverage matrix and median file
    samples.present <- sort(intersect(names(cov.meds), colnames(cov)))
    if(is.null(keep.samples)){
      keep.samples <- samples.present
    }else{
      keep.samples <- intersect(keep.samples, samples.present)
    }
    cn <- do.call("cbind", lapply(keep.samples, function(sid){
      2 * cov[, sid] / cov.meds[sid]
    }))
    colnames(cn) <- keep.samples
    cov <- as.data.frame(cbind(cov[, 1:3], cn))
  }

  return(cov)
}


######################
# Plotting functions #
######################
# Convert x, y pairs into step function
step.function <- function(x, y){
  x.step <- c(x[1], PedSV::stretch.vector(x, 2))
  y.step <- c(y[1], y[1], PedSV::stretch.vector(y, 2))[1:length(x.step)]
  data.frame("x" = x.step, y=y.step)
}

# Add shaded stepwise normal population variability in coverage
shade.middle.quantile <- function(cov, lower, upper, color){
  plot.vals <- step.function(c(cov$start, rev(cov$start)),
                             c(apply(cov[, -c(1:3)], 1, quantile, prob=lower, na.rm=T),
                               apply(cov[, -c(1:3)], 1, quantile, prob=upper, na.rm=T)))
  polygon(plot.vals, border=NA, col=color)
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Plot read depth for a single locus")
parser$add_argument("--bincov", metavar=".bed", type="character",
                    help="Coverage file to be visualized", required=TRUE)
parser$add_argument("--cov-medians", metavar=".tsv", type="character",
                    help="Two-column tsv mapping sample ID to median coverage", required=TRUE)
parser$add_argument("--metadata", metavar=".tsv", type="character",
                    help="Sample metadata .tsv", required=TRUE)
parser$add_argument("--sample-id", metavar="string", type="character",
                    help="IDs for samples to highlight on panel.")
parser$add_argument("--subset-samples", metavar=".tsv", type="character",
                    help="List of samples to subset [default: use all samples]")
parser$add_argument("--sv-interval", metavar="chrom:start-end", type="character",
                    help="SV coordinates to highlight")
parser$add_argument("--sv-label", metavar="string", type="character", default="SV",
                    help="Label for --sv-interval [default: 'SV']")
parser$add_argument("--gene-features", metavar=".bed", type="character",
                    help="BED file of gene features to plot.")
parser$add_argument("--gene-label", metavar="string", type="character",
                    help="Label for gene track [default: 'Gene']")
parser$add_argument("--no-parents", action="store_true",
                    help="Disable plotting of parents for complete trios")
parser$add_argument("--outfile", metavar="path", type="character", required=TRUE,
                    help="Path to output .pdf")
args <- parser$parse_args()

# # DEV
# args <- list("bincov" = "~/scratch/All_20_Batches.chr2.final_cleanup_DUP_chr2_812.chr2_14753993_17131946.rd.bed.gz",
#              "cov_medians" = "~/scratch/PedSV.v2.5.3.median_coverage.tsv.gz",
#              "metadata" = "~/scratch/PedSV.v2.5.3.cohort_metadata.w_control_assignments.tsv.gz",
#              "sample_id" = "PT_V1Q9W1NW",
#              "subset_samples" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_v2_callset_generation/v2.5.3/PedSV.v2.5.3.full_cohort_w_relatives_1000G.samples.list",
#              "sv_interval" = "chr2:15546644-16339295",
#              "sv_label" = "793kb Duplication",
#              "gene_features" = "~/scratch/MYCN.features.bed",
#              "gene_label" = "MYCN",
#              "no_parents" = FALSE,
#              "outfile" = "~/scratch/MYCN.DUP.rdviz.test.pdf")

# Load sample metadata
meta <- load.sample.metadata(args$metadata)

# Read & normalize read depth data
cov <- load.coverage(args$bincov, args$cov_medians)

# Check if both parents are present; if so, prepend these to args$sample_id
samples.to.plot <- args$sample_id
sample.lwds <- rep(3, times=length(samples.to.plot))
sample.colors <- cancer.colors[metadata.cancer.label.map[meta[samples.to.plot, "disease"]]]
if(!args$no_parents){
  for(sid in args$sample_id){
    fam.id <- meta[sid, "family_id"]
    if(is.na(fam.id)){
      next
    }
    parents <- rownames(meta)[which(meta$family_id == fam.id & rownames(meta) != sid)]
    n.parents <- length(parents)
    if(n.parents == 2 & all(parents %in% colnames(cov))){
      samples.to.plot <- c(parents, samples.to.plot)
      sample.lwds <- c(2, 2, sample.lwds)
      sample.colors <- c(sex.colors[meta[parents, "inferred_sex"]])
    }
  }
}

# After looking for parents, subset coverage data to --subset-samples plus any relevant parents
if(!is.null(args$subset_samples)){
  keep.samples <- NULL
  if(!is.null(args$subset_samples)){
    keep.samples <- read.table(args$subset_samples, header=F)[, 1]
  }
  cov <- cov[, c("chrom", "start", "end",
                 union(keep.samples, samples.to.plot))]
}


# Set hard-coded plot parameters
track.hex <- 0.075
gene.color <- EAS.colors["dark2"]

# Prepare plot area
pdf(args$out, height=2.75, width=3.5)
xlims <- range(cov[, c("start", "end")])
ylims.rd <- c(min(c(1, cov[, samples.to.plot], na.rm=T)),
              max(c(3, cov[, samples.to.plot], na.rm=T)) + 0.1)
n.extra.tracks <- length(unlist(args[c("sv_interval", "gene_features")]))
track.breaks <- ylims.rd[2] + (0:n.extra.tracks * track.hex * diff(ylims))
ylims <- c(ylims.rd[1], max(track.breaks))
prep.plot.area(xlims=xlims, ylims=ylims, yaxs="r", parmar=c(0.1, 2.4, 2.75, 0.5))

# Add horizontal lines at integer copy states
abline(h=0:10, lty=5, col="gray75")

# Add background shading corresponding to population averages
shade.middle.quantile(cov, 0.025, 0.975, "gray92")
shade.middle.quantile(cov, 0.25, 0.75, "gray85")

# Add lines for each highlight sample colored by cancer type
sapply(samples.to.plot, function(sid){
  points(step.function(cov$start, cov[, sid]), type="l", lwd=2,
         col=cancer.colors[metadata.cancer.label.map[meta[sid, "disease"]]])
})

# Prep top tracks
rect(xleft=-10e10, xright=10e10, ybottom=ylims.rd[2], ytop=par("usr")[4],
     border=NA, col="white", xpd=T)
# abline(h=ylims.rd[2], col="gray90")
track.content.buffer <- 0.01 * diff(par("usr")[3:4])

# Add SV call track
track.count <- 0
if(!is.null(args$sv_interval)){
  track.count <- track.count + 1
  sv.coords <- as.numeric(unlist(strsplit(unlist(strsplit(args$sv_interval, split=":"))[2], split="-")))
  rect(xleft=sv.coords[1], xright=sv.coords[2],
       ybottom=track.breaks[track.count]+track.content.buffer,
       ytop=track.breaks[track.count+1]-track.content.buffer,
       col="black")
  text(x=mean(sv.coords), y=mean(track.breaks[track.count + 0:1]),
       labels=args$sv_label, cex=4/6, col="white")
}

# Add gene body track
if(!is.null(args$gene_features)){
  track.count <- track.count + 1

  # Load gene features
  gdf <- read.table(args$gene_features, sep="\t", header=F)
  colnames(gdf) <- c("chrom", "start", "end", "feature")

  # Plot gene features
  tx.or.gene <- gdf[which(gdf$feature %in% c("transcript", "gene")), ]
  segments(x0=tx.or.gene$start, x1=tx.or.gene$end,
           y0=mean(track.breaks[track.count + 0:1]),
           y1=mean(track.breaks[track.count + 0:1]),
           lwd=2, col=gene.color, lend="butt")
  exons <- gdf[which(gdf$feature %in% c("exon", "UTR", "CDS")), ]
  rect(xleft=exons$start, xright=exons$end,
       ybottom=track.breaks[track.count + 1] - track.content.buffer,
       ytop=track.breaks[track.count] + track.content.buffer,
       border=gene.color, col=gene.color)

  # Add gene track label
  if(is.null(args$gene_label)){
    gene.label.font <- 1
    gene.label <- "Genes"
  }else{
    gene.label.font <- 3
    gene.label <- args$gene_label
  }
  if(any(gdf$start < xlims[1]) | is.null(args$gene_label)){
    axis(2, at=mean(track.breaks[track.count + 0:1]), labels=gene.label, tick=F,
         line=-1, cex=5/6, font=gene.label.font, col.axis=gene.color, las=2)
  }else{
    text(x=min(gdf$start), y=mean(track.breaks[track.count + 0:1]), pos=2,
         labels=gene.label, font=gene.label.font, col=gene.color, xpd=T)
  }
}

# Add X axis
if(diff(xlims) <= 100000){
  x.labels <- paste(round(axTicks(3) / 1000, 2), "kb", sep="")
}else{
  x.labels <- paste(round(axTicks(3) / 1000000, 2), "Mb", sep="")
}
clean.axis(3, at=axTicks(3), infinite=TRUE, labels=x.labels, title.line=0.25)

# Add Y axis
y.ax.at <- 0:floor(ylims.rd[2])
clean.axis(2, at=c(0, ylims.rd[2]), tck=0, labels=NA)
clean.axis(2, at=y.ax.at)
mtext("Copy Number", side=2, at=mean(ylims.rd), line=1.25)

# Add idiogram
clean.axis(3, at=xlims[1] + (diff(par("usr")[1:2]) * xlims/contig.lengths[cov$chrom[1]]),
           line=2, tck=-0.01, labels=NA, infinite=T)
clean.axis(3, infinite=T, labels=NA, line=0)
axis(3, line=0.6, at=par("usr")[1]-(0.025*diff(par("usr")[1:2])),
     labels=cov$chrom[1], xpd=T, tick=F, hadj=1, cex.axis=5/6)

# Close pdf
dev.off()

