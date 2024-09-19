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
require(beeswarm, quietly=TRUE)
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
  x.step <- c(x[1], RLCtools::stretch.vector(x, 2))
  y.step <- c(y[1], y[1], RLCtools::stretch.vector(y, 2))[1:length(x.step)]
  data.frame("x" = x.step, y=y.step)
}

# Add shaded stepwise normal population variability in coverage
shade.middle.quantile <- function(cov, lower, upper, color, method="quantile"){
  if(method == "quantile"){
    plot.vals <- step.function(c(cov$start, rev(cov$start)),
                               c(apply(cov[, -c(1:3)], 1, quantile, prob=lower, na.rm=T),
                                 apply(cov[, -c(1:3)], 1, quantile, prob=upper, na.rm=T)))
  }else if(method == "zscore"){
    plot.vals <- step.function(c(cov$start, rev(cov$start)),
                               c(apply(cov[, -c(1:3)], 1, function(v){mean(v) + (lower * sd(v))}),
                               rev(apply(cov[, -c(1:3)], 1, function(v){mean(v) + (upper * sd(v))}))))
  }
  polygon(plot.vals, border=NA, col=color)
}

# Add gene features at a prespecified location
add.gene.features <- function(gdf, y.mid, y.lims, gene.color, alpha=1){
  tx.or.gene <- gdf[which(gdf$feature %in% c("transcript", "gene")), ]
  segments(x0=tx.or.gene$start, x1=tx.or.gene$end,
           y0=y.mid, y1=y.mid, lwd=2, lend="butt",
           col=adjustcolor(gene.color, alpha=alpha))
  exons <- gdf[which(gdf$feature %in% c("exon", "UTR", "CDS")), ]
  rect(xleft=exons$start, xright=exons$end,
       ybottom=min(y.lims), ytop=max(y.lims),
       border=adjustcolor(gene.color, alpha=alpha),
       col=adjustcolor(gene.color, alpha=alpha))
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
                    action="append",
                    help="IDs for samples to highlight on panel.")
parser$add_argument("--subset-samples", metavar=".tsv", type="character",
                    help="List of samples to subset [default: use all samples]")
parser$add_argument("--sv-interval", metavar="chrom:start-end", type="character",
                    help="SV coordinates to highlight")
parser$add_argument("--sv-label", metavar="string", type="character", default="SV",
                    help="Label for --sv-interval [default: 'SV']")
parser$add_argument("--highlight-gene-features", metavar=".bed", type="character",
                    help="BED file of gene features to plot prominently.")
parser$add_argument("--highlight-gene-label", metavar="string", type="character",
                    help="Label for --highlight-gene-features [default: 'Gene']")
parser$add_argument("--background-gene-features", metavar=".bed", type="character",
                    help="BED file of gene features to plot subtly.")
parser$add_argument("--no-parents", action="store_true",
                    help="Disable plotting of parents for complete trios")
parser$add_argument("--no-idiogram", action="store_true",
                    help="Disable plotting of idiogram skeleton on top X axis")
parser$add_argument("--outfile", metavar="path", type="character", required=TRUE,
                    help="Path to output .pdf")
args <- parser$parse_args()

# # DEV
# args <- list("bincov" = "/Users/ryan/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_MS/PedSV_figures/other_data/All_20_Batches.chr2.final_cleanup_DUP_chr2_812.chr2_15046644_16839295.rd.w_PT_EJW8WH3G.bed.gz",
#              "cov_medians" = "~/scratch/PedSV.v2.5.3.median_coverage.tsv.gz",
#              "metadata" = "~/scratch/PedSV.v2.5.3.cohort_metadata.w_control_assignments.tsv.gz",
#              "sample_id" = "PT_V1Q9W1NW",
#              "subset_samples" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_v2_callset_generation/v2.5.3/PedSV.v2.5.3.full_cohort_w_relatives_1000G.samples.list",
#              "sv_interval" = "chr2:15546644-16339295",
#              "sv_label" = "793kb duplication",
#              "highlight_gene_features" = "~/scratch/MYCN.features.bed",
#              "highlight_gene_label" = "MYCN",
#              "background_gene_features" = NULL,
#              "no_parents" = FALSE,
#              "no_idiogram" = FALSE,
#              "outfile" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_MS/PedSV_figures/other_data/MYCN.de_novo.DUP.rdviz.pdf")
# args <- list("bincov" = "~/scratch/All_20_Batches.chr3.final_cleanup_DUP_chr3_397.chr3_12470403_12885616.rd.bed.gz",
#              "cov_medians" = "~/scratch/PedSV.v2.5.3.median_coverage.tsv.gz",
#              "metadata" = "~/scratch/PedSV.v2.5.3.cohort_metadata.w_control_assignments.tsv.gz",
#              "sample_id" = c("PT_5CZHK9CR", "PT_ZT2NW6WA", "SJ071781"),
#              "subset_samples" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_v2_callset_generation/v2.5.3/PedSV.v2.5.3.full_cohort_w_relatives_1000G.samples.list",
#              "sv_interval" = "chr3:12589035-12766984",
#              "sv_label" = "178kb duplication",
#              "highlight_gene_features" = "~/scratch/RAF1_TMEM40.features.bed.gz",
#              "highlight_gene_label" = "RAF1",
#              "background_gene_features" = NULL,
#              "no_parents" = TRUE,
#              "no_idiogram" = FALSE,
#              "outfile" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_MS/PedSV_figures/other_data/RAF1_TMEM40.DUP.rdviz.pdf")
# args <- list("bincov" = "~/scratch/YL_SV_v1.2_DUP_chr1_1376.chr1_144833951_146729681.rd.bed.gz",
#              "cov_medians" = "~/scratch/YL.SV.v1.2.median_coverage.tsv.gz",
#              "metadata" = "~/scratch/YL.SV.v1.2.analysis_metadata.tsv.gz",
#              "sample_id" = c("TPMCCDG12268", "NSLC-0233"),
#              "subset_samples" = "~/scratch/YL.SV.v1.2.analysis.samples.list",
#              "sv_interval" = "chr1:145375588-146188044",
#              "sv_label" = "812 kb duplication",
#              "highlight_gene_features" = NULL,
#              "highlight_gene_label" = NULL,
#              "background_gene_features" = NULL,
#              "no_parents" = TRUE,
#              "no_idiogram" = FALSE,
#              "outfile" = "~/scratch/YL.SV.dev.rdviz.pdf")


# Load sample metadata
meta <- load.sample.metadata(args$metadata)

# Read & normalize read depth data
cov <- load.coverage(args$bincov, args$cov_medians)

# Check if both parents are present; if so, prepend these to args$sample_id
samples.to.plot <- intersect(args$sample_id, colnames(cov))
sample.colors <- sapply(cancer.palettes[metadata.cancer.label.map[meta[samples.to.plot, "disease"]]],
                        function(p){p["dark1"]})
if(!args$no_parents & "family_id" %in% colnames(meta)){
  for(sid in args$sample_id){
    fam.id <- meta[sid, "family_id"]
    if(is.na(fam.id)){
      next
    }
    parents <- rownames(meta)[which(meta$family_id == fam.id & rownames(meta) != sid)]
    n.parents <- length(parents)
    if(n.parents == 2 & all(parents %in% colnames(cov))){
      child.pheno <- metadata.cancer.label.map[meta[sid, "disease"]]
      samples.to.plot <- c(parents, samples.to.plot)
      sample.colors <- c(rep(cancer.palettes[[child.pheno]]["light2"], 2), sample.colors)
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
                 intersect(colnames(cov), union(keep.samples, samples.to.plot)))]
}

# If region is on sex chromosome and all --sample-id are the same sex,
# subset only to samples with matching sex
if(cov[1, 1] %in% c("chrX", "chrY")){
  target.sexes <- table(meta[args$sample_id, "inferred_sex"])
  if(length(target.sexes) == 1){
    target.sex <- names(target.sexes)
    sex.matches <- rownames(meta)[which(meta$inferred_sex == target.sex)]
    cov <- cov[, c("chrom", "start", "end",
                   intersect(colnames(cov), union(sex.matches, samples.to.plot)))]
  }
}

# Get coverage median for each sample over the SV interval (used for margin swarmplot)
if(!is.null(args$sv_interval)){
  sv.coords <- as.numeric(unlist(strsplit(unlist(strsplit(args$sv_interval, split=":"))[2], split="-")))
  sample.pch <- c("MALE"=15, "FEMALE"=19)[meta[samples.to.plot, "inferred_sex"]]
  sv.cov.meds <- sapply(samples.to.plot, function(sid){
    median(cov[which(cov$start <= sv.coords[2] & cov$end >= sv.coords[1]), sid], na.rm=T)
  })
}

# Set hard-coded plot parameters
track.hex <- 0.11
gene.color <- EAS.colors["dark2"]
background.gene.color <- EAS.colors["light2"]
sample.lwds <- rep(1 + (0.8 / (length(samples.to.plot) ^ (1/3))),
                   times=length(samples.to.plot))

# Prepare plot area
pdf(args$out,
    height=2.75-(0.25*as.numeric(args$no_idiogram)),
    width=3.5-(0.25*as.numeric(args$no_idiogram)))
xlims <- range(cov[, c("start", "end")])
ylims.rd <- c(min(c(1, unlist(cov[, samples.to.plot]), na.rm=T)) - 0.1,
              max(c(3, unlist(cov[, samples.to.plot]), na.rm=T)) + 0.1)
n.extra.tracks <- length(args$sv_interval) + min(c(1, length(unlist(args[c("highlight_gene_features", "background_gene_features")]))))
track.breaks <- ylims.rd[2] + (0:n.extra.tracks * track.hex * diff(ylims.rd))
ylims <- c(ylims.rd[1], max(track.breaks) + (0.015*diff(ylims.rd)))
prep.plot.area(xlims=xlims, ylims=ylims,
               parmar=if(args$no_idiogram){c(0.35, 3.3, 2.25, 0.5)}else{c(0.35, 3.3, 2.75, 1.2)})
track.content.buffer <- 0.07 * diff(par("usr")[3:4])

# Add background shading corresponding to population averages
shade.middle.quantile(cov, 0.025, 0.975, "gray90")
# shade.middle.quantile(cov, 0.25, 0.75, "gray85")
# cov.z <- seq(qnorm(0.995), 0, length.out=20)
# for(i in 1:length(cov.z)){
#   shade.middle.quantile(cov, -cov.z[i], cov.z[i],
#                         color=paste("gray", 98-i, sep=""), method="zscore")
# }
# cov.q <- seq(0.5-0.025, 0, length.out=20)
# for(i in 1:length(cov.q)){
#   shade.middle.quantile(cov, 0.5-cov.q[i], 0.5+cov.q[i],
#                         color=paste("gray", 96-i, sep=""))
# }

# Add horizontal lines at integer copy states
segments(x0=xlims[1], x1=xlims[2], y0=0:max(ylims.rd), y1=0:max(ylims.rd),
         lty=5, col="gray75", xpd=T)

# Add lines for each highlight sample colored by cancer type
sapply(1:length(samples.to.plot), function(i){
  sid <- samples.to.plot[i]
  points(step.function(cov$start, cov[, sid]), type="l",
         lwd=sample.lwds[i], col=sample.colors[i])
})

# Prep top tracks
rect(xleft=-10e10, xright=10e10, ybottom=ylims.rd[2], ytop=par("usr")[4],
     border=NA, col="white", xpd=T)
# abline(h=ylims.rd[2], col="gray90")

# Add SV call track
track.count <- 0
if(!is.null(args$sv_interval)){
  track.count <- track.count + 1
  # rect(xleft=sv.coords[1], xright=sv.coords[2],
  #      ybottom=track.breaks[track.count]+track.content.buffer,
  #      ytop=track.breaks[track.count+1]-track.content.buffer,
  #      col="black")
  segments(x0=sv.coords[1], x1=sv.coords[2],
           y0=track.breaks[track.count], y1=track.breaks[track.count],
           lwd=2)
  text(x=mean(sv.coords), y=mean(track.breaks[track.count + 0:1]),
       labels=args$sv_label, cex=4.5/6)
}

# Add gene body track
if(!is.null(args$highlight_gene_features)){
  track.count <- track.count + 1

  # Load & plot background gene features
  if(!is.null(args$background_gene_features)){
    b.gdf <- read.table(args$background_gene_features, sep="\t", header=F)
    colnames(b.gdf) <- c("chrom", "start", "end", "feature")
    add.gene.features(b.gdf, y.mid=mean(track.breaks[track.count + 0:1]),
                      y.lims=c(track.breaks[track.count + 1] - track.content.buffer,
                               track.breaks[track.count] + track.content.buffer),
                      background.gene.color, alpha=0.5)
  }

  # Load & plot highlight gene features
  gdf <- read.table(args$highlight_gene_features, sep="\t", header=F)
  colnames(gdf) <- c("chrom", "start", "end", "feature")
  add.gene.features(gdf, y.mid=mean(track.breaks[track.count + 0:1]),
                    y.lims=c(track.breaks[track.count + 1] - track.content.buffer,
                             track.breaks[track.count] + track.content.buffer),
                    gene.color, alpha=1)

  # Add gene track label
  if(is.null(args$highlight_gene_label)){
    gene.label.font <- 1
    gene.label <- "Genes"
    gene.label.color <- background.gene.color
  }else{
    gene.label.font <- 4
    gene.label <- args$highlight_gene_label
    gene.label.color <- gene.color
  }
  if(any(gdf$start < xlims[1]) | is.null(args$highlight_gene_label)){
    axis(2, at=mean(track.breaks[track.count + 0:1]), labels=gene.label, tick=F,
         line=-0.9, cex=5/6, font=gene.label.font, col.axis=gene.label.color, las=2)
  }else{
    text(labels=gene.label, x=min(gdf$start),
         y=mean(track.breaks[track.count + 0:1])-(0.005*diff(par("usr")[3:4])),
         pos=2, xpd=T, font=gene.label.font, col=gene.label.color, cex=5/6)
  }
}

# Add X axis
x.labels <- unique(sort(round(axTicks(3) / 1000000, 3)))
max.x.ticks <- 5
if(length(x.labels) > max.x.ticks){
  x.labels <- x.labels[c(TRUE, FALSE)]
}
x.ax.at <- 1000000 * x.labels
if(args$no_idiogram){
  x.title <- paste(cov[1, 1], "coordinate (Mb)")
}else{
  x.title <- NULL
  axis(3, at=par("usr")[2]+(0.01*diff(par("usr")[1:2])), tick=F, line=-1.4,
       cex.axis=5/6, labels="Mb", hadj=0, xpd=T)
}
clean.axis(3, at=x.ax.at, infinite=TRUE, labels=x.labels, title=x.title, title.line=0.25)

# Add left margin swarmplot of median coverages per sample
if(!is.null(args$sv_interval)){
  points(beeswarm(sv.cov.meds, add=T, method="compactswarm", corral="random",
                  at=par("usr")[1]-(0.11*diff(par("usr")[1:2])),
                  corralWidth=diff(c(0.07, 0.15)*diff(par("usr")[1:2]))),
         xpd=T, pch=sample.pch, xpd=T, col=sample.colors,
         cex=1.1/(length(samples.to.plot)^(1/3)))
}

# Add Y axis
y.ax.at <- 0:floor(ylims.rd[2])
clean.axis(2, at=c(0, ylims.rd[2]), tck=0, labels=NA)
clean.axis(2, at=y.ax.at)
mtext("Copy number", side=2, at=mean(ylims.rd), line=1.9)

# Add idiogram
if(!args$no_idiogram){
  clean.axis(3, at=xlims[1] + (diff(par("usr")[1:2]) * xlims/contig.lengths[cov$chrom[1]]),
             line=2, tck=-0.01, labels=NA, infinite=T)
  clean.axis(3, infinite=T, labels=NA, line=0)
  axis(3, line=0.6, at=par("usr")[1]-(0.025*diff(par("usr")[1:2])),
       labels=cov$chrom[1], xpd=T, tick=F, hadj=1, cex.axis=5/6)
}

# Close pdf
dev.off()

