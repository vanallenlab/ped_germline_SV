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

# Declare local constants
svtypes <- c("DEL", "DUP", "CNV", "INS", "INV", "CPX", "CTX")
sv.dens.bw <- c("DEL" = 3/4,
                "DUP" = 1,
                "CNV" = 0.5,
                "INS" = 10,
                "INV" = 2/3,
                "CPX" = 1/3)


##################
# Data functions #
##################
# Collect a table of counts per SV type stratified by frequency bins
get.counts.table <- function(bed, af.field="AF", ac.field="AC"){
  sapply(rev(svtypes), function(svtype){
    if(svtype == "CNV"){
      rev(c(0, 0, length(which(bed$SVTYPE == svtype))))
    }else{
      rev(c(length(which(bed$SVTYPE == svtype & bed[, ac.field] <= 1)),
            length(which(bed$SVTYPE == svtype & bed[, ac.field] > 1 & bed[, af.field] < 0.01)),
            length(which(bed$SVTYPE == svtype & bed[, af.field] >= 0.01))))
    }
  })
}

# Collect a list of log10-scaled SV size density functions for plotting
get.svlen.densities <- function(bed){
  svtypes.sub <- rev(setdiff(svtypes, "CTX"))
  data <- lapply(svtypes.sub, function(svtype){
    density(log10(bed$SVLEN[bed$SVTYPE == svtype]), adjust=sv.dens.bw[svtype])
  })
  names(data) <- sv.abbreviations[svtypes.sub]
  return(data)
}

# Tabulate SV counts by type and frequency
tabulate.counts <- function(bed){
  counts <- as.data.frame(t(get.counts.table(bed)))
  colnames(counts) <- c("common_af_ge_1pct", "rare_af_lt_1pct", "singleton_ac_1")
  counts$ALL <- apply(counts, 1, sum)
  counts$rare_af_lt_1pct <- counts$rare_af_lt_1pct + counts$singleton_ac_1
  counts["ALL", ] <- apply(counts, 2, sum)
  res <- do.call("rbind", lapply(rev(rownames(counts)), function(svtype){
    do.call("rbind", lapply(rev(colnames(counts)), function(freq){
      c(svtype, freq, counts[svtype, freq])
    }))
  }))
  res <- as.data.frame(res)
  colnames(res) <- c("#svtype", "freq", "N")
  return(res[which(res$N > 0), ])
}

# Tabulate median SV sizes by type and frequency
tabulate.sizes <- function(bed, af.field="AF", ac.field="AC"){
  freq.idxs <- list("ALL" = 1:nrow(bed),
                    "singleton_ac_1" = which(bed[, ac.field] <= 1),
                    "rare_af_lt_1pct" = which(bed[, af.field] < 0.01),
                    "common_af_ge_1pct" = which(bed[, af.field] >= 0.01))
  res <- lapply(c("ALL", setdiff(svtypes, "CTX")), function(svtype){
    sv.idxs <- if(svtype == "ALL"){1:nrow(bed)}else{which(bed$SVTYPE == svtype)}
    do.call("rbind", lapply(names(freq.idxs), function(freq){
      hit.idxs <- intersect(freq.idxs[[freq]], sv.idxs)
      if(length(hit.idxs) > 0){c(svtype, freq, median(bed$SVLEN[hit.idxs]))}
    }))
  })
  df <- as.data.frame(do.call("rbind", res))
  colnames(df) <- c("#svtype", "freq", "N")
  return(df)
}


######################
# Plotting functions #
######################
# Plot a horizontal stacked barplot of SVs colored by frequency bin
plot.count.bars <- function(bed, af.field="AF", ac.field="AC", greyscale=TRUE,
                            bar.buffer=0.15, count.label.xadj=0.075){
  # Get counts matrix
  counts <- get.counts.table(bed, af.field, ac.field)

  # Get plot dimensions
  xlims <- c(0, max(apply(counts, 2, sum)))
  ylims <- c(-2, ncol(counts))
  label.xadj <- count.label.xadj * diff(xlims)

  # Prep plot area
  PedSV::prep.plot.area(xlims, ylims, parmar=c(0.25, 3.5, 0.1, 2.5), xaxs="i", yaxs="i")
  segments(x0=xlims[1], x1=xlims[1], y0=0, y1=ylims[2], col="gray85", xpd=T)

  # Add bars & cap labels
  sapply(1:ncol(counts), function(i){
    x.stops <- c(0, cumsum(counts[, i]))
    bar.pal <- sv.palettes[[colnames(counts)[i]]]
    if(greyscale){
      bar.pal <- hex2grey(bar.pal)
    }
    rect(xleft=x.stops[-length(x.stops)], xright=x.stops[-c(1)],
         ybottom=i-1+bar.buffer, ytop=i-bar.buffer,
         col=bar.pal[c("dark2", "main", "light2")],
         border=bar.pal[c("dark2", "main", "light2")])
    n.total <- x.stops[length(x.stops)]
    rect(xleft=0, xright=n.total, ybottom=i-1+bar.buffer, ytop=i-bar.buffer,
         col=NA, xpd=T)
    if(n.total > 10000){
      count.label <- paste(round(n.total/1000, 1), "k", sep="")
    }else{
      count.label <- prettyNum(n.total, big.mark=",")
    }
    text(x=x.stops[length(x.stops)]-label.xadj, y=i-0.5, labels=count.label,
         pos=4, cex=5/6, xpd=T)
  })

  # Add left axis
  axis(2, at=(1:ncol(counts))-0.5, las=2, tick=F, line=-0.8,
       labels=sv.abbreviations[colnames(counts)])

  # Add legend
  legend.x.at <- c(-0.4, 0.85, -0.4)*diff(xlims)
  legend.y.at <- c(-0.5, -0.5, -1.5)-0.2
  points(x=legend.x.at, y=legend.y.at, xpd=T, pch=15, cex=1.3,
         col=hex2grey(DEL.colors[c("dark2", "main", "light2")]))
  text(x=legend.x.at, y=legend.y.at, pos=4, xpd=T,
       labels=c("Common", "Rare", "Singleton"))
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Plot SV counts and sizes")
parser$add_argument("bed", metavar=".tsv", type="character",
                    help="SV sites .bed")
parser$add_argument("--af-field", default="AF", metavar="string", type="character",
                    help="column header to use for AF-related analyses")
parser$add_argument("--ac-field", default="AC", metavar="string", type="character",
                    help="column header to use for AF-related analyses")
parser$add_argument("--out-prefix", metavar="path", type="character", required=TRUE,
                    help="path/prefix for all output files")
args <- parser$parse_args()

# # DEV:
# args <- list("bed" = "~/scratch/PedSV.v1.1.validation_cohort.analysis_samples.wAFs.bed.gz",
#              "af_field" = "AF",
#              "ac_field" = "AC",
#              "out_prefix" = "~/scratch/PedSV.dev")

# Load BED
bed <- PedSV::load.sv.bed(args$bed)

# Plot stacked bar colored by frequency and type
pdf(paste(args$out_prefix, "sv_site_counts.pdf", sep="."),
    height=1.7, width=2)
plot.count.bars(bed, args$af_field, args$ac_field)
dev.off()

# Ridgeplot of SV sizes
pdf(paste(args$out_prefix, "sv_size_distribs.pdf", sep="."),
    height=1.85, width=2.2)
ridgeplot(get.svlen.densities(bed), xlims=log10(c(10, 5000000)), x.axis=FALSE,
          fill=hex2grey(DEL.colors[["light2"]]),
          border=hex2grey(DEL.colors[["dark1"]]), border.lwd=1.25,
          parmar=c(2.2, 3.5, 0.1, 0.1))
clean.axis(1, at=log10(logscale.major.bp),
           labels=logscale.major.bp.labels[seq(1, length(logscale.major.bp), 2)],
           labels.at=log10(logscale.major.bp)[seq(1, length(logscale.major.bp), 2)],
           label.line=-0.9, title.line=0.2, title=bquote("SV Size" ~ (log[10])))
dev.off()

# Table of SV counts by type and frequency
write.table(tabulate.counts(bed),
            paste(args$out_prefix, "sv_counts_by_freq.tsv", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)

# Table of median SV sizes by type and frequency
write.table(tabulate.sizes(bed),
            paste(args$out_prefix, "sv_sizes_by_freq.tsv", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)
