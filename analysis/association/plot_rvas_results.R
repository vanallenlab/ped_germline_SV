#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Visualize summary statistics from gene-based rare SV burden tests


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(PedSV, quietly=TRUE)
require(argparse, quietly=TRUE)
PedSV::load.constants("all")

csq.labels <- c("LOF" = "Loss-of-Function",
                "GeneDisruptive" = "Gene-Disruptive",
                "CG" = "Copy-Gain")


# Helper function to plot QQ
# Eventually this should be moved out to PedSVR library
plot.qq <- function(p, smallest.p=NULL, cutoff=NULL, highlights=NULL,
                    cutoff.color="red",
         highlight.color="#4EE69A",
         highlight.name="Positive Controls",
         print.stats=T, echo.lambdas=F,
         legend=T, xmax=NULL, ymax=NULL, reflection=F,
         label.cex=1, pt.color="grey15", pt.cex=0.2,
         parmar=c(2.1, 3.1, 0.6, 1)){

  # Dummy function to plot N/A QQ for analyses with no cases
  empty.qq <- function(){
    par(mar=parmar, bty="n")
    plot(x=0:1, y=0:1, type="n", xaxt="n", yaxt="n",
         xlab="", ylab="", xaxs="i", yaxs="i")
    axis(1, at=axTicks(1), labels=NA, tck=-0.02)
    mtext(1, text=expression(Expected ~ ~-log[10](italic(p))),
          line=1.2, cex=label.cex)
    axis(2, at=axTicks(2), labels=NA, tck=-0.02)
    mtext(2, text=expression(Observed ~ ~-log[10](italic(p))),
          line=1.5, cex=label.cex)
    text(x=0.5,y=0.5,labels="N/A")
  }

  # Main QQ function
  if(all(p==1 | is.na(p) | is.nan(p) | is.infinite(p))){
    empty.qq()
  }else{
    colors <- rep(pt.color, times=length(p))
    if(!is.null(highlights)){
      hits <- sort(unique(unlist(apply(highlights, 1, function(coords){
        which(stats$chr == as.character(coords[1])
              & stats$pos <= as.numeric(coords[3])
              & stats$pos >= as.numeric(coords[2]))
      }))))
      colors[hits] <- highlight.color
    }

    if (!is.numeric(p)){
      stop("P values must be numeric.")
    }
    keep.idxs <- which(!is.na(p) & !is.nan(p) & !is.null(p) &
                         is.finite(p) & p <= 1 & p >= 0)
    p <- p[keep.idxs]
    colors <- colors[keep.idxs]

    new.order <- order(p)
    p <- p[new.order]
    colors <- colors[new.order]
    hits <- which(colors == highlight.color)

    expected <- ppoints(length(p))
    qqconf <- function (p.expected, reflection=F){
      n <- length(p.expected)
      mpts <- matrix(nrow=n * 2, ncol=2)
      for (i in seq(from=1, to=n)) {
        mpts[i, 1] <- -log10((i - 0.5)/n)
        mpts[i, 2] <- -log10(qbeta(0.975, i, n - i))
        mpts[n * 2 + 1 - i, 1] <- -log10((i - 0.5)/n)
        mpts[n * 2 + 1 - i, 2] <- -log10(qbeta(0.025, i, n -
                                                 i))
      }
      mpts <- as.data.frame(mpts)
      if(reflection == T){
        mpts[, 2] <- -mpts[, 2]
      }
      return(mpts)
    }
    conf.int <- qqconf(expected, reflection)

    lambda <- dchisq(median(p), df=1)/dchisq(median(expected), df=1)
    if(!is.null(highlights)){
      expected.hits <- ppoints(length(p[hits]))
      lambda.hits <- dchisq(median(p[hits]), df=1)/dchisq(median(expected.hits), df=1)
      lambda.ratio <- lambda.hits / lambda
      if(echo.lambdas==T){
        cat(paste(round(lambda, 4), "\t",
                  round(lambda.hits, 4), "\t",
                  round(lambda.ratio, 4), "\n",
                  sep=""))
      }
    }

    if(is.null(cutoff)){
      cutoff <- 0.05/length(p)
    }

    p.rounded <- FALSE
    if(!is.null(smallest.p)){
      if(any(p < smallest.p)){
        p[which(p < smallest.p)] <- smallest.p
        p.rounded <- TRUE
      }
    }

    if(is.null(ymax)){
      maxp <- max(-log10(p[which(!(is.infinite(-log10(p))))]))
      ymax <- max(c(maxp, -log10(cutoff) + 2))
    }

    expected <- -log10(expected)
    if(is.null(xmax)){
      xmax <- max(expected, na.rm=T)
    }

    if(reflection == F){
      p <- -log10(p)
      log.cutoff <- -log10(cutoff)
      ab.end <- 1
      x.ax <- 1
      mars <- c(2.1, 3.1, 0.5, 1)
      lambda.ypos <- c(0.9, 0.825, 0.75)
      x.title <- T
    }else{
      p <- log10(p)
      log.cutoff <- log10(cutoff)
      ab.end <- -1
      ymax <- -ymax
      x.ax <- 3
      mars <- c(0.5, 3.1, 0.6, 1)
      lambda.ypos <- rev(c(0.9, 0.825, 0.75))
      x.title <- F
    }

    par(mar=parmar, bty="n")
    plot(x=expected, y=p, type="n", xaxt="n", yaxt="n",
         xlab="", ylab="", xaxs="i", yaxs="i",
         xlim=c(0, 1.1 * xmax), ylim=range(c(0, 1.1 * ymax)))
    polygon(x=conf.int[, 1], y=conf.int[, 2], col="gray90", border=NA)
    abline(0, ab.end, col="gray50")
    if(!is.null(cutoff)){
      abline(h=log.cutoff, col=cutoff.color, lty=2)
    }

    if (print.stats == T){
      xpos.adjust <- 0.025
      stat.label.color <- "gray20"
      if(reflection == T){
        ypos.ref <- par("usr")[3]
      }else{
        ypos.ref <- par("usr")[4]
      }
      if(is.null(highlights)){
        text(par("usr")[1] - (xpos.adjust * (par("usr")[2] - par("usr")[1])),
             lambda.ypos[1] * ypos.ref, pos=4, font=2, cex=0.8, col=stat.label.color,
             labels=bquote(lambda ~ .(paste("=", sprintf("%.2f", lambda), sep=""))))
      }else{
        text(par("usr")[1] - (xpos.adjust * (par("usr")[2] - par("usr")[1])),
             lambda.ypos[1] * ypos.ref, pos=4, font=2, cex=0.8, col=stat.label.color,
             labels=bquote(lambda[italic("All")] ~ .(paste("=", sprintf("%.2f", lambda), sep=""))))
        text(par("usr")[1] - (xpos.adjust * (par("usr")[2] - par("usr")[1])),
             lambda.ypos[2] * ypos.ref, pos=4, font=2, cex=0.8, col=stat.label.color,
             labels=bquote(lambda[italic("Cyan")] ~ .(paste("=", sprintf("%.2f", lambda.hits), sep=""))))
        text(par("usr")[1] - (xpos.adjust * (par("usr")[2] - par("usr")[1])),
             lambda.ypos[3] * ypos.ref, pos=4, font=2, cex=0.8, col=stat.label.color,
             labels=bquote(lambda[italic("C")] / lambda[italic("A")] ~ .(paste("=", sprintf("%.2f", lambda.ratio), sep=""))))
      }
    }

    points(x=expected, y=p, pch=19, col=colors, cex=pt.cex)

    if(!is.null(highlights)){
      if(any(colors == highlight.color)){
        points(x=expected[which(colors == highlight.color)],
               y=p[which(colors == highlight.color)],
               col=colors[which(colors == highlight.color)],
               pch=19, cex=pt.cex)
      }
      if(legend==T){
        legend("left", pch=19, col=highlight.color, cex=0.5,
               legend=paste(highlight.name, " (n=",
                            prettyNum(nrow(highlights), big.mark=","),
                            ")", sep=""))
      }
    }

    axis(x.ax, at=c(0, 10e10), labels=NA, tck=0)
    axis(x.ax, at=axTicks(1), tick=F, line=-1.1, cex.axis=0.75, labels=abs(axTicks(1)))
    if(x.title == T){
      mtext(x.ax, text=expression(Expected ~ ~-log[10](italic(p))),
            line=1.15, cex=label.cex)
    }

    if(reflection == T){
      y.at <- seq(0, floor(par("usr")[3]), by=floor(par("usr")[3]/6))
    }else{
      y.at <- seq(0, ceiling(par("usr")[4]), by=ceiling(par("usr")[4]/6))
    }
    if(p.rounded){
      y.at <- y.at[-length(y.at)]
      axis(2, at=-log10(smallest.p), labels=NA, tck=-0.02)
      axis(2, at=-log10(smallest.p), tick=F, las=2, line=-0.6, cex.axis=0.75,
           labels=bquote("" >= .(-log10(smallest.p))))
    }
    axis(2, at=c(0, 10e10), labels=NA, tck=0)
    axis(2, at=y.at, labels=NA, tck=-0.02)
    axis(2, at=y.at, cex.axis=0.75, tick=F, line=-0.6, las=2, labels=abs(y.at))
    mtext(2, text=expression(Observed ~ ~-log[10](italic(p))),
          line=1.5, cex=label.cex)
  }
}




###########
# RScript #
###########
# TODO: update this to be scriptable
# # Parse command line arguments and options
# parser <- ArgumentParser(description="Case:control SV burden tests")
# parser$add_argument("--bed", metavar=".tsv", type="character", action="append",
#                     help=paste("SV sites .bed. Can be supplied one or more",
#                                "times. Must be supplied at least once. Order",
#                                "must match --ad."), required=TRUE)
# parser$add_argument("--ad", metavar=".tsv", type="character", action="append",
#                     help=paste("Allele dosage .bed. Can be supplied one or more",
#                                "times. Must be supplied at least once. Order must",
#                                "match --bed."), required=TRUE)
# parser$add_argument("--metadata", metavar=".tsv", type="character",
#                     help="sample metadata .tsv", required=TRUE)
# parser$add_argument("--gene-lists", metavar=".tsv", type="character",
#                     help="Two-column .tsv of gene lists (set name and path to gene list)")
# parser$add_argument("--genomic-disorder-hits", metavar=".txt", type="character",
#                     help=".txt of variant IDs to be treated as genomic disorder CNVs")
# parser$add_argument("--subset-samples", metavar=".tsv", type="character",
#                     help="list of samples to subset [default: use all samples]")
# parser$add_argument("--exclude-variants", metavar=".txt", type="character",
#                     help=paste("list of variant IDs to exclude from analysis.",
#                                "Not recommended for most use-cases.",
#                                "[default: use all variants]"))
# parser$add_argument("--af-field", metavar="string", type="character", action="append",
#                     help=paste("Column header to use for AF-related analyses.",
#                                "Can be supplied once to be applied uniformly to",
#                                "all values of --bed, or can be supplied an",
#                                "equal number of times as --bed to apply one",
#                                "value per --bed. [default: AF]"))
# parser$add_argument("--ac-field", metavar="string", type="character", action="append",
#                     help=paste("Column header to use for AC-related analyses.",
#                                "Same requirements as --af-field. [default: AC]"))
# parser$add_argument("--out-prefix", metavar="path", type="character", required=TRUE,
#                     help="path/prefix for all output files")
# args <- parser$parse_args()

# DEV: load local sumstats for development
ss <- read.table("~/scratch/PedSV.v2.5.2.sv_rvas.sumstats.bed.gz",
                 header=T, sep="\t", comment.char="")

# Iterate over cancers
sapply(setdiff(names(cancer.colors), "control"), function(cancer){
  # Iterate over consequences
  sapply(unique(ss$consequence), function(csq){
    pvals <- 10^-as.numeric(ss[which(ss$consequence == csq),
                           paste(cancer, "neglog10_p", sep=".")])
    png(paste("~/scratch/", cancer, ".", csq, ".qq.png", sep=""),
        height=300*2.75, width=300*2.75, res=300)
    plot.qq(pvals, cutoff=0.05 / length(pvals), pt.color=cancer.colors[cancer],
            parmar=c(2.1, 3.1, 2, 1))
    mtext(paste(cancer.names.long[cancer], " vs. Control\n",
          "Rare ", csq.labels[csq], " SVs"), 3)
    dev.off()
  })
})


