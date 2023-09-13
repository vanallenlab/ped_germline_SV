#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Visualize global burden results from precomputed summary statistics
# See global_burden_tests.R for the generation of these statistics


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(PedSV, quietly=TRUE)
require(beeswarm, quietly=TRUE)
PedSV::load.constants("all")


##################
# Data functions #
##################
# Clean a single set of summary statistics for a single hypothesis for all four cancers
clean.sumstats.singleHyp <- function(ss, hypothesis, stat.prefix=NULL){
  tdf <- ss[which(ss$hypothesis == hypothesis),
            c("disease", "coefficient", "std.err", "P.value")]
  tdf$ci_diff <- (qnorm(0.975) * tdf$std.err)
  tdf$std.err <- NULL
  if(is.null(stat.prefix)){
    colnames(tdf) <- c("disease", "or", "p", "ci_diff")
  }else{
    colnames(tdf) <- c("disease", paste(stat.prefix, c("or", "p", "ci_diff"), sep="."))
  }
  return(tdf)
}

# Alias for helper function to column-wise join rare, very rare, and singleton stats
merge.by.freq <- function(ss.list){
  sdf <- Reduce(f=function(df1, df2){merge(df1, df2, by="disease")}, x=ss.list)
  rownames(sdf) <- sdf$disease
  sdf$disease <- NULL
  sdf[rev(setdiff(names(cancer.colors), "control")), ]
}


######################
# Plotting functions #
######################
# Plot effect sizes versus frequency bins for all four cancers
plot.or.by.freq <- function(plot.stats, title=NULL, connect.cancers=c(),
                            shaded.pancan.ci=FALSE, max.log.y=4){
  # Get plot values
  or.range <- range(plot.stats[, grep("\\.or$", colnames(plot.stats))], na.rm=T)
  ymin <- max(c(-max.log.y, min(c(0, or.range[1]-(diff(or.range)/4)))))
  ymax <- min(c(max.log.y, max(c(0, or.range[2]+(diff(or.range)/4)))))

  # Prep plot area
  prep.plot.area(xlims=c(0.25, 2.75), ylims=c(ymin, ymax), parmar=c(2.1, 2.6, 1.1, 0.1))
  axis(1, at=0.5:2.5, tick=F, line=-0.9, labels=c("AF<1%", "AF<0.1%", "AC=1"),
       cex.axis=5/6)
  mtext("SV Frequency", 1, 1)
  clean.axis(2, at=log(2^(-10:10)), labels=2^(-10:10), title="Odds Ratio", infinite=T)
  mtext(title, 3, 0)

  # Infer location of points when swarmed
  coords <- do.call("cbind", lapply(1:3, function(i){
    freq <- names(freq.names)[i]
    coords <- beeswarm(plot.stats[, paste(freq, "or", sep=".")], do.plot=F,
                       at=i-0.5, pch=23, corral="wrap", corralWidth=0.4)
    colnames(coords) <- paste(freq, colnames(coords), sep=".")
    coords[, 1:2]
  }))
  rownames(coords) <- rownames(plot.stats)

  # Add shaded area for pan-cancer CI
  if(shaded.pancan.ci){
    ci.x <- as.numeric(coords["pancan", grep("\\.x$", colnames(coords))])
    ci.y <- as.numeric(coords["pancan", grep("\\.y$", colnames(coords))])
    ci.y.all <- c(ci.y + plot.stats["pancan", grep("\\.ci_diff$", colnames(plot.stats))],
                  rev(ci.y - plot.stats["pancan", grep("\\.ci_diff$", colnames(plot.stats))]))
    polygon(x=c(ci.x, rev(ci.x)), y=ci.y.all,
            border=NA, bty="n", col="white")
    polygon(x=c(ci.x, rev(ci.x)), y=ci.y.all,
            border=NA, bty="n", col=adjustcolor(cancer.colors["pancan"], alpha=0.15))
  }

  # Add line for null (OR=1)
  abline(h=0, lty=5)

  # Add lines for individual cancers
    for(cancer in connect.cancers){
      segments(x0=as.numeric(coords[cancer, grep("\\.x$", colnames(coords))][1:2]),
               x1=as.numeric(coords[cancer, grep("\\.x$", colnames(coords))][2:3]),
               y0=as.numeric(coords[cancer, grep("\\.y$", colnames(coords))][1:2]),
               y1=as.numeric(coords[cancer, grep("\\.y$", colnames(coords))][2:3]),
               lwd=if(cancer=="pancan"){3}else{2}, col=cancer.colors[cancer])
  }

  # Add all point estimates last
  sapply(names(freq.names), function(freq){
    points(coords[, paste(freq, c("x", "y"), sep=".")], pch=23,
           bg=cancer.colors[rownames(plot.stats)], col="black")
  })
}


#################
# Rscript block #
#################
# Two simple command-line arguments
args <- commandArgs(trailingOnly=T)

# #DEV:
# args <- c("~/scratch/PedSV.v2.1.case_control.dev.global_burden_tests.tsv",
#           "~/scratch/PedSVv.2.1.case_control.dev")

# Load sumstats
ss <- read.table(args[1], header=T, sep="\t", comment.char="", check.names=F)
colnames(ss)[1] <- gsub("^#", "", colnames(ss)[1])

# Plot effect sizes of large SVs
sapply(c("", "DEL", "DUP", "INV", "CPX"), function(suffix){
  plot.stats <- merge.by.freq(lapply(names(freq.names), function(freq){
    hyp <- gsub("\\.$", "", paste("large", freq, suffix, sep="."))
    clean.sumstats.singleHyp(ss, hyp, freq)
  }))
  svtype <- if(suffix == ""){"SV"}else{suffix}
  sv.name <- if(suffix == ""){"SVs"}else{sv.abbreviations[suffix]}
  pdf(paste(args[2], ".large_", svtype, ".or_by_freq.pdf", sep=""),
      height=2.25, width=2.25)
  plot.or.by.freq(plot.stats, title=paste("Large (>1Mb)", sv.name),
                  shaded.pancan.ci=T, connect.cancers="pancan")
  dev.off()
})

# Plot effect sizes of gene-disruptive SVs
sapply(c("", "LoF", "CG", "IED", "nonLoF_disruptive"), function(suffix){
  plot.stats <- merge.by.freq(lapply(names(freq.names), function(freq){
    if(suffix == ""){
      hyp <- gsub("\\.$", "", paste(freq, "gene_disruptive", sep="."))
    }else{
      hyp <- gsub("\\.$", "", paste(freq, "gene_disruptive", paste(freq, suffix, "SVs", sep="_"), sep="."))
    }
    clean.sumstats.singleHyp(ss, hyp, freq)
  }))
  out.tag <- if(suffix == ""){"all"}else{suffix}
  title <- if(suffix == ""){"Gene-disruptive SVs"}else if(suffix == "nonLoF_disruptive"){"Non-LoF Genic SVs"}else{paste(suffix, "SVs")}
  pdf(paste(args[2], ".gene_disruptive_", out.tag, ".or_by_freq.pdf", sep=""),
      height=2.25, width=2.25)
  plot.or.by.freq(plot.stats, title=title, shaded.pancan.ci=T, connect.cancers="pancan")
  dev.off()
})

