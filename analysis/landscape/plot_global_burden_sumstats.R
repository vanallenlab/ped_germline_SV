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
  tdf$lower <- exp(tdf$coefficient + (qnorm(0.025) * tdf$std.err))
  tdf$upper <- exp(tdf$coefficient + (qnorm(0.975) * tdf$std.err))
  tdf$coefficient <- exp(tdf$coefficient)
  tdf$std.err <- NULL
  if(is.null(stat.prefix)){
    colnames(tdf) <- c("disease", "or", "p", "lower", "upper")
  }else{
    colnames(tdf) <- c("disease", paste(stat.prefix, c("or", "p", "lower", "upper"), sep="."))
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
                            shaded.pancan.ci=FALSE, min.y.symmetry=0.1){
  # Get plot values
  or.range <- range(plot.stats[, grep("\\.or$", colnames(plot.stats))], na.rm=T)
  ymin <- min(c(1-min.y.symmetry, or.range[1]-(diff(or.range)/4)))
  ymax <- max(c(1+min.y.symmetry, or.range[2]+(diff(or.range)/4)))

  # Prep plot area
  prep.plot.area(xlims=c(0.25, 2.75), ylims=c(ymin, ymax), parmar=c(2.1, 2.6, 1.1, 0.1), yaxs="r")
  axis(1, at=0.5:2.5, tick=F, line=-0.9, labels=c("AF<1%", "AF<0.1%", "AC=1"),
       cex.axis=5/6)
  mtext("SV frequency", 1, 1)
  y.ax.len <- diff(par("usr")[3:4])
  clean.axis(2, title="Odds ratio", infinite=T)
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
    ci.y.all <- as.numeric(c(plot.stats["pancan", grep("\\.lower$", colnames(plot.stats))],
                             rev(plot.stats["pancan", grep("\\.upper$", colnames(plot.stats))])))
    polygon(x=c(ci.x, rev(ci.x)), y=ci.y.all,
            border=NA, bty="n", col="white")
    polygon(x=c(ci.x, rev(ci.x)), y=ci.y.all,
            border=NA, bty="n", col=adjustcolor(cancer.colors["pancan"], alpha=0.15))
  }

  # Add line for null (OR=1)
  abline(h=1, lty=5)

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
    for(cancer in rownames(plot.stats)){
      points(coords[cancer, paste(freq, c("x", "y"), sep=".")],
             pch=if(cancer == "pancan"){23}else{18},
             cex=if(cancer == "pancan"){1.1}else{0.9},
             bg=if(cancer == "pancan"){cancer.colors[cancer]}else{NA},
             col=if(cancer == "pancan"){"black"}else{cancer.colors[cancer]})
    }
  })
}


#################
# Rscript block #
#################
# Two simple command-line arguments
args <- commandArgs(trailingOnly=T)

# # DEV:
# args <- c("/Users/ryan/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_MS/PedSV_figures/PedSV.v2.5.3_analysis_outputs/stats/PedSV.v2.5.3.global_burden_tests.tsv.gz",
#           "~/scratch/PedSV.v.2.5.3.full_cohort.dev")


# Load sumstats
ss <- load.burden.sumstats(args[1])


# Set standard plot parameters
lineplot.height <- 2.25
lineplot.width <- 2.25
doublebar.height <- 2.25
doublebar.width <- 3.5


# Plot effect sizes of large SVs
sapply(c("", "DEL", "DUP", "INV", "CPX", "unbalanced"), function(suffix){
  plot.stats <- merge.by.freq(lapply(names(freq.names), function(freq){
    hyp <- gsub("\\.$", "", paste("large", freq, suffix, sep="."))
    clean.sumstats.singleHyp(ss, hyp, freq)
  }))
  svtype <- if(suffix == ""){"SV"}else{suffix}
  sv.name <- if(suffix == ""){"SVs"}else{sv.abbreviations[suffix]}
  pdf(paste(args[2], ".large_", svtype, ".or_by_freq.pdf", sep=""),
      height=lineplot.height, width=lineplot.width)
  plot.or.by.freq(plot.stats,
                  title=if(suffix == "unbalanced"){"Unbalanced SVs >1Mb"}else{paste("Large (>1Mb)", tolower(sv.name))},
                  shaded.pancan.ci=T, connect.cancers="pancan")
  dev.off()
})


# Horizontal double-wide barplot of rare unbalanced vs. balanced for main figure
pdf(paste(args[2], ".large_rare_unbal_vs_bal.double_bars.pdf", sep=""),
    height=doublebar.height, width=doublebar.width)
doublewide.barplot.by.phenotype(plot.df.l=stats2barplotdf(ss[which(ss$hypothesis == "large.rare.unbalanced"), ], ci.mode="binomial"),
                                plot.df.r=stats2barplotdf(ss[which(ss$hypothesis == "large.rare.balanced"), ], ci.mode="binomial"),
                                left.axis.units="percent", title="Samples w/rare SV >1Mb",
                                label.l="Unbalanced SVs", label.r="Balanced SVs")
dev.off()


# Horizontal double-wide barplot of rare autosomal vs. allosomal unbalanced SVs for supp figure
pdf(paste(args[2], ".large_rare_unbal_auto_vs_allo.double_bars.pdf", sep=""),
    height=doublebar.height, width=doublebar.width)
doublewide.barplot.by.phenotype(plot.df.l=stats2barplotdf(ss[which(ss$hypothesis == "large.rare.unbalanced.autosomal_only"), ], ci.mode="binomial"),
                                plot.df.r=stats2barplotdf(ss[which(ss$hypothesis == "large.rare.unbalanced.allosomal_only"), ], ci.mode="binomial"),
                                left.axis.units="percent", title="Samples w/unbal. SV >1Mb",
                                label.l="Autosomal", label.r="Allosomal")
dev.off()


# Horizontal double-wide barplot of rare unbalanced SVs in males vs. females for supp figure
pdf(paste(args[2], ".large_rare_unbal_male_vs_female.double_bars.pdf", sep=""),
    height=doublebar.height, width=doublebar.width)
doublewide.barplot.by.phenotype(plot.df.l=stats2barplotdf(ss[which(ss$hypothesis == "large.rare.unbalanced.MALE_only"), ], ci.mode="binomial"),
                                plot.df.r=stats2barplotdf(ss[which(ss$hypothesis == "large.rare.unbalanced.FEMALE_only"), ], ci.mode="binomial"),
                                left.axis.units="percent", title="Samples w/unbal. SV >1Mb",
                                label.l="Males (XY)", label.r="Females (XX)")
dev.off()


# Horizontal double-wide barplot of rare unbalanced autosomal SVs in males vs. females for supp figure
pdf(paste(args[2], ".large_rare_unbal_male_vs_female.autosomes_only.double_bars.pdf", sep=""),
    height=doublebar.height, width=doublebar.width)
doublewide.barplot.by.phenotype(plot.df.l=stats2barplotdf(ss[which(ss$hypothesis == "large.rare.unbalanced.autosomal_only.MALE_only"), ], ci.mode="binomial"),
                                plot.df.r=stats2barplotdf(ss[which(ss$hypothesis == "large.rare.unbalanced.autosomal_only.FEMALE_only"), ], ci.mode="binomial"),
                                left.axis.units="percent", title="Rare autosomal SV >1Mb",
                                label.l="Males (XY)", label.r="Females (XX)")
dev.off()


# Horizontal double-wide barplot of rare unbalanced autosomal SVs in *EUROPEAN* males vs. females for supp figure
pdf(paste(args[2], ".large_rare_unbal_male_vs_female.autosomes_only.EUR_only.double_bars.pdf", sep=""),
    height=doublebar.height, width=doublebar.width)
doublewide.barplot.by.phenotype(plot.df.l=stats2barplotdf(ss[which(ss$hypothesis == "large.rare.unbalanced.autosomal_only.EUR_only.MALE_only"), ], ci.mode="binomial"),
                                plot.df.r=stats2barplotdf(ss[which(ss$hypothesis == "large.rare.unbalanced.autosomal_only.EUR_only.FEMALE_only"), ], ci.mode="binomial"),
                                left.axis.units="percent", title="Rare autosomal SV >1Mb",
                                label.l="Eur. males (XY)", label.r="Eur. females (XX)")
dev.off()


# Horizontal double-wide barplot of rare unbalanced autosomal SVs after excluding COSMIC + CPGs in males vs. females for supp figure
pdf(paste(args[2], ".large_rare_unbal_male_vs_female.autosomes_only.no_COSMIC_no_CPG.double_bars.pdf", sep=""),
    height=doublebar.height, width=doublebar.width)
doublewide.barplot.by.phenotype(plot.df.l=stats2barplotdf(ss[which(ss$hypothesis == "large.rare.unbalanced.autosomal_only.no_COSMIC_no_CPG.MALE_only"), ], ci.mode="binomial"),
                                plot.df.r=stats2barplotdf(ss[which(ss$hypothesis == "large.rare.unbalanced.autosomal_only.no_COSMIC_no_CPG.FEMALE_only"), ], ci.mode="binomial"),
                                left.axis.units="percent", title="Rare autosomal SV >1Mb",
                                label.l="Males (XY)", label.r="Females (XX)")
dev.off()


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
  title <- if(suffix == ""){"Gene-disruptive SVs"}else if(suffix == "nonLoF_disruptive"){"Non-LoF genic SVs"}else{paste(suffix, "SVs")}
  pdf(paste(args[2], ".gene_disruptive_", out.tag, ".or_by_freq.pdf", sep=""),
      height=lineplot.height, width=lineplot.width)
  plot.or.by.freq(plot.stats, title=title, shaded.pancan.ci=T, connect.cancers="pancan")
  dev.off()
})


# Horizontal double-wide barplot of rare LoF vs. CG for main figure
pdf(paste(args[2], ".singleton_lof_vs_cg.double_bars.pdf", sep=""),
    height=doublebar.height, width=doublebar.width)
doublewide.barplot.by.phenotype(plot.df.l=stats2barplotdf(ss[which(ss$hypothesis == "singleton.gene_disruptive.singleton_LoF_SVs"), ]),
                                plot.df.r=stats2barplotdf(ss[which(ss$hypothesis == "singleton.gene_disruptive.singleton_CG_SVs"), ]),
                                title="Singleton SVs per sample",
                                label.l="Loss-of-function (LoF)", label.r="Gene copy gain (CG)",
                                group.label.cex=5/6)
dev.off()


# Horizontal double-wide barplot of rare & singleton gene-disruptive SVs for supp figure
pdf(paste(args[2], ".single_gene_disruptive.rare_and_singleton.double_bars.pdf", sep=""),
    height=doublebar.height, width=doublebar.width)
doublewide.barplot.by.phenotype(plot.df.l=stats2barplotdf(ss[which(ss$hypothesis == "rare.single_gene_disruptive"), ]),
                                plot.df.r=stats2barplotdf(ss[which(ss$hypothesis == "singleton.single_gene_disruptive"), ]),
                                title="Single gene-disruptive SVs",
                                label.l="Rare (AF<1%)", label.r="Singleton (AC=1)")
dev.off()


# Horizontal double-wide barplot of rare gene-disruptive SVs in males vs. females for supp figure
pdf(paste(args[2], ".rare_gene_disruptive.male_vs_female.double_bars.pdf", sep=""),
    height=doublebar.height, width=doublebar.width)
doublewide.barplot.by.phenotype(plot.df.l=stats2barplotdf(ss[which(ss$hypothesis == "rare.gene_disruptive.MALE_only"), ]),
                                plot.df.r=stats2barplotdf(ss[which(ss$hypothesis == "rare.gene_disruptive.FEMALE_only"), ]),
                                title="Gene-disruptive SVs", y.title.line=-0.2,
                                label.l="Male (XY)", label.r="Female (XX)")
dev.off()


# Horizontal double-wide barplot of COSMIC & CPG gene-disruptive SVs for main figure
# Note: custom height for rare only to fit with rest of this specific figure
pdf(paste(args[2], ".rare_cosmic_and_cpg_disruptive.double_bars.pdf", sep=""),
    height=2.75, width=3.75)
doublewide.barplot.by.phenotype(plot.df.l=stats2barplotdf(ss[which(ss$hypothesis == "rare.cosmic.gene_disruptive"), ]),
                                plot.df.r=stats2barplotdf(ss[which(ss$hypothesis == "rare.cpg.gene_disruptive"), ]),
                                title="Samples with rare\ngene-disruptive SV", label.l="COSMIC Tier 1",
                                label.r="Germline CPGs", left.axis.units="percent",
                                group.label.cex=5/6, parmar=c(1.15, 3.25, 0.75, 0.25))
dev.off()
pdf(paste(args[2], ".singleton_cosmic_and_cpg_disruptive.double_bars.pdf", sep=""),
    height=doublebar.height, width=doublebar.width)
doublewide.barplot.by.phenotype(plot.df.l=stats2barplotdf(ss[which(ss$hypothesis == "singleton.cosmic.gene_disruptive"), ]),
                                plot.df.r=stats2barplotdf(ss[which(ss$hypothesis == "singleton.cpg.gene_disruptive"), ]),
                                title="Samples with singleton\ngene-disruptive SV", label.l="COSMIC Tier 1",
                                label.r="Germline CPGs", left.axis.units="percent",
                                group.label.cex=5/6, parmar=c(1.15, 3.25, 0.75, 0.25))
dev.off()

