#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2024-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Plot allele frequency correlations between PedSV callset and gnomAD


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(PedSV, quietly=TRUE)
require(argparse, quietly=TRUE)
PedSV::load.constants("all")

# Declare local constants
gnomad.pop.map <- c("AFR" = "AFR",
                    "AMR" = "AMR",
                    "EAS" = "EAS",
                    "EUR" = "NFE",
                    "SAS" = "SAS")


######################
# Plotting functions #
######################
# Plot AF correlation vs. gnomAD for a single population
plot.af.comparison <- function(bed, pop, cohort, pt.cex=0.1, bandwidth=2, alpha=1){
  # Get plot data
  bed <- filter.bed(bed, query="", autosomal=TRUE, pass.only=TRUE)
  x <- bed[, gsub("^_", "", paste(cohort, pop, "AF", sep="_"))]
  y <- bed[, paste("gnomad_v3.1_sv", gnomad.pop.map[pop], "AF", sep="_")]
  keepers <- which(!is.na(x) & !is.na(y) & x > 0)
  x <- log10(x[keepers]); y <- log10(y[keepers])
  x.min <- min(x)
  y.min <- min(y[which(!is.infinite(y))])
  ax.lims <- c(min(c(x.min, y.min)), log10(1))

  # Prep plot area
  prep.plot.area(xlims=ax.lims, ylims=ax.lims, parmar=c(2.5, 3.4, 1, 0.1),
                 xaxs="r", yaxs="r")
  abline(0, 1, col="gray30")
  rect(xleft=par("usr")[1], xright=x.min, ybottom=par("usr")[3], ytop=par("usr")[4],
       col="gray90", border=NA)
  rect(xleft=par("usr")[1], xright=x.min, ybottom=par("usr")[3], ytop=par("usr")[4],
       col="white", border=NA, density=4)
  text(x=mean(c(par("usr")[1], x.min)), y=mean(par("usr")[3:4]),
       labels="Below AF resolution", col="gray65", cex=5/6,
       srt=(180 / pi) * atan(diff(par("usr")[3:4])/(x.min-par("usr")[1])))
  label.denoms <- gsub("000$", "k", gsub("000000$", "M", as.character(ceiling(1/logscale.major))))
  ax.labels <- gsub("^1:1$", "1", paste("1", label.denoms, sep=":"))
  clean.axis(1, at=log10(logscale.major), title="Allele freq. (this study)",
             labels=ax.labels)
  clean.axis(2, at=log10(logscale.major), title="Allele freq. (gnomAD)",
             labels=ax.labels, title.line=1.3)
  mtext(3, text=pop.names.long[pop], font=2, line=0)

  # Add points colored by SV type to make common insertion
  plot.df <- color.points.by.density(x, y, bandwidth=bandwidth,
                                     palette=colorRampPalette(pop.palettes[[pop]][4:2])(256))
  # plot.df <- data.frame("x"=x, "y"=y, "col"=as.character(pop.colors[pop]))
  points(plot.df$x, plot.df$y, col=adjustcolor(plot.df$col, alpha=alpha), pch=19, cex=pt.cex, xpd=T)

  # Add correlation coefficient below title
  cor.stats <- cor.test(10^x, 10^y, use="complete.obs")
  r <- cor.stats$estimate^2
  mtext(3, text=bquote(italic(R)^2 == .(formatC(round(r, 3), digits=3))),
        line=-1, cex=5.5/6)
  return(cor.stats)
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Compare SV AFs to gnomAD")
parser$add_argument("bed", metavar=".tsv", type="character",
                    help="SV sites .bed")
parser$add_argument("--cohort-prefix", default="", metavar="string", type="character",
                    help="String prefix to append to frequency columns")
parser$add_argument("--out-prefix", metavar="path", type="character", required=TRUE,
                    help="Path/prefix for all output files")
args <- parser$parse_args()

# # DEV:
# args <- list("bed" = "~/scratch/PedSV.v2.5.3.full_cohort.analysis_samples.sites.bed.gz",
#              "cohort_prefix" = "",
#              "out_prefix" = "~/scratch/PedSV.v2.5.3.dev.full_cohort")


# Infer frequency columns to use
drop.cohort.freqs <- c()
if(!is.null(args$cohort_prefix)){
  if(args$cohort_prefix == "trio"
     & length(grep("case_control", c(args$af_field, args$ac_field), fixed=T)) == 0){
    drop.cohort.freqs <- c("case_control")
  }else if(args$cohort_prefix == "case_control"
           & length(grep("trio", c(args$af_field, args$ac_field), fixed=T)) == 0){
    drop.cohort.freqs <- c("trio")
  }
}

# Load BED
bed <- PedSV::load.sv.bed(args$bed, drop.cohort.frequencies=drop.cohort.freqs,
                          keep.all.pop.frequencies=TRUE)

# Plot one AF correlation per population
pops.in.bed <- unique(sapply(colnames(bed)[grep("_AF$", colnames(bed))],
                             function(col.name){unlist(strsplit(col.name, split="_"))[1]}))
for(pop in intersect(names(pop.colors), pops.in.bed)){
  col.prefix <- gsub("^_", "", paste(args$cohort_prefix, pop, sep="_"))
  if(!col.prefix %in% colnames(bed)){
    col.prefix <- pop
  }
  if(paste(col.prefix, "AF", sep="_") %in% colnames(bed)){
    cat(paste("gnomAD AF comparison for ", pop.names.long[pop], ":\n", sep=""))
    tiff(paste(args$out_prefix, pop, "vs_gnomad.tiff", sep="."),
         height=1300, width=1300, res=400)
    cor.stats <- plot.af.comparison(bed, pop, args$cohort_prefix)
    dev.off()
    cat(paste("R2 =", cor.stats$estimate^2, "\nP =", cor.stats$p.value, "\n"))
  }
}

