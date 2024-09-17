#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Plot genome-wide SV site summaries for a single cohort


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
                "INV" = 3/4,
                "CPX" = 1/3)


##################
# Data functions #
##################
# Collect a table of counts per SV type stratified by frequency bins
get.counts.table <- function(bed, af.field="POPMAX_AF", ac.field="AC"){
  sapply(rev(svtypes), function(svtype){
    if(svtype=="CNV"){
      rev(c(0, 0, length(which(bed$SVTYPE == "CNV"))))
    }else{
      single.idx <- filter.bed(bed, query=paste(svtype, "singleton", sep="."),
                               af.field, ac.field, return.idxs=TRUE)
      rare.idx <- filter.bed(bed, query=paste(svtype, "rare", sep="."),
                             af.field, ac.field, return.idxs=TRUE)
      all.idx <- filter.bed(bed, query=svtype, af.field, ac.field, return.idxs=TRUE)
      rev(c(length(single.idx),
            length(setdiff(rare.idx, single.idx)),
            length(setdiff(all.idx, rare.idx))))
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

# Gather Hardy-Weinberg data
make.hwe.mat <- function(bed, prefix=NULL){
  # Get column IDs
  n.genos <- gsub("^_", "", paste(prefix, "N_BI_GENOS", sep="_"))
  n.homref <- gsub("^_", "", paste(prefix, "N_HOMREF", sep="_"))
  n.het <- gsub("^_", "", paste(prefix, "N_HET", sep="_"))
  n.homalt <- gsub("^_", "", paste(prefix, "N_HOMALT", sep="_"))

  # Restrict to autosomal, biallelic variants with at least one non-ref call
  sub.dat <- bed[which(bed$chrom %in% paste("chr", 1:22, sep="") & bed$FILTER == "PASS"), ]
  sub.dat <- sub.dat[which(apply(sub.dat[, c(n.het, n.homalt)], 1, sum) > 0), ]
  HWE.mat <- data.frame("AA"=as.numeric(sub.dat[, n.homref]),
                        "AB"=as.numeric(sub.dat[, n.het]),
                        "BB"=as.numeric(sub.dat[, n.homalt]))
  HWE.mat[complete.cases(HWE.mat), ]
}

# Gather data for large rare SV summary table
get.large.sv.summary.data <- function(bed, ad, meta, cancers){
  # Get data for plotting
  plot.df <- do.call("rbind", lapply(c("DEL", "MONO", "DUP", "POLY", "CPX", "INV", "CTX"), function(svtype){
    if(svtype %in% c("MONO", "POLY")){
      sid.with.aneu <- rownames(meta)[which(meta$any_aneuploidy)]
      auto.ad <- compress.ad.matrix(apply(meta[, paste("chr", 1:22, "_CopyNumber", sep="")], 1, round) - 2, action="any")
      auto.ad[which(!names(auto.ad) %in% sid.with.aneu)] <- 0
      allo.ad <- apply(meta[, paste("chr", c("X", "Y"), "_CopyNumber", sep="")], 1, function(v){round(sum(v))}) - 2
      allo.ad[which(!names(allo.ad) %in% sid.with.aneu)] <- 0
      if(svtype == "MONO"){
        auto.ad <- abs(sapply(auto.ad, function(v){min(c(v, 0))}))
        allo.ad <- abs(sapply(allo.ad, function(v){min(c(v, 0))}))
      }else{
        auto.ad <- abs(sapply(auto.ad, function(v){max(c(v, 0))}))
        allo.ad <- abs(sapply(allo.ad, function(v){max(c(v, 0))}))
      }
      n.auto <- sum(auto.ad)
      n.allo <- sum(allo.ad)
      n.all <- n.auto + n.allo
      carriers <- auto.ad + allo.ad
    }else{
      lr.bed <- filter.bed(bed, query=paste(svtype, "large.rare", sep="."))
      n.all <- nrow(lr.bed)
      n.allo <- length(which(lr.bed$chrom %in% c("chrX", "chrY")))
      n.auto <- n.all - n.allo
      carriers <- query.ad.from.sv.bed(ad, lr.bed, action="any")
    }
    counts <- as.vector(sapply(cancers, function(pheno){
      sids <- rownames(meta)[which(metadata.cancer.label.map[meta$disease] == pheno)]
      sids <- intersect(names(carriers), sids)
      n.carriers <- sum(carriers[sids], na.rm=T)
      c(length(sids), n.carriers, n.carriers / length(sids))
    }))
    c(svtype, n.all, n.auto, n.allo, counts)
  }))
  plot.df <- as.data.frame(plot.df)
  colnames(plot.df) <- c("svtype", "n.all", "n.auto", "n.allo",
                         sapply(cancers, function(pheno){paste(pheno, c("all", "carriers", "freq"), sep=".")}))
  return(plot.df)
}


######################
# Plotting functions #
######################
# Plot a horizontal stacked barplot of SVs colored by frequency bin
plot.count.bars <- function(counts, greyscale=TRUE, parmar=c(0.25, 5.25, 0.1, 2.65),
                            bar.buffer=0.15, count.label.xadj=0.05){
  # Get plot dimensions
  xlims <- c(0, max(apply(counts, 2, sum)))
  ylims <- c(0, ncol(counts))
  label.xadj <- count.label.xadj * diff(xlims)

  # Prep plot area
  prep.plot.area(xlims, ylims, parmar=parmar, xaxs="i", yaxs="i")
  segments(x0=xlims[1], x1=xlims[1], y0=0, y1=ylims[2], col="gray85", xpd=T)

  # Add bars & cap labels
  sapply(1:ncol(counts), function(i){
    x.stops <- c(0, cumsum(counts[, i]))
    bar.pal <- sv.palettes[["DEL"]] # Not matched to SV type here to ensure consistent greyscale for all SVs
    if(greyscale){
      bar.pal <- hex2grey(bar.pal)
    }
    rect(xleft=x.stops[-length(x.stops)], xright=x.stops[-c(1)],
         ybottom=i-1+bar.buffer, ytop=i-bar.buffer,
         col=rev(bar.pal[c("dark2", "main", "light2")]),
         border=rev(bar.pal[c("dark2", "main", "light2")]))
    n.total <- x.stops[length(x.stops)]
    rect(xleft=0, xright=n.total, ybottom=i-1+bar.buffer, ytop=i-bar.buffer,
         col=NA, xpd=T)
    if(n.total > 10000){
      count.label <- paste(round(n.total/1000, 1), "k")
    }else{
      count.label <- prettyNum(n.total, big.mark=",")
    }
    text(x=x.stops[length(x.stops)]-label.xadj, y=i-0.5, labels=count.label,
         pos=4, cex=5/6, xpd=T)
  })

  # Add left axis
  axis(2, at=(1:ncol(counts))-0.5, las=2, tick=F, line=-0.8,
       labels=sv.names[colnames(counts)])

  # Add legend
  legend.x.at <- rep(0.7, 3)*diff(xlims)
  legend.y.at <- c(0.05, 0.15, 0.25)*diff(ylims)
  points(x=legend.x.at, y=legend.y.at, xpd=T, pch=15, cex=1.3,
         col=hex2grey(DEL.colors[c("dark2", "main", "light2")]))
  text(x=legend.x.at, y=legend.y.at-(0.01*diff(ylims)), pos=4, xpd=T,
       labels=rev(c("Common", "Rare", "Singleton")))
}

# Hardy-Weinberg ternary plot
plot.HWE <- function(bed, pop=NULL, title=NULL, full.legend=F, lab.cex=1,
                     colors=c("#4DAC26", "#81F850", "#AC26A1")){
  require(HardyWeinberg, quietly=T)

  #Gather HW p-values & colors
  HWE.mat <- make.hwe.mat(bed, prefix=pop)
  if(nrow(HWE.mat) > 0){
    HW.p <- HardyWeinberg::HWChisqStats(X=HWE.mat, x.linked=F, pvalues=T)
    HW.cols <- rep(colors[1], times=length(HW.p))
    HW.cols[which(HW.p<0.05)] <- colors[2]
    HW.cols[which(HW.p<0.05/length(HW.p))] <- colors[3]
  }

  # Generate HW plot frame
  par(mar=c(1, 1, 1, 1), bty="n")
  plot(x=1.15*c(-1/sqrt(3), 1/sqrt(3)), y=c(-0.15, 1.15), type="n",
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
  segments(x0=c(-1/sqrt(3), 0, 1/sqrt(3)),
           x1=c(0, 1/sqrt(3), -1/sqrt(3)),
           y0=c(0, 1, 0), y1=c(1, 0, 0))
  if(nrow(HWE.mat) > 0){
    HWTernaryPlot(X=HWE.mat, n=max(HWE.mat, na.rm=T), newframe=F,
                  vbounds=F, mafbounds=F,
                  region=1, vertexlab=NA,
                  alpha=0.05,
                  curvecols=c(colors[1:2], NA, NA), pch=NA)
  }

  #Add axes
  text(x=c(-1/sqrt(3), 1/sqrt(3)), y=0, labels=c("Ref/Ref", "Alt/Alt"),
       pos=1, cex=0.8, xpd=T, font=2)
  text(x=0, y=1, labels="Ref/Alt", pos=3, cex=0.8, xpd=T, font=2)

  #Finish HW plot
  if(nrow(HWE.mat) > 0){
    HWTernaryPlot(X=HWE.mat, n=max(HWE.mat, na.rm=T), newframe=F,
                  vbounds=F, mafbounds=F,
                  region=1, vertexlab=NA,
                  alpha=0.03/nrow(HWE.mat),
                  curvecols=c(colors[1], colors[3], NA, NA),
                  pch=21, cex=0.3, signifcolour=F, markercol=NA,
                  markerbgcol=adjustcolor(HW.cols, alpha=0.25))
  }
  segments(x0=c(-1/sqrt(3), 0, 1/sqrt(3)),
           x1=c(0, 1/sqrt(3), -1/sqrt(3)),
           y0=c(0, 1, 0), y1=c(1, 0, 0))

  #Add legend
  if(nrow(HWE.mat) > 0){
    n.pass <- length(which(HW.p>=0.05))
    cat(paste("PASS: ", n.pass/length(HW.p), "\n", sep=""))
    n.nom <- length(which(HW.p<0.05 & HW.p>=0.05/nrow(HWE.mat)))
    cat(paste("NOMINAL FAILS: ", n.nom/length(HW.p), "\n", sep=""))
    n.bonf <- length(which(HW.p<0.05/nrow(HWE.mat)))
    cat(paste("BONFERRONI FAILS: ", n.bonf/length(HW.p), "\n", sep=""))
    n.samples <- max(apply(HWE.mat, 1, sum), na.rm=T)
    n.SV <- nrow(HWE.mat)
    legend("right", pch=19, col=colors, pt.cex=1.3,
           legend=c(paste(round(100*(n.pass/nrow(HWE.mat)), 0), "%", sep=""),
                    paste(round(100*(n.nom/nrow(HWE.mat)), 0), "%", sep=""),
                    paste(round(100*(n.bonf/nrow(HWE.mat)), 0), "%", sep="")),
           bty="n", bg=NA, cex=0.7)
  }else{
    cat(paste("No SVs with AC>0 found for population '", pop, "'\n", sep=""))
    n.samples <- 0
    n.SV <- 0
  }
  text(x=par("usr")[2], y=par("usr")[4]-(0.2*(par("usr")[4]-par("usr")[3])),
       pos=2, cex=0.7, labels=paste(title, "\n \n ", sep=""), font=2)
  text(x=par("usr")[2], y=par("usr")[4]-(0.2*(par("usr")[4]-par("usr")[3])),
       pos=2, cex=0.7,
       labels=paste(" \n", prettyNum(n.samples, big.mark=","),
                    " Samples\n ", sep=""))
  text(x=par("usr")[2], y=par("usr")[4]-(0.2*(par("usr")[4]-par("usr")[3])),
       pos=2, cex=0.7, labels=paste(" \n \n", prettyNum(n.SV, big.mark=","),
                                    " SVs", sep=""))
}

# Plot a summary table graphic of large, rare SVs
plot.large.sv.summary <- function(bed, ad, meta, schematic.wex=12, count.wex=2,
                                  freq.wex=3, svtype.x.adj=0.5, w.buffer=0.25,
                                  title.cex=5/6, content.cex=4/6){
  # Get data for plotting
  cancers <- intersect(names(cancer.colors), metadata.cancer.label.map[meta$disease])
  plot.df <- get.large.sv.summary.data(bed, ad, meta, cancers)

  # Set plot parameters
  n.rows <- nrow(plot.df)
  n.cancers <- length(cancers)
  x.breaks <- c(0, schematic.wex,
                schematic.wex + (count.wex * 1:3),
                schematic.wex + (3 * count.wex) + (freq.wex * 1:n.cancers))

  # Prepare plot area
  prep.plot.area(xlims=range(x.breaks), ylims=c(n.rows, 0),
                 parmar=c(0.1, 0.1, 2.1, 0.1))

  # Shade relevant columns for each disease
  sapply(1:n.cancers, function(i){
    rect(xleft=x.breaks[4+i]+w.buffer, xright=x.breaks[5+i]-w.buffer,
         ybottom=par("usr")[3], ytop=par("usr")[4], bty="n", border=NA,
         col=cancer.palettes[[cancers[i]]]["light3"])
  })

  # Add row lines in light grey
  segments(x0=x.breaks[1]+w.buffer, x1=x.breaks[length(x.breaks)]-w.buffer,
           y0=1:n.rows, y1=1:n.rows, col="gray95")

  # Add column headers
  clean.axis(3, at=x.breaks[1:2] + (c(-1, 1) * w.buffer),
             tck=0, labels=NA, label.line=-0.9, cex.axis=1)
  text(x=0-svtype.x.adj, y=-0.5, xpd=T, pos=4, labels="SV type", cex=title.cex)
  clean.axis(3, tck=0, labels=c(NA, "Total", "Auto.", "Allo.", NA),
             label.line=-0.9, cex.axis=content.cex,
             at=c(x.breaks[2]+w.buffer,
                     sapply(2:4, function(x){mean(x.breaks[c(x, x+1)])}),
                  x.breaks[5]-w.buffer))
  axis(3, at=mean(x.breaks[c(2, 5)]), labels="Rare SVs >1 Mb", tick=F, cex.axis=title.cex)
  clean.axis(3, tck=0, labels=c(NA, cancer.names.vshort[cancers], NA),
             label.line=-0.9, cex.axis=content.cex,
             at=c(x.breaks[5]+w.buffer,
                  sapply(5:(length(x.breaks)-1), function(x){mean(x.breaks[c(x, x+1)])}),
                  x.breaks[length(x.breaks)]-w.buffer))
  axis(3, at=mean(x.breaks[c(5, length(x.breaks))]), labels="Carrier rate",
       tick=F, cex.axis=title.cex)

  # Add content for each row
  sapply(1:n.rows, function(i){
    if(plot.df$svtype[i] == "MONO"){
      sv.label <- "Monosomy"
    }else if(plot.df$svtype[i] == "POLY"){
      sv.label <- "Polysomy"
    }else{
      sv.label <- sv.names[plot.df$svtype[i]]
    }
    text(x=x.breaks[1]-svtype.x.adj, y=i-0.5, pos=4, labels=sv.label, cex=5/6)
    text(x=sapply(2:4, function(x){mean(x.breaks[c(x, x+1)])}),
         y=i-0.5, labels=plot.df[i, 2:4], cex=content.cex)
    freq.labels <- sapply(cancers, function(cancer){
      paste(round(100 * as.numeric(plot.df[i, paste(cancer, "freq", sep=".")]), 2), "%", sep="")
    #         "%\n(",
    #     prettyNum(plot.df[i, paste(cancer, "carriers", sep=".")], big.mark=","),
    #         ")", sep="")
    })
    text(x=sapply(5:(length(x.breaks)-1), function(x){mean(x.breaks[c(x, x+1)])}),
         y=i-0.5, cex=content.cex, labels=freq.labels)
  })
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Plot SV counts and sizes")
parser$add_argument("--bed", metavar=".tsv", type="character", required=TRUE,
                    help="SV sites .bed")
parser$add_argument("--ad-matrix", metavar=".tsv", type="character", required=TRUE,
                    help="Allele dosage .bed.")
parser$add_argument("--metadata", metavar=".tsv", type="character",
                    help="sample metadata .tsv", required=TRUE)
parser$add_argument("--subset-samples", metavar=".tsv", type="character",
                    help="list of samples to subset [default: use all samples]")
parser$add_argument("--cohort-prefix", default="", metavar="string", type="character",
                    help="String prefix to append to frequency columns")
parser$add_argument("--af-field", default="AF", metavar="string", type="character",
                    help=paste("Column header to use for AF-related analyses.",
                               "Overrides --cohort-prefix."))
parser$add_argument("--ac-field", default="AC", metavar="string", type="character",
                    help=paste("Column header to use for AF-related analyses.",
                               "Overrides --cohort-prefix."))
parser$add_argument("--out-prefix", metavar="path", type="character", required=TRUE,
                    help="Path/prefix for all output files")
args <- parser$parse_args()

# # DEV:
# args <- list("bed" = "~/scratch/PedSV.v2.5.4.full_cohort_w_relatives_1000G.sites.bed.gz",
#              "ad_matrix" = "~/scratch/PedSV.v2.5.4.full_cohort_w_relatives_1000G.allele_dosages.bed.gz",
#              "metadata" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_v2_callset_generation/v2.5.4/PedSV.v2.5.4.cohort_metadata.w_control_assignments.tsv.gz",
#              # "subset_samples" = "~/scratch/PedSV.v2.5.3.final_analysis_cohort.samples.list",
#              "subset_samples" = NULL,
#              "cohort_prefix" = "full_cohort",
#              "af_field" = "POPMAX_AF",
#              "ac_field" = "AC",
#              "out_prefix" = "~/scratch/PedSV.v2.5.4.full_cohort_w_relatives")
# args <- list("bed" = "~/scratch/YL-gatsv-v1-allBatches.annotated.samples_excluded.bed.gz",
#              "metadata" = "/Users/ryan/Desktop/Collins/VanAllen/jackie_younglung/younglung_metadata/YL.SV.v1.2.analysis_metadata.tsv.gz",
#              "subset_samples" = "~/Desktop/Collins/VanAllen/jackie_younglung/data/ancestry_and_relatedness/YL.SV.v1.2.analysis.samples.list",
#              "cohort_prefix" = "",
#              "af_field" = "POPMAX_AF",
#              "ac_field" = "AC",
#              "out_prefix" = "~/scratch/YL.v1.2")


# Infer frequency columns to use
if(is.null(args$af_field)){
  args$af_field <- paste(args$cohort_prefix, "AF", sep="_")
}
if(is.null(args$ac_field)){
  args$ac_field <- paste(args$cohort_prefix, "AC", sep="_")
}
if(args$cohort_prefix == "trio"
   & length(grep("case_control", c(args$af_field, args$ac_field), fixed=T)) == 0){
  drop.cohort.freqs <- c("case_control")
}else if(args$cohort_prefix == "case_control"
         & length(grep("trio", c(args$af_field, args$ac_field), fixed=T)) == 0){
  drop.cohort.freqs <- c("trio")
}else{
  drop.cohort.freqs <- c()
}

# Load BED
bed <- PedSV::load.sv.bed(args$bed, drop.cohort.frequencies=drop.cohort.freqs,
                          keep.all.pop.frequencies=TRUE)

# Load metadata
keepers <- NULL
if(!is.null(args$subset_samples)){
  keepers <- read.table(args$subset_samples, header=F)[, 1]
}
meta <- load.sample.metadata(args$metadata, keep.samples=keepers,
                             reassign.parents=FALSE)

# Plot stacked bar colored by frequency and type
counts.by.freq <- get.counts.table(bed, args$af_field, args$ac_field)
pdf(paste(args$out_prefix, "sv_site_counts.pdf", sep="."),
    height=1.5, width=2.7)
plot.count.bars(counts.by.freq)
dev.off()

# Ridgeplot of SV sizes
pdf(paste(args$out_prefix, "sv_size_distribs.pdf", sep="."),
    height=1.7, width=2.2)
ridgeplot(get.svlen.densities(bed), xlims=log10(c(10, 5000000)), x.axis=FALSE,
          fill=hex2grey(DEL.colors[["light2"]]),
          border=hex2grey(DEL.colors[["dark1"]]), border.lwd=1.25,
          parmar=c(2.2, 3.5, 0.1, 0.1))
clean.axis(1, at=log10(logscale.major.bp),
           labels=logscale.major.bp.labels[seq(1, length(logscale.major.bp), 2)],
           labels.at=log10(logscale.major.bp)[seq(1, length(logscale.major.bp), 2)],
           label.line=-0.9, title.line=0.2, title=bquote("SV size" ~ (log[10])))
dev.off()

# Table of SV counts by type and frequency
write.table(tabulate.counts(bed),
            paste(args$out_prefix, "sv_counts_by_freq.tsv", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)

# Table of median SV sizes by type and frequency
write.table(tabulate.sizes(bed),
            paste(args$out_prefix, "sv_sizes_by_freq.tsv", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)

# HWE ternary plots per population
pops.in.bed <- unique(sapply(colnames(bed)[grep("_AF$", colnames(bed))],
                             function(col.name){unlist(strsplit(col.name, split="_"))[1]}))
for(pop in intersect(names(pop.colors), pops.in.bed)){
  col.prefix <- gsub("^_", "", paste(args$cohort_prefix, pop, sep="_"))
  if(!col.prefix %in% colnames(bed)){
    col.prefix <- pop
  }
  if(paste(col.prefix, "AF", sep="_") %in% colnames(bed)){
    cat(paste("HWE for ", pop, ":\n", sep=""))
    tiff(paste(args$out_prefix, pop, "HWE.tiff", sep="."),
         height=1300, width=1300, res=400)
    plot.HWE(bed, pop=col.prefix, title=pop.names.short[pop], full.legend=T)
    dev.off()
  }
}

# Table of very large (>1Mb) SVs
pdf(paste(args$out_prefix, "large_rare_SV_table.pdf", sep="."), height=2.5, width=4.3)
plot.large.sv.summary(bed, args$ad_matrix, meta)
dev.off()


