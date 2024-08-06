#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Analysis of de novo SV properties in cancer trios for manuscript


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(PedSV, quietly=TRUE)
require(DescTools, quietly=TRUE)
PedSV::load.constants("all")


##################
# Data functions #
##################
# Load list of de novo SVs and add data from sites BED
load.denovo.svs <- function(dn.tsv, bed, meta){
  df <- read.table(dn.tsv, header=F, sep="\t")
  colnames(df) <- c("VID", "sample", "origin")
  df <- merge(df, meta, all.x=T, all.y=F, by.x="sample", by.y="row.names", sort=F)
  df <- merge(df, bed, all.x=T, all.y=F, by.x="VID", by.y="row.names", sort=F)
  return(df)
}

# Load published control de novo counts and compute rates/CI
load.control.rates <- function(counts.tsv){
  df <- read.table(counts.tsv, header=T, sep="\t")
  rates <- PoissonCI(x=df$denovo, n=df$trios)
  colnames(rates) <- c("estimate", "lower95", "upper95")
  df <- as.data.frame(cbind(df, rates))
  oth.rates <- as.vector(apply(rates[which(df$svtype %in% c("INV", "CPX")), ], 2, sum))
  df <- as.data.frame(rbind(df[which(df$svtype %in% c("DEL", "DUP", "INS")), ],
                            c("OTH", max(df$trios),
                              sum(df$denovo[which(df$svtype %in% c("INV", "CPX"))]),
                              oth.rates)))
  df[, -c(1)] <- apply(df[, -c(1)], 2, as.numeric)
  return(df)
}

# Load gnomAD de novo mutation rates
load.gnomad.rates <- function(gnomad.tsv){
  df <- read.table(gnomad.tsv, header=T, sep="\t")
  oth.rates <- as.vector(apply(df[which(df$svtype %in% c("INV", "CPX")), 2:4], 2, sum))
  df <- as.data.frame(rbind(df[which(df$svtype %in% c("DEL", "DUP", "INS")), ],
                            c("OTH", oth.rates)))
  df[, c("svtype", "estimate", "lower95", "upper95")]
}

# Run all pairwise case-control rate comparisons
all.case.control.comparisons <- function(denovos, meta, control.counts){
  control.n <- max(control.counts$trios)
  res <- do.call("rbind", lapply(metadata.cancer.label.map[unique(denovos$disease)],
                                 function(cancer){
                                   dn.sub <- denovos[which(metadata.cancer.label.map[denovos$disease] == cancer), ]
                                   case.n <- length(get.eligible.samples(meta, cancer)$cases)
                                   do.call("rbind", lapply(c("ALL", "DEL", "DUP", "INS", "OTH"), function(svtype){
                                     if(svtype == "ALL"){
                                       case.k <- nrow(dn.sub)
                                       control.k <- sum(control.counts$denovo)
                                     }else{
                                       control.k <- control.counts$denovo[which(control.counts$svtype == svtype)]
                                       if(svtype == "OTH"){
                                         case.k <- length(which(!(dn.sub$SVTYPE %in% c("DEL", "DUP", "INS"))))
                                       }else{
                                         case.k <- length(which(dn.sub$SVTYPE == svtype))
                                       }
                                     }
                                     res <- fisher.test(matrix(c(control.n-control.k, case.n-case.k,
                                                                 control.k, case.k),
                                                               byrow=T, nrow=2))
                                     c(cancer, svtype, case.k, control.k, as.numeric(unlist(res))[c(1, 4, 2:3)])
                                   }))
                                 }))
  res <- as.data.frame(res)
  colnames(res) <- c("cancer", "svtype", "denovo.case", "denovo.control", "pvalue",
                     "OR", "OR_lower", "OR_upper")
  res[order(res$pvalue), ]
}

# Compute de novo rates for a single cancer type
calc.denovo.rates <- function(denovos, meta, cancer="pancan"){
  elig.cases <- get.eligible.samples(meta, cancer)$cases
  n.children <- length(elig.cases)
  denovos <- denovos[which(denovos$sample %in% elig.cases), ]
  k <- c(sum(denovos$SVTYPE == "DEL", na.rm=T),
         sum(denovos$SVTYPE == "DUP", na.rm=T),
         sum(denovos$SVTYPE == "INS", na.rm=T),
         sum(!(denovos$SVTYPE %in% c("DEL", "DUP", "INS")), na.rm=T))
  rates <- as.data.frame(PoissonCI(x=k, n=n.children))
  colnames(rates) <- c("estimate", "lower95", "upper95")
  rates$svtype <- c("DEL", "DUP", "INS", "OTH")
  rates[, c("svtype", "estimate", "lower95", "upper95")]
}

# Explode table of SV frequency vs. genic context
get.freq.by.context <- function(freq.idxs, bed){
  # Get tiers of consequence
  context.idxs <- list("gene_disruptive" = filter.bed(bed, "gene_disruptive", return.idxs=TRUE))
  context.idxs[["other_coding"]] <- setdiff(filter.bed(bed, "coding", return.idxs=TRUE), unlist(context.idxs))
  context.idxs[["intronic"]] <- setdiff(filter.bed(bed, "intronic", return.idxs=TRUE), unlist(context.idxs))
  context.idxs[["intergenic"]] <- setdiff(filter.bed(bed, "intergenic", return.idxs=TRUE), unlist(context.idxs))

  # Stratify consequence by freqs
  as.data.frame(do.call("rbind", lapply(freq.idxs, function(f.idxs){
    sapply(context.idxs, function(c.idxs){length(intersect(f.idxs, c.idxs))})
  })))
}


######################
# Plotting functions #
######################
# Plot de novo rates in our data vs. published expectations
plot.denovo.rates <- function(denovos, control.rates, gnomad.rates,
                              parmar=c(1, 2.8, 0.25, 0)){
  # Simplify rates for insertions and inv/cpx/other
  ews.rates <- calc.denovo.rates(denovos, meta, "EWS")
  nbl.rates <- calc.denovo.rates(denovos, meta, "NBL")
  plot.data <- list("EWS" = ews.rates,
                    "NBL" = nbl.rates,
                    "control" = control.rates,
                    "gnomad" = gnomad.rates)
  x.buffer <- seq(-0.9, -0.1, length.out=length(plot.data) + 2)[-c(1, length(plot.data) + 2)]

  # Get plot dimensions
  n.x <- nrow(ews.rates)
  n.groups <- length(plot.data)
  ymax.ci <- max(as.numeric(unlist(do.call("rbind", plot.data)[, 2:4]), na.rm=T))
  ymax.meanplus <- (5 * max(as.numeric(unlist(do.call("rbind", plot.data)[, 2]), na.rm=T))) / 3
  ymax <- min(c(ymax.ci, ymax.meanplus))

  # Prepare plot area
  prep.plot.area(xlims=c(0, n.x), ylims=c(0, ymax), parmar=parmar)
  sapply(1:n.x, function(x){
    axis(1, at=x+x.buffer[c(1, n.groups)], tck=0, labels=NA,
         col=cancer.colors[["control"]])
    axis(1, at=x-0.5, tick=F, line=-0.9, labels=sv.abbreviations[ews.rates$svtype][x])
  })

  clean.axis(2, title=bquote("Rate of" ~ italic("se novo") ~ "SVs"),
             title.line=0.8, infinite.positive=TRUE)

  # Add data for each group
  sapply(1:n.groups, function(i){
    # Get color
    g.name <- names(plot.data)[i]
    if(g.name %in% names(cancer.colors)){
      color <- cancer.colors[g.name]
    }else if(g.name == "gnomad"){
      color <- "#C5BAA7"
    }

    # Add bars & points
    segments(x0=(1:n.x)+x.buffer[i], x1=(1:n.x)+x.buffer[i],
             y0=as.numeric(plot.data[[i]]$lower95),
             y1=as.numeric(plot.data[[i]]$upper95),
             lwd=3, lend="butt", col=color)
    points(x=(1:n.x)+x.buffer[i], y=plot.data[[i]]$estimate,
           pch=19, col=color)
  })
}

# Plot a horizontal stacked barplot of coding contexts per frequency bin
plot.context.by.freq <- function(counts, group.names=NULL,
                                 bar.buffer=0.15, count.label.xadj=0.04,
                                 parmar=c(0.25, 3.5, 2.5, 0.25)){
  # Reorder counts matrix
  counts <- counts[nrow(counts):1, ]

  # Normalize
  props <- as.data.frame(t(apply(counts, 1, function(v){v / sum(v)})))[, 1:3]

  # Get plot dimensions & other properties
  xlims <- c(0, max(apply(props, 1, sum)))
  ylims <- c(-2, nrow(props))
  label.xadj <- count.label.xadj * diff(xlims)
  if(!is.null(group.names)){
    group.names <- rev(group.names)
  }else{
    group.names <- rownames(counts)
  }
  bar.pal <- context.colors[colnames(props)]

  # Prep plot area
  prep.plot.area(xlims, ylims, parmar=parmar, xaxs="i", yaxs="i")
  segments(x0=xlims[1], x1=xlims[1], y0=0, y1=ylims[2], col="gray85", xpd=T)

  # Add bars & cap labels
  sapply(1:nrow(props), function(i){
    x.stops <- unlist(c(0, cumsum(as.numeric(props[i, ]))))

    rect(xleft=x.stops[-length(x.stops)], xright=x.stops[-c(1)],
         ybottom=i-1+bar.buffer, ytop=i-bar.buffer,
         col=bar.pal, border="white", lwd=0.5)
    n.total <- x.stops[length(x.stops)]
    rect(xleft=0, xright=n.total, ybottom=i-1+bar.buffer, ytop=i-bar.buffer,
         col=NA, xpd=T)
    text(x=x.stops[length(x.stops)-1]-label.xadj, y=i-0.5, pos=4, cex=5/6, xpd=T,
         labels=paste(round(100 * x.stops[length(x.stops)-1], 0), "%", sep=""))
  })

  # Add left & vertical axes
  axis(2, at=(1:nrow(props))-0.5, las=2, tick=F, line=-0.8, labels=group.names)
  clean.axis(3, label.units="percent", title="Proportion of SVs",
             title.line=0, label.line=-0.8, infinite.positive=TRUE)

  # Add legend
  legend.x.at <- c(0, 0, 0)*diff(xlims)
  legend.y.at <- c(-0.3, -0.9, -1.5)-0.2
  points(x=legend.x.at, y=legend.y.at, xpd=T, pch=15, cex=1.3, col=bar.pal)
  text(x=legend.x.at, y=legend.y.at, pos=4, xpd=T,
       labels=c("Gene-disruptive", "Other coding", "Intronic"))
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Plot de novo SV properties")
parser$add_argument("--bed", metavar=".tsv", type="character", required=TRUE,
                    help="SV sites .bed")
parser$add_argument("--metadata", metavar=".tsv", type="character",
                    help="sample metadata .tsv", required=TRUE)
parser$add_argument("--de-novos", metavar=".tsv", type="character", required=TRUE,
                    help=paste("three-column .tsv of variant_id, proband_id, ",
                               "denovo|mosaic for confirmed de novo SV events"))
parser$add_argument("--include-families", metavar=".txt", type="character",
                    help="list of families to be included in de novo statistics",
                    required=TRUE)
parser$add_argument("--control-counts", metavar=".tsv", type="character", required=TRUE,
                    help=paste(".tsv of de novo counts for published controls"))
parser$add_argument("--gnomad-rates", metavar=".tsv", type="character", required=TRUE,
                    help=paste(".tsv of estimated SV mutation rates from gnomAD"))
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
# args <- list("bed" = "~/scratch/PedSV.v2.5.3.full_cohort.analysis_samples.sites.bed.gz",
#              "metadata" = "~/scratch/PedSV.v2.5.3.cohort_metadata.w_control_assignments.tsv.gz",
#              "de_novos" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/de_novo_analysis/PedSV.v2.5.3.confirmed_de_novo_events.tsv",
#              "include_families" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/de_novo_analysis/PedSV.v2.5.3.complete_trios.no_outliers.list",
#              "control_counts" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/de_novo_analysis/refs/belyeu_denovo_counts.tsv",
#              "gnomad_rates" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/de_novo_analysis/refs/gnomad_sv_mutation_rates.tsv",
#              "cohort_prefix" = "",
#              "af_field" = "POPMAX_AF",
#              "ac_field" = "AC",
#              "out_prefix" = "~/scratch/PedSV.v2.5.3.dev.complete_trios")

# Infer frequency columns to use
if(is.null(args$af_field)){
  args$af_field <- paste(args$cohort_prefix, "AF", sep="_")
}
if(is.null(args$ac_field)){
  args$ac_field <- paste(args$cohort_prefix, "AC", sep="_")
}

# Load BED
bed <- PedSV::load.sv.bed(args$bed,
                          drop.cohort.frequencies=c("trio", "case_control"))

# Load metadata
meta <- load.sample.metadata(args$metadata, reassign.parents=FALSE)
fams <- read.table(args$include_families, header=F)[, 1]
meta <- meta[which(meta$family_id %in% fams), ]

# Load & format de novo SVs
denovos <- load.denovo.svs(args$de_novos, bed, meta)

# NBL vs. EWS comparison (all de novo SVs)
nbl.n <- length(which(meta$disease == "neuroblastoma" & meta$proband))
nbl.dn <- length(which(denovos$disease == "neuroblastoma"))
ews.n <- length(which(meta$disease == "ewing" & meta$proband))
ews.dn <- length(which(denovos$disease == "ewing"))
fisher.test(matrix(c(ews.n-ews.dn, nbl.n-nbl.dn, ews.dn, nbl.dn),
                   byrow=T, nrow=2))

# Load published control rates
control.data <- load.control.rates(args$control_counts)
control.counts <- control.data[, c("svtype", "trios", "denovo")]
control.rates <- control.data[, c("svtype", "estimate", "lower95", "upper95")]
gnomad.rates <- load.gnomad.rates(args$gnomad_rates)

# Run all pairwise rate comparisons between cancer + control
cc.stats <- all.case.control.comparisons(denovos, meta, control.counts)
print(cc.stats)
write.table(cc.stats, paste(args$out_prefix, "denovo_SVs", "case_control_comparisons",
                            "vs_published_controls", "tsv", sep="."),
            col.names=T, row.names=F, quote=F, sep="\t")

# Plot comparison of rates for selected SV types
pdf(paste(args$out_prefix, "de_novo.rates.pdf", sep="."),
    height=1.8, width=3)
plot.denovo.rates(denovos, control.rates, gnomad.rates)
dev.off()

# Plot size comparison of de novo SVs vs. rare & common SVs (excluding insertions)
freq.idxs <- list("denovo" = which(rownames(bed) %in% denovos$VID))
freq.idxs[["singleton"]] <- setdiff(filter.bed(bed, "singleton", return.idxs=TRUE),
                                    unlist(freq.idxs))
freq.idxs[["rare"]] <- setdiff(filter.bed(bed, "rare", return.idxs=TRUE),
                               unlist(freq.idxs))
freq.idxs[["common"]] <- setdiff(filter.bed(bed, "", return.idxs=TRUE),
                                 unlist(freq.idxs))
size.data <- lapply(freq.idxs, function(idxs){
  log10(bed[setdiff(idxs, which(bed$SVTYPE %in% c("INS", "CTX"))), "SVLEN"])
})
pdf(paste(args$out_prefix, "de_novo.sizes.pdf", sep="."),
    height=1.9, width=2.5)
ridgeplot(lapply(rev(size.data), density, adjust=0.75),
          names=rev(c("", "Singleton", "Rare", "Common")), xlims=log10(c(10, 25000000)),
          fill=rev(c("black", hex2grey(DEL.colors[c("dark2", "main", "light2")]))),
          border=rev(c("black", rep(hex2grey(DEL.colors[["dark1"]]), 3))),
          border.lwd=1.25, x.axis=FALSE, parmar=c(2.2, 3.75, 0.25, 0.5))
axis(2, at=3.5, labels=bquote(italic("De novo")), tick=F, las=2, line=-0.8)
clean.axis(1, at=log10(logscale.major.bp),
           labels=logscale.major.bp.labels[seq(1, length(logscale.major.bp), 2)],
           labels.at=log10(logscale.major.bp)[seq(1, length(logscale.major.bp), 2)],
           label.line=-0.9, title.line=0.2, title=bquote("SV size" ~ (log[10])))
dev.off()
RLCtools::format.pval(wilcox.test(size.data[["denovo"]], unlist(size.data[-1]))$p.value)
RLCtools::format.pval(wilcox.test(size.data[["denovo"]], size.data[["singleton"]])$p.value)
RLCtools::format.pval(wilcox.test(size.data[["denovo"]], size.data[["rare"]])$p.value)
RLCtools::format.pval(wilcox.test(size.data[["denovo"]], size.data[["common"]])$p.value)

# Comparison of SV frequencies vs. genic context (coding, intronic, intergenic)
freq.by.context.counts <- get.freq.by.context(freq.idxs, bed)
RLCtools::format.pval(chisq.test(t(freq.by.context.counts[c("denovo", "singleton"), ]))$p.value)
RLCtools::format.pval(chisq.test(t(freq.by.context.counts[c("denovo", "rare"), ]))$p.value)
RLCtools::format.pval(chisq.test(t(freq.by.context.counts[c("denovo", "common"), ]))$p.value)
pdf(paste(args$out_prefix, "de_novo.contexts.pdf", sep="."),
    height=1.9, width=2.4)
plot.context.by.freq(freq.by.context.counts,
                     group.names=c("", "Singleton", "Rare", "Common"),
                     parmar=c(0, 3.75, 2, 0.5))
axis(2, at=3.5, labels=bquote(italic("De novo")), tick=F, las=2, line=-0.8)
dev.off()



