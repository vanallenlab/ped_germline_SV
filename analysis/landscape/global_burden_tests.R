#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Conduct genome-wide burden tests for all cancer types in one (or more) cohorts


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(Hmisc, quietly=TRUE)
require(PedSV, quietly=TRUE)
require(stats, quietly=TRUE)
require(survival, quietly=TRUE)
PedSV::load.constants("all")


##################
# Data functions #
##################
# Load precomputed burden stats from .tsv
load.precomputed.burden.stats <- function(tsv.in, keep.samples=NULL){
  df <- read.table(tsv.in, sep="\t", header=T, comment.char="")
  rnames <- df[, 1]
  df <- apply(df[, -1], 2, as.numeric)
  df[which(is.na(df))] <- 0
  rownames(df) <- rnames
  if(!is.null(keep.samples)){
    df <- df[intersect(rownames(df), keep.samples), ]
  }
  as.data.frame(df)
}


# Wrapper to format AD values as a single vector for a given query
# Note: will return list of data frames, one per cohort, if action=="verbose"
get.ad.values <- function(data, query, action, af.fields, ac.fields,
                          keep.idx.list=NULL, autosomal=TRUE){
  query.parts <- unlist(strsplit(query, split=".", fixed=T))
  res <- lapply(1:length(data), function(i){
    info <- data[[i]]
    query.bed <- filter.bed(info$bed, query, af.fields[i], ac.fields[i],
                            keep.idx=keep.idx.list[[i]], autosomal=autosomal)
    if("genes_disrupted" %in% query.parts){
      weights <- apply(query.bed[, c("PREDICTED_LOF", "PREDICTED_COPY_GAIN",
                                     "PREDICTED_INTRAGENIC_EXON_DUP")],
                       1, function(vals){sum(sapply(vals, length), na.rm=T)})
    }else if("genomic_imbalance" %in% query.parts){
      weights <- calc.genomic.imbalance(query.bed)
      names(weights) <- rownames(query.bed)
      # Set reciprocal translocations as an arbitrarily large size for this analysis
      ctx.idx <- which(query.bed$SVTYPE == "CTX")
      weights[ctx.idx] <- 50000000
    }else{
      weights <- NULL
    }
    query.ad.from.sv.bed(info$ad, query.bed, action, weights)
  })
  if(is.data.frame(res[[1]])){
    return(res)
  }else{
    unlist(res)
  }
}


# Wrapper to execute a single burden association test
burden.test <- function(data, query, meta, ad.vals, family,
                        af.fields, ac.field, extra.terms=NULL){
  # Check if cohort term is necessary
  if(length(data) > 1){
    if(is.null(extra.terms)){
      extra.terms <- c("cohort")
    }else{
      extra.terms <- unique(c(extra.terms, "cohort"))
    }
  }

  # Transform genomic imbalance terms to log10 for more intuitive scaling
  if("genomic_imbalance" %in% unlist(strsplit(query, split="."))){
    ad.vals[which(ad.vals < 1)] <- 1
    ad.vals <- log10(ad.vals)
  }

  # Run one comparison for each cancer type with at least one sample present in data
  cancers <- setdiff(unique(meta[names(ad.vals), "disease"]), "control")
  cancers <- c("pancan", metadata.cancer.label.map[cancers[which(!is.na(cancers))]])
  res <- lapply(cancers, function(cancer){
    # Get list of eligible samples
    elig.sids <- get.eligible.samples(meta, cancer)
    case.ids <- elig.sids[["cases"]]
    control.ids <- elig.sids[["controls"]]

    # Run burden test
    X <- ad.vals[intersect(names(ad.vals), c(case.ids, control.ids))]
    Y <- get.phenotype.vector(case.ids, control.ids)
    overlapping.ids <- intersect(names(X[which(!is.na(X))]), names(Y))
    case.ids.no.na <- intersect(case.ids, overlapping.ids)
    control.ids.no.na <- intersect(control.ids, overlapping.ids)
    X <- X[overlapping.ids]
    Y <- Y[overlapping.ids]
    n.carrier.ctrl <- length(which(Y[which(X > 0 & !is.na(X))] == 0))
    n.carrier.case <- length(which(Y[which(X > 0 & !is.na(X))] == 1))
    c(cancer, length(case.ids.no.na), n.carrier.case, mean(X[case.ids], na.rm=T),
      sd(X[case.ids], na.rm=T), length(control.ids.no.na), n.carrier.ctrl,
      mean(X[control.ids], na.rm=T), sd(X[control.ids], na.rm=T),
      pedsv.glm(meta, X, Y, family=family, extra.terms=extra.terms))
  })

  # Format and return data.frame of test statistics
  res.df <- as.data.frame(do.call("rbind", res), row.names=NULL)
  colnames(res.df) <- c("disease", "n.case", "n.case.carrier", "case.mean",
                        "case.stdev", "n.control", "n.control.carrier",
                        "control.mean", "control.stdev", "coefficient",
                        "std.err", "test.statistic", "P.value", "test")
  res.df$hypothesis <- query
  res.df$model <- if(family$family=="binomial"){"logit"}else{"linear"}
  res.df[, c("hypothesis", "disease", "model", "n.case", "n.case.carrier",
             "case.mean", "case.stdev", "n.control", "n.control.carrier",
             "control.mean", "control.stdev", "coefficient", "std.err",
             "test.statistic", "P.value", "test")]
}


# Supercategory burden test for a category to prepare for rapid re-testing of subsets
supercategory.burden.test <- function(data, query, meta, action, family,
                                      af.fields, ac.fields, keep.idx.list=NULL,
                                      extra.terms=NULL, autosomal=TRUE,
                                      precomp.ad.vals=NULL, aneuploidies=TRUE){
  if(is.null(precomp.ad.vals)){
    ad.df <- get.ad.values(data, query=query, action="verbose", af.fields, ac.fields,
                           keep.idx.list=keep.idx.list, autosomal=autosomal)
    ad.df <- lapply(ad.df, function(df){df[, intersect(colnames(df), rownames(meta))]})
    ad.vals <- unlist(lapply(ad.df, compress.ad.matrix, action=action))
  }else{
    ad.df <- NULL
    ad.vals <- precomp.ad.vals
  }
  query.parts <- unlist(strsplit(query, split=".", fixed=T))
  if(length(intersect(query.parts, c("notsmall", "large", "karyotypic"))) > 0
     & length(intersect(names(sv.colors), query.parts)) == 0
     & !("balanced" %in% query.parts)
     & aneuploidies){
    if(autosomal | "MALE_only" %in% query.parts | "FEMALE_only" %in% query.parts){
      has.aneu <- intersect(rownames(meta)[which(meta$autosomal_aneuploidy)], names(ad.vals))
    }else{
      has.aneu <- intersect(rownames(meta)[which(meta$any_aneuploidy)], names(ad.vals))
    }
    if(length(has.aneu) > 0){
      for(sid in has.aneu){
        if(action == "any"){
          ad.vals[sid] <- 1
        }else if(action %in% c("count", "sum")){
          ad.vals[sid] <- if(is.na(ad.vals[sid])){1}else{ad.vals[sid] + 1}
        }else if(action == "max" & "genomic_imbalance" %in% query.parts){
          # For analyses of dosage imbalance, set aneuploidy samples as having an arbitrarily large SV (>50Mb)
          ad.vals[sid] <- 50000000
        }
      }
    }
  }
  new.stats <- burden.test(data, query=query, meta=meta, ad.vals=ad.vals,
                           family=family, extra.terms=extra.terms)
  return(list("ad.df" = ad.df, "new.stats" = new.stats))
}


# Main wrapper function to run one category of burden test for each SV type
# Note: sv.subsets must be a list of three-element lists, where each inner list
# is of the format (prefix, SV IDs to include, plot title)
main.burden.wrapper <- function(data, query, meta, action, af.fields, ac.fields,
                                sv.subsets, all.stats, out.prefix, keep.idx.list=NULL,
                                extra.terms=NULL, main.title="Values", barplot.height=2.5,
                                barplot.width=3, barplot.units=NULL,
                                custom.hypothesis=NULL, autosomal=TRUE,
                                aneuploidies=TRUE){
  # Get other parameters dictated by value of action
  family <- binomial()
  if(action == "any"){
    ci.mode <- "binomial"
  }else if(action %in% c("sum", "count")){
    ci.mode <- "normal"
  }

  # Gather verbose AD matrix for reuse later, and test all SVs before splitting by type
  supercategory.res <- supercategory.burden.test(data, query=query, meta=meta,
                                                 action=action, family=family,
                                                 af.fields=af.fields, ac.fields=ac.fields,
                                                 keep.idx.list=keep.idx.list,
                                                 extra.terms=extra.terms,
                                                 autosomal=autosomal,
                                                 aneuploidies=aneuploidies)
  new.stats <- supercategory.res[["new.stats"]]
  if(!is.null(custom.hypothesis)){
    new.stats$hypothesis <- custom.hypothesis
  }
  all.stats <- rbind(all.stats, new.stats)
  pdf(paste(out.prefix, "ALL.pdf", sep="."),
      height=barplot.height, width=barplot.width)
  barplot.by.phenotype(stats2barplotdf(new.stats, ci.mode), title=main.title,
                       top.axis.units=barplot.units)
  dev.off()
  ad.dfs <- supercategory.res[["ad.df"]]
  if(is.null(sv.subsets)){
    return(all.stats)
  }
  for(subset.info in sv.subsets){
    if(any(subset.info[[2]] %in% as.vector(unlist(sapply(ad.dfs, function(df){rownames(df)}))))){
      ad.vals <- unlist(lapply(ad.dfs, function(df){
        keep.rows <- intersect(rownames(df), subset.info[[2]])
        if(length(keep.rows) > 0){
          compress.ad.matrix(df[keep.rows, ], action=action)
        }else{
          res.df <- data.frame(matrix(0, ncol=ncol(df), nrow=1))
          colnames(res.df) <- colnames(df)
          df
        }
      }))
      new.stats <- burden.test(data, query, meta, ad.vals, family, af.fields, ac.fields, extra.terms)
      if(!is.null(custom.hypothesis)){
        new.stats$hypothesis <- paste(custom.hypothesis, subset.info[[1]], sep=".")
      }else{
        new.stats$hypothesis <- paste(new.stats$hypothesis, subset.info[[1]], sep=".")
      }
      all.stats <- rbind(all.stats, new.stats)
      pdf(paste(out.prefix, subset.info[[1]], "pdf", sep="."),
          height=barplot.height, width=barplot.width)
      barplot.by.phenotype(stats2barplotdf(new.stats, ci.mode), title=subset.info[[3]],
                           top.axis.units=barplot.units)
      dev.off()
    }
  }
  return(all.stats)
}


######################
# Plotting functions #
######################
# Custom one-off plot function for case:control CNV burden by sex & size bin
plot.cnv.size.bins.supp <- function(ss, cnv, y.lims, title=NULL, bin.space=2,
                                    sex.space=0.15, parmar=c(2, 2.5, 1, 0.5)){
  # Get universal Y axis limits for both CNV types
  ymax <- max(c(log(3), 1.15 * max(as.numeric(ss$coefficient))))
  ymin <- min(c(log(0.5), 1.15 * min(as.numeric(ss$coefficient))))
  ylims <- c(ymin, ymax)

  # Reformat plot data
  plot.df <- ss[grepl(cnv, ss$hypothesis),
                c("hypothesis", "disease", "control.mean", "control.stdev",
                  "coefficient", "std.err", "P.value")]
  plot.df[, c("coefficient", "std.err")] <- apply(plot.df[, c("coefficient", "std.err")], 2, as.numeric)
  plot.df$lower_ci <- plot.df$coefficient + (qnorm(0.025) * plot.df$std.err)
  plot.df$upper_ci <- plot.df$coefficient + (qnorm(0.975) * plot.df$std.err)
  # plot.df[, c("coefficient", "std.err", "lower_ci", "upper_ci")] <-
  #   apply(plot.df[, c("coefficient", "std.err", "lower_ci", "upper_ci")], 2, exp)

  # Get other plot parameters
  bonf.cutoff <- 0.05 / nrow(ss)
  size.bins <- unique(sapply(strsplit(plot.df$hypothesis, split=".", fixed=T), function(p){p[2]}))
  bin.labels <- c("large" = ">1 Mb",
                  "100kb_to_1Mb" = "100 kb - 1 Mb",
                  "10kb_to_100kb" = "10 - 100 kb",
                  "1kb_to_10kb" = "1 - 10 kb")
  n.bins <- length(size.bins)
  n.cancers <- length(unique(plot.df$disease))
  bin.ends <- seq(n.cancers + bin.space, (n.cancers + bin.space) * n.bins,
                  by=n.cancers + bin.space) - bin.space
  bin.starts <- bin.ends - n.cancers
  xlims <- range(c(bin.starts, bin.ends))

  # Prep plot area
  prep.plot.area(xlims, ylims, parmar=parmar, xaxs="r")
  abline(h=log(1), lty=5)
  mtext(3, line=0, text=title)

  # Add points
  sapply(1:n.bins, function(b){
    sapply(2:1, function(s){
      sex <- c("MALE", "FEMALE")[s]
      idxs <- grep(paste(size.bins[b], cnv, sex, sep="."), plot.df$hypothesis)
      sapply(1:length(idxs), function(k){
        x.k <- bin.starts[b] + (`^`(-1, s) * sex.space) + k - 1
        y.k <- plot.df[idxs[k], "coefficient"]
        lci.k <- plot.df[idxs[k], "lower_ci"]
        uci.k <- plot.df[idxs[k], "upper_ci"]
        p.k <- as.numeric(plot.df[idxs[k], "P.value"])
        cancer.k <- plot.df[idxs[k], "disease"]
        pch.k <- c(22, 21)[s]
        if(p.k <= 0.005){
          fill.k <- col.k <- cancer.palettes[[cancer.k]][["dark1"]]
          border.k <- cancer.palettes[[cancer.k]][["dark2"]]
        }else if(p.k <= 0.05){
          fill.k <- col.k <-  cancer.palettes[[cancer.k]][["light2"]]
          border.k <- cancer.palettes[[cancer.k]][["light1"]]
        }else if(p.k > 0.05){
          border.k <-  cancer.palettes[[cancer.k]][["light2"]]
          col.k <-  cancer.palettes[[cancer.k]][["light3"]]
          fill.k <- "white"
        }
        segments(x0=x.k, x1=x.k, y0=lci.k, y1=uci.k, lend="round", col=col.k)
        points(x=x.k, y=y.k, pch=pch.k, col=border.k, bg=fill.k)
      })
    })
  })

  # Add bin labels and group bars
  segments(x0=bin.starts, x1=bin.ends, y0=par("usr")[3], y1=par("usr")[3],
           col="gray70", xpd=T)
  sv.avg.labels <- sapply(size.bins, function(sb){
    avg <- as.numeric(plot.df[intersect(grep(sb, plot.df$hypothesis),
                                        which(plot.df$disease == "pancan")),
                              "control.mean"])
    as.character(round(mean(avg), max(c(1, ceiling(-log10(avg))))))
  })
  x.labels <- paste(bin.labels[size.bins], "\n(", sv.avg.labels, " / genome)", sep="")
  sapply(1:n.bins, function(x){
    axis(1, at=((bin.starts + bin.ends) / 2)[x], tick=F, line=-1.5,
         labels=x.labels[x], padj=1, cex.axis=5/6)
  })

  # Add Y axis
  clean.axis(2, at=log(2^(-6:6)), labels=c(paste(2, -6:-3, sep="^"), 2^(-2:6)),
             parse.labels=TRUE, infinite=T, title="Odds ratio")
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Case:control SV burden tests")
parser$add_argument("--bed", metavar=".tsv", type="character", action="append",
                    help=paste("SV sites .bed. Can be supplied one or more",
                               "times. Must be supplied at least once. Order",
                               "must match --ad."), required=TRUE)
parser$add_argument("--ad", metavar=".tsv", type="character", action="append",
                    help=paste("Allele dosage .bed. Can be supplied one or more",
                               "times. Must be supplied at least once. Order must",
                               "match --bed."), required=TRUE)
parser$add_argument("--metadata", metavar=".tsv", type="character",
                    help="sample metadata .tsv", required=TRUE)
parser$add_argument("--gene-lists", metavar=".tsv", type="character",
                    help="Two-column .tsv of gene lists (set name and path to gene list)")
parser$add_argument("--genomic-disorder-hits", metavar=".txt", type="character",
                    help=".txt of variant IDs to be treated as genomic disorder CNVs")
parser$add_argument("--subset-samples", metavar=".tsv", type="character",
                    help="list of samples to subset [default: use all samples]")
parser$add_argument("--exclude-variants", metavar=".txt", type="character",
                    help=paste("list of variant IDs to exclude from analysis.",
                               "Not recommended for most use-cases.",
                               "[default: use all variants]"))
parser$add_argument("--precomputed-burden-stats", metavar=".tsv", type="character",
                    help=paste(".tsv of precomputed burden stats per sample for",
                               "tests too laborious to perform in-memory in R"))
parser$add_argument("--af-field", metavar="string", type="character", action="append",
                    help=paste("Column header to use for AF-related analyses.",
                               "Can be supplied once to be applied uniformly to",
                               "all values of --bed, or can be supplied an",
                               "equal number of times as --bed to apply one",
                               "value per --bed. [default: AF]"))
parser$add_argument("--ac-field", metavar="string", type="character", action="append",
                    help=paste("Column header to use for AC-related analyses.",
                               "Same requirements as --af-field. [default: AC]"))
parser$add_argument("--out-prefix", metavar="path", type="character", required=TRUE,
                    help="path/prefix for all output files")
args <- parser$parse_args()

# # DEV:
# args <- list("bed" = c("~/scratch/PedSV.v2.5.3.case_control_cohort.analysis_samples.sites.bed.gz"),
#              "ad" = c("~/scratch/PedSV.v2.5.3.case_control_cohort.analysis_samples.allele_dosages.bed.gz"),
#              "metadata" = "~/scratch/PedSV.v2.5.3.cohort_metadata.w_control_assignments.tsv.gz",
#              "gene_lists" = NULL,
#              "genomic_disorder_hits" = NULL,
#              "subset_samples" = "~/scratch/PedSV.v2.5.3.case_control_analysis_cohort.samples.list",
#              "precomputed_burden_stats" = "~/scratch/test.m.tsv",
#              "exclude_variants" = NULL,
#              "af_field" = "POPMAX_AF",
#              "ac_field" = "AC",
#              "out_prefix" = "~/scratch/PedSV.v2.5.3.case_control.dev")
# args <- list("bed" = c("~/scratch/PedSV.v2.5.3.trio_cohort.analysis_samples.sites.bed.gz"),
#              "ad" = c("~/scratch/PedSV.v2.5.3.trio_cohort.analysis_samples.allele_dosages.bed.gz"),
#              "metadata" = "~/scratch/PedSV.v2.5.3.cohort_metadata.w_control_assignments.tsv.gz",
#              "gene_lists" = NULL,
#              "genomic_disorder_hits" = NULL,
#              "subset_samples" = "~/scratch/PedSV.v2.5.3.trio_analysis_cohort.samples.list",
#              "precomputed_burden_stats" = NULL,
#              "exclude_variants" = NULL,
#              "af_field" = "POPMAX_AF",
#              "ac_field" = "AC",
#              "out_prefix" = "~/scratch/PedSV.v2.5.3.trio.dev")
# args <- list("bed" = c("~/scratch/PedSV.v2.5.4.full_cohort.analysis_samples.sites.bed.gz"),
#              "ad" = c("~/scratch/PedSV.v2.5.4.full_cohort.analysis_samples.allele_dosages.bed.gz"),
#              "metadata" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_v2_callset_generation/v2.5.4/PedSV.v2.5.4.cohort_metadata.w_control_assignments.tsv.gz",
#              "gene_lists" = "~/scratch/PedSV.gene_lists.tsv",
#              "genomic_disorder_hits" = NULL,
#              "subset_samples" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_v2_callset_generation/v2.5.4/PedSV.v2.5.4.final_analysis_cohort.samples.list",
#              "precomputed_burden_stats" = "~/scratch/PedSV.v2.5.4.full_cohort_w_relatives.precomputed_burden_stats.tsv.gz",
#              "exclude_variants" = NULL,
#              "af_field" = "POPMAX_AF",
#              "ac_field" = "AC",
#              "out_prefix" = "~/scratch/PedSV.v2.5.4.full_cohort.dev")

# Load BEDs and pair AD paths with each
data <- lapply(1:length(args$bed), function(i){
  if(is.null(args$exclude_variants)){
    bed <- load.sv.bed(args$bed[i])
  }else{
    load.sv.bed(args$bed[i], drop.vids=args$exclude_variants)
  }
  list("bed" = bed, "ad" = args$ad[i])
})

# Load metadata
keepers <- NULL
if(!is.null(args$subset_samples)){
  keepers <- read.table(args$subset_samples, header=F)[, 1]
}
meta <- load.sample.metadata(args$metadata, keep.samples=keepers,
                             reassign.parents=FALSE, annotate.aneuploidy.bp=TRUE)

# Load precomputed burden stats (if provided) and subset to samples in meta
precomp.stats <- NULL
if(!is.null(args$precomputed_burden_stats)){
  precomp.stats <- load.precomputed.burden.stats(args$precomputed_burden_stats,
                                                 rownames(meta))
}

# Format AF/AC fields
if(is.null(args$af_field)){
  af.fields <- rep("AF", times=length(data))
}else if(length(args$af_field) == 1){
  af.fields <- rep(args$af_field, times=length(data))
}else if(length(args$af_field) != length(data)){
  stop("Incorrect --af-field specification. See --help.")
}else{
  af.fields <- args$af_field
}
if(is.null(args$ac_field)){
  ac.fields <- rep("AC", times=length(data))
}else if(length(args$ac_field) == 1){
  ac.fields <- rep(args$ac_field, times=length(data))
}else if(length(args$ac_field) != length(data)){
  stop("Incorrect --ac-field specification. See --help.")
}else{
  ac.fields <- args$ac_field
}

# Read gene lists, if optioned
if(!is.null(args$gene_lists)){
  gene.lists <- load.gene.lists(args$gene_lists)
}

# Prepare data.frame for collecting test statistics
empty.stats.df <- data.frame("hypothesis"=character(0), "disease"=character(0),
                             "model"=character(0), "n.case"=numeric(0),
                             "case.mean"=numeric(0), "n.control"=numeric(0),
                             "control.mean"=numeric(0), "coefficient"=numeric(0),
                             "std.err"=numeric(0), "test.statistic"=numeric(0),
                             "P.value"=numeric(0), "test"=character(0))
all.stats <- empty.stats.df

# Set plot dimensions
barplot.height <- 0.5 + (length(unique(meta$disease)) / 4)
barplot.width <- 3
swarmplot.height <- 2.25
swarmplot.width <- 0.75 + (length(unique(meta$disease)) / 2)


# Count of large autosomal variants per genome by SV type
sv.subsets <- lapply(c("DEL", "DUP", "INV", "CPX"), function(svtype){
  list(svtype,
       as.character(unlist(sapply(data, function(info){
         rownames(info$bed)[which(info$bed$SVTYPE == svtype)]
       }))),
       paste("Samples with", tolower(sv.abbreviations[svtype]), ">1 Mb"))
})
sv.subsets <- c(sv.subsets,
                list(list("CTX",
                          as.character(unlist(sapply(data, function(info){
                            rownames(info$bed)[which(info$bed$SVTYPE == "CTX")]
                          }))),
                          "Samples with reciprocal translocations")))
all.stats <- main.burden.wrapper(data, query="large", meta, action="any",
                                 af.fields, ac.fields, sv.subsets=sv.subsets, all.stats,
                                 paste(args$out_prefix, "large_sv_per_genome.by_cancer", sep="."),
                                 main.title="Samples with any SV >1 Mb",
                                 barplot.height=barplot.height, barplot.width=barplot.width,
                                 barplot.units="percent")


# Carrier rate of rare/vrare/singleton large (1 Mb) & notsmall (50 kb) variants per genome by SV type
for(freq in c("rare", "vrare", "singleton")){
  for(size in c("large", "notsmall")){
    sv.subsets <- lapply(c("DEL", "DUP", "INV", "CPX"), function(svtype){
      list(svtype,
           as.character(unlist(sapply(data, function(info){
             rownames(info$bed)[which(info$bed$SVTYPE == svtype)]
           }))),
           paste("Samples with", tolower(freq.names[freq]), tolower(sv.abbreviations[svtype]),
                 if(svtype != "CTX"){if(size == "large"){">1 Mb"}else{">100 kb"}}))
    })
    all.stats <- main.burden.wrapper(data, query=paste(size, freq, sep="."), meta, action="any",
                                     af.fields, ac.fields, sv.subsets=sv.subsets, all.stats,
                                     paste(args$out_prefix, paste(freq, size, "sv_per_genome.by_cancer", sep="_"), sep="."),
                                     main.title=paste("Samples with", tolower(freq.names[freq]), "SV",
                                                      if(size == "large"){">1 Mb"}else{">100 kb"}),
                                     barplot.height=barplot.height, barplot.width=barplot.width,
                                     barplot.units="percent", autosomal=FALSE)
    # Add one separate test for all large *unbalanced* SVs
    # (DEL + DUP + aneuploidy + qualifying CPX)
    all.stats <- main.burden.wrapper(data, query=paste(size, freq, "unbalanced", sep="."), meta, action="any",
                                     af.fields, ac.fields, sv.subsets=NULL, all.stats,
                                     paste(args$out_prefix, paste(freq, size, "unbalanced_sv_per_genome.by_cancer", sep="_"), sep="."),
                                     main.title=paste("Pct. w/", tolower(freq.names[freq]),
                                                      "unbalanced SV",
                                                      if(size == "large"){">1 Mb"}else{">100 kb"}),
                                     barplot.height=barplot.height, barplot.width=barplot.width,
                                     barplot.units="percent", autosomal=FALSE)
    # Add one separate test for all large *unbalanced* SVs for each sex
    # (DEL + DUP + aneuploidy + qualifying CPX)
    for(sex in c("MALE", "FEMALE")){
      all.stats <- main.burden.wrapper(data, query=paste(size, freq, "unbalanced",
                                                         paste(sex, "only", sep="_"), sep="."),
                                       meta[which(meta$inferred_sex == sex), ], action="any",
                                       af.fields, ac.fields, sv.subsets=NULL, all.stats,
                                       paste(args$out_prefix, paste(freq, size, "unbalanced_sv_per_genome.by_cancer", sep="_"),
                                             paste(sex, "only", sep="_"), sep="."),
                                       main.title=paste("Pct. w/", tolower(freq.names[freq]),
                                                        "unbalanced SV",
                                                        if(size == "large"){">1 Mb"}else{">100 kb"}),
                                       barplot.height=barplot.height, barplot.width=barplot.width,
                                       barplot.units="percent", autosomal=FALSE)
    }
    # Add one separate test for all large *unbalanced* *autosomal* SVs
    # (DEL + DUP + aneuploidy + qualifying CPX)
    all.stats <- main.burden.wrapper(data, query=paste(size, freq, "unbalanced", sep="."), meta, action="any",
                                     af.fields, ac.fields, sv.subsets=NULL, all.stats,
                                     paste(args$out_prefix, paste(freq, size, "unbalanced_sv_per_genome.by_cancer.autosomal_only", sep="_"), sep="."),
                                     main.title=paste("Pct. w/", tolower(freq.names[freq]),
                                                      "unbalanced SV",
                                                      if(size == "large"){">1 Mb"}else{">100 kb"}),
                                     barplot.height=barplot.height, barplot.width=barplot.width,
                                     barplot.units="percent", autosomal=TRUE,
                                     custom.hypothesis=paste(size, freq, "unbalanced", "autosomal_only", sep="."))
    # Add one separate test for all large *unbalanced* *autosomal* SVs for each sex
    # (DEL + DUP + aneuploidy + qualifying CPX)
    for(sex in c("MALE", "FEMALE")){
      all.stats <- main.burden.wrapper(data, query=paste(size, freq, "unbalanced",
                                                         paste(sex, "only", sep="_"), sep="."),
                                       meta[which(meta$inferred_sex == sex), ], action="any",
                                       af.fields, ac.fields, sv.subsets=NULL, all.stats,
                                       paste(args$out_prefix, paste(freq, size, "unbalanced_sv_per_genome.by_cancer.autosomal_only", sep="_"),
                                             paste(sex, "only", sep="_"), sep="."),
                                       main.title=paste("Pct. w/", tolower(freq.names[freq]),
                                                        "unbalanced SV",
                                                        if(size == "large"){">1 Mb"}else{">100 kb"}),
                                       barplot.height=barplot.height, barplot.width=barplot.width,
                                       barplot.units="percent", autosomal=TRUE,
                                       custom.hypothesis=paste(size, freq, "unbalanced", "autosomal_only",
                                                               paste(sex, "only", sep="_"), sep="."))
    }
    # Add one separate test for all large *unbalanced* *allosomal* SVs
    # (DEL + DUP + aneuploidy + qualifying CPX)
    all.stats <- main.burden.wrapper(data, query=paste(size, freq, "unbalanced", sep="."), meta, action="any",
                                     af.fields, ac.fields, sv.subsets=NULL, all.stats,
                                     keep.idx.list=lapply(data, function(d){which(d$bed$chrom %in% c("chrX", "chrY"))}),
                                     paste(args$out_prefix, paste(freq, size, "unbalanced_sv_per_genome.by_cancer.allosomal_only", sep="_"), sep="."),
                                     main.title=paste("Pct. w/", tolower(freq.names[freq]),
                                                      "unbalanced SV",
                                                      if(size == "large"){">1 Mb"}else{">100 kb"}),
                                     barplot.height=barplot.height, barplot.width=barplot.width,
                                     barplot.units="percent", autosomal=FALSE,
                                     custom.hypothesis=paste(size, freq, "unbalanced", "allosomal_only", sep="."))
    # Add one separate test for all large *balanced* SVs
    # (INS + INV + qualifying CPX)
    all.stats <- main.burden.wrapper(data, query=paste(size, freq, "balanced", sep="."), meta, action="any",
                                     af.fields, ac.fields, sv.subsets=NULL, all.stats,
                                     paste(args$out_prefix, paste(freq, size, "balanced_sv_per_genome.by_cancer", sep="_"), sep="."),
                                     main.title=paste("Pct. w/", tolower(freq.names[freq]),
                                                      "balanced SV",
                                                      if(size == "large"){">1 Mb"}else{">100 kb"}),
                                     barplot.height=barplot.height, barplot.width=barplot.width,
                                     barplot.units="percent", autosomal=FALSE)
  }
  # Add one separate test for all huge (>5Mb) *unbalanced* SVs
  # (DEL + DUP + aneuploidy + qualifying CPX)
  all.stats <- main.burden.wrapper(data, query=paste("karyotypic", freq, "unbalanced", sep="."), meta, action="any",
                                   af.fields, ac.fields, sv.subsets=NULL, all.stats,
                                   paste(args$out_prefix, paste(freq, "karyotypic_unbalanced_sv_per_genome.by_cancer", sep="_"), sep="."),
                                   main.title=paste("Pct. w/", tolower(freq.names[freq]),
                                                    "unbalanced SV >5Mb"),
                                   barplot.height=barplot.height, barplot.width=barplot.width,
                                   barplot.units="percent", autosomal=FALSE)
}


# Rerun selected tests of large unbalanced SVs after restricting to Europeans only
for(size in c("large")){
  for(freq in c("rare")){
    # All large *unbalanced* *autosomal* SVs
    # (DEL + DUP + aneuploidy + qualifying CPX)
    all.stats <- main.burden.wrapper(data, query=paste(size, freq, "unbalanced", sep="."),
                                     meta[which(meta$inferred_ancestry == "EUR"), ], action="any",
                                     af.fields, ac.fields, sv.subsets=NULL, all.stats,
                                     paste(args$out_prefix, paste(freq, size, "unbalanced_sv_per_genome.by_cancer.autosomal_only.EUR_only", sep="_"), sep="."),
                                     main.title=paste("Pct. w/", tolower(freq.names[freq]),
                                                      "unbalanced SV",
                                                      if(size == "large"){">1 Mb"}else{">100 kb"}),
                                     barplot.height=barplot.height, barplot.width=barplot.width,
                                     barplot.units="percent", autosomal=TRUE,
                                     custom.hypothesis=paste(size, freq, "unbalanced", "autosomal_only", "EUR_only", sep="."))

    # Add one separate test for all large *unbalanced* *autosomal* SVs for each sex
    # (DEL + DUP + qualifying CPX; *no* aneuploidy)
    for(sex in c("MALE", "FEMALE")){
      all.stats <- main.burden.wrapper(data, query=paste(size, freq, "unbalanced",
                                                         paste(sex, "only", sep="_"), sep="."),
                                       meta[which(meta$inferred_sex == sex & meta$inferred_ancestry == "EUR"), ], action="any",
                                       af.fields, ac.fields, sv.subsets=NULL, all.stats,
                                       paste(args$out_prefix, paste(freq, size, "unbalanced_sv_per_genome.by_cancer.autosomal_only.EUR_only", sep="_"),
                                             paste(sex, "only", sep="_"), sep="."),
                                       main.title=paste("Pct. w/", tolower(freq.names[freq]),
                                                        "unbalanced SV",
                                                        if(size == "large"){">1 Mb"}else{">100 kb"}),
                                       barplot.height=barplot.height, barplot.width=barplot.width,
                                       barplot.units="percent", autosomal=TRUE,
                                       custom.hypothesis=paste(size, freq, "unbalanced", "autosomal_only", "EUR_only",
                                                               paste(sex, "only", sep="_"), sep="."))
    }
  }
}


# Rerun selected tests of large unbalanced SVs after excluding COSMIC/CPG
if(!is.null(args$gene_lists)){
  if(length(intersect(names(gene.lists), c("COSMIC", "CPG"))) == 2){

    # Exclude all SVs with a predicted coding consequence on any COSMIC gene or CPG
    x.genes <- sort(unique(unlist(gene.lists[c("COSMIC", "CPG")])))
    x.csqs <- c("PREDICTED_COPY_GAIN", "PREDICTED_DUP_PARTIAL", "PREDICTED_INTRAGENIC_EXON_DUP",
                "PREDICTED_LOF", "PREDICTED_PARTIAL_EXON_DUP")
    no.xgenes <- lapply(data, function(l){
      which(apply(l$bed[, x.csqs], 1, function(gl){
        length(intersect(unlist(gl), x.genes)) == 0
      }))
    })

    # Set test parameters
    size <- "large"
    freq <- "rare"

    # All large *unbalanced* *autosomal* SVs
    # (DEL + DUP + qualifying CPX; *no* aneuploidy)
    all.stats <- main.burden.wrapper(data, query=paste(size, freq, "unbalanced", sep="."),
                                     meta, action="any",
                                     af.fields, ac.fields, sv.subsets=NULL, all.stats,
                                     keep.idx.list=no.xgenes,
                                     paste(args$out_prefix, paste(freq, size, "unbalanced_sv_per_genome.by_cancer.autosomal_only.no_COSMIC_no_CPG", sep="_"), sep="."),
                                     main.title=paste("Pct. w/", tolower(freq.names[freq]),
                                                      "unbalanced SV",
                                                      if(size == "large"){">1 Mb"}else{">100 kb"}),
                                     barplot.height=barplot.height, barplot.width=barplot.width,
                                     barplot.units="percent", autosomal=TRUE, aneuploidies=FALSE,
                                     custom.hypothesis=paste(size, freq, "unbalanced", "autosomal_only", "no_COSMIC_no_CPG", sep="."))

    # Add one separate test for all large *unbalanced* *autosomal* SVs for each sex
    # (DEL + DUP + qualifying CPX; *no* aneuploidy)
    for(sex in c("MALE", "FEMALE")){
      all.stats <- main.burden.wrapper(data, query=paste(size, freq, "unbalanced",
                                                         paste(sex, "only", sep="_"), sep="."),
                                       meta[which(meta$inferred_sex == sex), ], action="any",
                                       af.fields, ac.fields, sv.subsets=NULL, all.stats,
                                       keep.idx.list=no.xgenes,
                                       paste(args$out_prefix, paste(freq, size, "unbalanced_sv_per_genome.by_cancer.autosomal_only.no_COSMIC_no_CPG", sep="_"),
                                             paste(sex, "only", sep="_"), sep="."),
                                       main.title=paste("Pct. w/", tolower(freq.names[freq]),
                                                        "unbalanced SV",
                                                        if(size == "large"){">1 Mb"}else{">100 kb"}),
                                       barplot.height=barplot.height, barplot.width=barplot.width,
                                       barplot.units="percent", autosomal=TRUE, aneuploidies=FALSE,
                                       custom.hypothesis=paste(size, freq, "unbalanced", "autosomal_only", "no_COSMIC_no_CPG",
                                                               paste(sex, "only", sep="_"), sep="."))
    }
  }
}


# Focused secondary analyses of largest rare unbalanced SV per genome by sex
for(freq in c("rare", "vrare", "singleton")){
  smooth.or.collection <- list("x" = list(), "or" = list(),
                               "or.lower" = list(), "or.upper" = list())
  for(sex in c("ALL", "MALE", "FEMALE")){
    # Set sex-specific values and parameters
    if(sex == "ALL"){
      sex.meta <- meta
      sex.title <- "all_sexes"
      km.y.title <- "Samples (%)"
    }else{
      sex.meta <- meta[which(meta$inferred_sex == sex), ]
      sex.title <- paste(sex, "only", sep="_")
      km.y.title <- paste(sex.names[sex], "samples (%)")
    }

    # Collect all SV data
    largest.sv.data <-
      supercategory.burden.test(data=data,
                                query=paste(freq, "unbalanced.notsmall.genomic_imbalance",
                                            sex.title, sep="."),
                                meta=sex.meta, family=binomial(), action="max",
                                af.fields=af.fields, ac.fields=ac.fields)
    all.stats <- rbind(all.stats, largest.sv.data$new.stats)

    # K-M style visualization of largest unbalanced SV vs. case status
    largest.sv <- apply(largest.sv.data$ad.df[[1]], 2, max, na.rm=T)
    cancers.km.layer.ordered <- c("control", "pancan")
    surv.models <- lapply(cancers.km.layer.ordered, function(cancer){
      elig.sids <- intersect(get.eligible.samples(sex.meta, cancer)$cases, names(largest.sv))
      survfit(Surv(log10(largest.sv[elig.sids]), rep(1, length(elig.sids))) ~ 1)
    })
    ylims <- c(0, max(sapply(surv.models, function(ss){
      max(ss$surv[which(ss$time>log10(50000))], na.rm=T)
    }), na.rm=T) + 0.025)

    pdf(paste(args$out_prefix, freq, "genomic_imbalance_km", sex.title, "pdf", sep="."),
        height=2.15, width=2.5)
    svlen.line.plot(x.svlen=lapply(surv.models, function(sm){sm$time}),
                    y.value=lapply(surv.models, function(sm){sm$surv}),
                    ci.lower=lapply(surv.models, function(sm){sm$lower}),
                    ci.upper=lapply(surv.models, function(sm){sm$upper}),
                    colors=cancer.colors[cancers.km.layer.ordered], ci.alpha=c(0, 0),
                    xlab=paste("Largest unbalanced\n", tolower(freq.names[freq]), "SV"),
                    xlab.line=1, ylab=km.y.title, xlim=log10(c(50000, 5000000)),
                    lwds=c(3, 3), x.axis.labels=c("50 kb", "500 kb", "5 Mb"),
                    ylims=ylims, y.axis.units="percent",
                    x.axis.labels.at=log10(c(50000, 500000, 5000000)),
                    parmar=c(3, 3, 0.25, 1))
    rect(xleft=log10(500000), xright=log10(5000000), ybottom=0, ytop=0.05,
         col=NA, lty=2, xpd=T)
    dev.off()

    # Small inset plot starting at 500 kb
    pdf(paste(args$out_prefix, freq, "genomic_imbalance_km.inset", sex.title,
              "pdf", sep="."),
        height=0.6*(12/5), width=0.65*(12/5))
    svlen.line.plot(x.svlen=lapply(surv.models, function(sm){sm$time}),
                    y.value=lapply(surv.models, function(sm){sm$surv}),
                    ci.lower=lapply(surv.models, function(sm){sm$lower}),
                    ci.upper=lapply(surv.models, function(sm){sm$upper}),
                    colors=cancer.colors[cancers.km.layer.ordered], ci.alpha=c(0, 0),
                    xlab=NA, ylab=NA, xlim=log10(c(500000, 5000000)),
                    ylims=c(0, 0.05), y.axis.units="percent",
                    lwds=c(rep(2, length(cancers.km.layer.ordered)-2), 3, 3),
                    x.axis.labels=c("1 Mb", "5Mb"),
                    x.axis.labels.at=log10(c(1000000, 5000000)), x.tck=-0.025,
                    parmar=c(1.15, 1.75, 0.4, 0.75))
    dev.off()

    # Burden test series every log-step
    tryCatch(
      {
        ad.df <- get.ad.values(data, query=paste(freq, "notsmall", "unbalanced", sex.title, sep="."),
                               action="verbose", af.fields=af.fields, ac.fields=ac.fields,
                               autosomal=FALSE)
        unbal.sv.size <- unlist(lapply(data, function(l){v <- calc.genomic.imbalance(l$bed); names(v) <- rownames(l$bed); return(v)}))
        size.burden.stats <- empty.stats.df
        size.or.x <- 10^seq(4, log10(5000000), 0.05)
        for(size in size.or.x){
          qual.sv.ids <- names(unbal.sv.size)[which(unbal.sv.size >= size)]
          ad.vals <- unlist(lapply(ad.df, compress.ad.matrix, action="any", keep.vids=qual.sv.ids))
          ad.vals[intersect(rownames(sex.meta)[which(sex.meta$any_aneuploidy)], names(ad.vals))] <- 1
          new.burden.stats <- burden.test(data, query=paste(freq, ".", size/1000, "kb", sep=""),
                                          meta=sex.meta, ad.vals=ad.vals, family=binomial(),
                                          af.fields=af.fields, ac.field=ac.fields)
          size.burden.stats <- rbind(size.burden.stats, new.burden.stats)
        }
        size.burden.y.ln <- as.numeric(size.burden.stats[which(size.burden.stats$disease == "pancan"), "coefficient"])
        size.burden.y.ln.se <- as.numeric(size.burden.stats[which(size.burden.stats$disease == "pancan"), "std.err"])
        size.burden.y <- exp(size.burden.y.ln)
        size.burden.y.smooth <- predict(smooth.spline(size.or.x, size.burden.y, spar=2), newdata=size.or.x)$y
        size.burden.ci.lower <- exp(size.burden.y.ln + (qnorm(0.025) * size.burden.y.ln.se))
        size.burden.ci.lower.smooth <- predict(smooth.spline(size.or.x, size.burden.ci.lower, spar=2), newdata=size.or.x)$y
        size.burden.ci.upper <- exp(size.burden.y.ln + (qnorm(0.975) * size.burden.y.ln.se))
        size.burden.ci.upper.smooth <- predict(smooth.spline(size.or.x, size.burden.ci.upper, spar=2), newdata=size.or.x)$y
        pdf(paste(args$out_prefix, freq, "genomic_imbalance_effect_size", sex.title, "pdf", sep="."),
            height=2.15, width=2)
        svlen.line.plot(x.svlen=list(log10(size.or.x)),
                        y.value=list(size.burden.y.smooth),
                        ci.lower=list(size.burden.ci.lower.smooth),
                        ci.upper=list(size.burden.ci.upper.smooth),
                        step=FALSE, colors=cancer.colors["pancan"], ci.alpha=0.25,
                        xlab=paste("Largest unbalanced\n", tolower(freq.names[freq]), "SV"),
                        xlab.line=1, ylab="Odds ratio", y.title.line=0.1,
                        xlim=log10(c(50000, 5000000)),
                        lwds=c(3, 3), x.axis.labels=c("50 kb", "500 kb", "5 Mb"),
                        ylims=c(min(c(0.8, min(size.burden.y.smooth))),
                                max(c(1.2, 1.1*size.burden.y.smooth))),
                        x.axis.labels.at=log10(c(50000, 500000, 5000000)),
                        parmar=c(3, 2, 0.25, 1))
        dev.off()

        # Add smoothed effect sizes to collector
        smooth.or.collection[["x"]][[sex]] <- log10(size.or.x)
        smooth.or.collection[["or"]][[sex]] <- size.burden.y.smooth
        smooth.or.collection[["or.lower"]][[sex]] <- size.burden.ci.lower.smooth
        smooth.or.collection[["or.upper"]][[sex]] <- size.burden.ci.upper.smooth
      },
      error=function(e){}
    )
  }

  # Plot effect size visualization with male/female split overlaid
  tryCatch(
    {
      pdf(paste(args$out_prefix, freq, "genomic_imbalance_effect_size",
                "sex_overlay", "pdf", sep="."),
          height=2.15, width=2)
      svlen.line.plot(x.svlen=rev(smooth.or.collection$x),
                      y.value=rev(smooth.or.collection$or),
                      ci.lower=rev(smooth.or.collection$or.lower),
                      ci.upper=rev(smooth.or.collection$or.upper),
                      step=FALSE, ci.alpha=rev(c(0.1, 0, 0)),
                      colors=rev(c(cancer.colors["pancan"], sex.colors["MALE"], sex.colors["FEMALE"])),
                      xlab=paste("Largest unbalanced\n", tolower(freq.names[freq]), "SV"),
                      xlab.line=1, ylab="Odds ratio", y.title.line=0.1,
                      xlim=log10(c(50000, 5000000)),
                      lwds=c(3, 3, 3), x.axis.labels=c("50 kb", "500 kb", "5 Mb"),
                      ylims=c(min(c(0.8, unlist(smooth.or.collection$or))),
                              max(c(1.2, 1.1*unlist(smooth.or.collection$or)))),
                      x.axis.labels.at=log10(c(50000, 500000, 5000000)),
                      parmar=c(3, 2, 0.25, 1))
      dev.off()
    },
    error=function(e){}
  )
}


# Focused analysis of singleton unbalanced SV burden by size range, sex, and SV type
# Collect all data
cnv.to.unbal.query.map <- c("DEL" = "copy_loss", "DUP" = "copy_gain")
cnv.size.bins <- c("large", "100kb_to_1Mb", "10kb_to_100kb", "1kb_to_10kb")
ss <- do.call("rbind", lapply(c("DEL", "DUP"), function(cnv){
  # Collect data for copy gain & copy loss separately
  do.call("rbind", lapply(cnv.size.bins, function(size){
    query <- paste("singleton", size, cnv.to.unbal.query.map[cnv], sep=".")
    ad.vals <- get.ad.values(data, query, action="count", af.fields=af.fields,
                             ac.fields=ac.fields, autosomal=TRUE)
    do.call("rbind", lapply(c("MALE", "FEMALE"), function(sex){
      burden.test(data, paste(query, ".", sex, "_only", sep=""),
                  meta[which(meta$inferred_sex == sex), ],
                  ad.vals, family=binomial(), af.fields=af.fields,
                  ac.field=ac.fields)
    }))
  }))
}))
all.stats <- rbind(all.stats, ss)
# Plot data for each of gain/loss in separate panel, but with matching Y axis
for(cnv in c("DEL", "DUP")){
  title <- paste("Autosomal singleton ", tolower(sv.names[cnv]), "s", sep="")
  pdf(paste(args$out_prefix, "singleton", cnv, "burden_by_size_range_and_sex",
            "pdf", sep="."), height=2.15, width=4.6)
  plot.cnv.size.bins.supp(ss, cnv.to.unbal.query.map[cnv], title=title)
  dev.off()
}


# Sum total of autosomal genomic imbalance per genome by SV type
if(!is.null(precomp.stats)){
  for(freq in c("rare", "vrare", "singleton")){
    for(cnv in c("DEL", "DUP", "CNV")){
      precomp.colname <- paste(freq, "imbalance", cnv, sep="_")
      if(precomp.colname %in% colnames(precomp.stats)){

        # Get precomputed values, which do not account for aneuploidies
        ad.vals.presex <- as.vector(precomp.stats[, precomp.colname])
        names(ad.vals.presex) <- rownames(precomp.stats)

        # Run association test (also split by sex)
        for(sex in c("ALL", "MALE", "FEMALE")){
          # Prepare sex-specific data and parameters
          if(sex == "ALL"){
            # Add aneuploidies only for model with both sexes
            auto.aneu.bp <- meta[names(ad.vals.presex), "autosomal_aneuploidy_bp"]
            sex.aneu.bp <- meta[names(ad.vals.presex), "sex_aneuploidy_bp"]
            if(cnv == "DEL"){
              auto.aneu.bp <- sapply(auto.aneu.bp, function(v){min(c(0, v))})
              sex.aneu.bp <- sapply(sex.aneu.bp, function(v){min(c(0, v))})
            }else if(cnv == "DUP"){
              auto.aneu.bp <- sapply(auto.aneu.bp, function(v){max(c(0, v))})
              sex.aneu.bp <- sapply(sex.aneu.bp, function(v){max(c(0, v))})
            }
            ad.vals <- ad.vals.presex + abs(auto.aneu.bp) + abs(sex.aneu.bp)
            sex.meta <- meta
            sex.title <- "all_sexes"

          }else{
            ad.vals <- ad.vals.presex
            sex.meta <- meta[which(meta$inferred_sex == sex), ]
            sex.title <- paste(sex, "only", sep="_")

          }
          ad.vals <- log10(ad.vals)
          ad.vals[which(is.infinite(ad.vals))] <- 0

          new.stats <- supercategory.burden.test(data, query=paste(paste("total", freq, cnv, "imbalance_per_genome", sep="_"),
                                                                   sex.title, sep="."),
                                                 sex.meta, action="sum", family=binomial(),
                                                 af.fields, ac.fields,
                                                 precomp.ad.vals=ad.vals)$new.stats
          all.stats <- rbind(all.stats, new.stats)

          # Plot swarms
          ylims <- c(min(ad.vals[which(!is.infinite(ad.vals))]), log10(5000000))
          sv.label <- if(cnv == "CNV"){"CNVs"}else{tolower(sv.abbreviations[cnv])}
          pdf(paste(args$out_prefix,
                    paste("total", freq, cnv, "imbalance_per_genome", sep="_"),
                    sex.title, "pdf", sep="."),
              height=swarmplot.height, width=swarmplot.width)
          swarmplot.by.phenotype(ad.vals, sex.meta, ylims=ylims,
                                 title=format.pval(as.numeric(new.stats$P.value)[1]),
                                 y.axis.title=paste("Nucleotides altered\nby",
                                                    tolower(freq.names[freq]),
                                                    sv.label),
                                 title.line=0.5, y.title.line=1.3,
                                 y.ticks=log10(logscale.major.bp),
                                 y.tick.labels=logscale.major.bp.labels,
                                 parse.labels=FALSE,
                                 parmar=c(3, 4.5, 1.5, 0.25))
          axis(2, at=log10(logscale.minor), labels=NA, tck=-0.0125)
          dev.off()
        }
      }
    }
  }
}


# Genomic disorder analysis, if optioned
if(!is.null(args$genomic_disorder_hits)){
  gd.ids <- unique(read.table(args$genomic_disorder_hits, header=F)[, 1])
  if(length(gd.ids) > 0){
    for(freq in c("rare", "vrare")){
      sv.subsets <- lapply(c("DEL", "DUP"), function(svtype){
        list(svtype,
             as.character(unlist(sapply(data, function(info){
               intersect(rownames(info$bed)[which(info$bed$SVTYPE == svtype)],
                         gd.ids)
             }))),
             paste("Samples w/", tolower(freq.names[freq]), "GD",
                   tolower(sv.abbreviations[svtype])))
      })
      keep.idx.list <- lapply(data, function(d){which(rownames(d$bed) %in% gd.ids)})
      all.stats <- main.burden.wrapper(data, query=freq, meta, action="any",
                                       af.fields, ac.fields, sv.subsets=sv.subsets, all.stats,
                                       out.prefix=paste(args$out_prefix, paste(freq, "genomic_disorders.by_cancer", sep="_"), sep="."),
                                       keep.idx.list=keep.idx.list,
                                       main.title=paste("Samples with", tolower(freq.names[freq]),
                                                        "Genomic Disorder"),
                                       barplot.height=barplot.height, barplot.width=barplot.width,
                                       barplot.units="percent",
                                       custom.hypothesis=paste("genomic_disorders", freq, sep="."))
    }
  }
}


# Total number of rare/vrare/singleton variants (of any size) per genome
if(!is.null(precomp.stats)){
  for(freq in c("rare", "vrare", "singleton")){
    precomp.colname <- paste(freq, "sv_count", sep="_")
    if(precomp.colname %in% colnames(precomp.stats)){
      ad.vals <- as.vector(precomp.stats[, precomp.colname])
      names(ad.vals) <- rownames(precomp.stats)

      # Perform this analysis within continental ancestry groups to
      # more explicitly control for pop strat
      cases.per.pop <- table(meta$inferred_ancestry[which(meta$disease != "control")])
      sapply(names(cases.per.pop)[which(cases.per.pop > 50)], function(pop){
        pop.meta <- meta[which(meta$inferred_ancestry == pop), ]
        pop.ad.vals <- ad.vals[intersect(names(ad.vals),
                                         rownames(meta)[which(meta$inferred_ancestry == pop)])]
        ylims <- round(quantile(pop.ad.vals, probs=c(0.005, 0.995), na.rm=T), 0)
        new.stats <- supercategory.burden.test(data, query=paste("all", freq, "per_genome", pop, "only", sep="_"),
                                               pop.meta, action="count", family=binomial(),
                                               af.fields, ac.fields,
                                               precomp.ad.vals=pop.ad.vals)$new.stats
        all.stats <- rbind(all.stats, new.stats)
        pdf(paste(args$out_prefix, paste("all", freq, "per_genome", pop, "only", sep="_"), "pdf", sep="."),
            height=swarmplot.height, width=swarmplot.width)
        swarmplot.by.phenotype(pop.ad.vals, pop.meta, ylims=ylims, title.line=0.5,
                               title=format.pval(as.numeric(new.stats$P.value)[1]),
                               y.axis.title=paste(freq.names[freq], "SVs"),
                               parmar=c(3, 4.5, 1.5, 0.25))
        dev.off()
      })
    }
  }
}


# Total number of all autosomal gene-disruptive SVs per genome, rendered as swarmplot
ad.vals <- get.ad.values(data, query="genes_disrupted", action="sum",
                         af.fields, ac.fields, autosomal=TRUE)
new.stats <- burden.test(data, "all_genes_disrupted",
                         meta, ad.vals, binomial(), af.fields, ac.fields)
all.stats <- rbind(all.stats, new.stats)
pdf(paste(args$out_prefix, "all_genes_disrupted_per_genome.by_cancer.pdf", sep="."),
    height=swarmplot.height+0.1, width=3)
swarmplot.by.phenotype(ad.vals, meta, ylims=quantile(ad.vals, probs=c(0, 0.999), na.rm=T), title.line=0.5,
                       title=format.pval(as.numeric(new.stats$P.value)[1]),
                       y.axis.title="Genes disrupted by all SVs", y.title.line=0.75,
                       y.ticks=c(0, 10e10),
                       parmar=c(3, 2.75, 1.6, 0.25))
clean.axis(2)
dev.off()



# Number of rare/vrare/singleton LoF, CG, IEDs per genome
for(freq in c("rare", "vrare", "singleton")){
  sv.subsets <- list(
    list(paste(freq, "LoF_SVs", sep="_"),
         as.character(unlist(sapply(data, function(info){
           rownames(info$bed)[which(sapply(info$bed$PREDICTED_LOF, length) > 0)]
         }))),
         paste(freq.names[freq], "LoF SVs per sample")),
    list(paste(freq, "LoF_DELs", sep="_"),
         as.character(unlist(sapply(data, function(info){
           rownames(info$bed)[which(sapply(info$bed$PREDICTED_LOF, length) > 0
                                    & info$bed$SVTYPE == "DEL")]
         }))),
         paste(freq.names[freq], "LoF dels. per sample")),
    list(paste(freq, "LoF_nonDELs", sep="_"),
         as.character(unlist(sapply(data, function(info){
           rownames(info$bed)[which(sapply(info$bed$PREDICTED_LOF, length) > 0
                                    & info$bed$SVTYPE != "DEL")]
         }))),
         paste(freq.names[freq], "non-del. LoF SVs")),
    list(paste(freq, "CG_SVs", sep="_"),
         as.character(unlist(sapply(data, function(info){
           rownames(info$bed)[which(sapply(info$bed$PREDICTED_COPY_GAIN, length) > 0)]
         }))),
         paste(freq.names[freq], "CG SVs per sample")),
    list(paste(freq, "IED_SVs", sep="_"),
         as.character(unlist(sapply(data, function(info){
           rownames(info$bed)[which(sapply(info$bed$PREDICTED_INTRAGENIC_EXON_DUP, length) > 0)]
         }))),
         paste(freq.names[freq], "IED SVs per sample")),
    list(paste(freq, "nonLoF_disruptive_SVs", sep="_"),
         as.character(unlist(sapply(data, function(info){
           rownames(info$bed)[which(sapply(info$bed$PREDICTED_LOF, length) == 0)]
         }))),
         paste(freq.names[freq], "non-LoF gene-disruptive SVs"))
  )
  for(sv.size.limit in c("", "no_large_unbalanced")){
    if(sv.size.limit == "no_large_unbalanced"){
      test.data <- lapply(data, function(l){
        list("bed" = l$bed[which(calc.genomic.imbalance(l$bed) < 1000000), ],
             "ad" = l[[2]])
      })
    }else{
      test.data <- data
    }
    query <- gsub("\\.$", "", gsub("..", ".", paste(freq, "gene_disruptive", sv.size.limit, sep="."), fixed=T))
    all.stats <- main.burden.wrapper(test.data, query=query, meta, action="count",
                                     af.fields, ac.fields, sv.subsets=sv.subsets, all.stats,
                                     paste(args$out_prefix, query, "by_cancer", sep="."),
                                     main.title=paste(freq.names[freq], "gene-disruptive SVs / sample"),
                                     barplot.height=barplot.height, barplot.width=barplot.width)
    # One for each stratified by sex, restricting to *autosomal* only
    for(sex in c("MALE", "FEMALE")){
      all.stats <- main.burden.wrapper(test.data, query=paste(query, paste(sex, "only", sep="_"), sep="."),
                                       meta[which(meta$inferred_sex == sex), ], action="count",
                                       af.fields, ac.fields, sv.subsets=sv.subsets, all.stats,
                                       paste(args$out_prefix, paste(freq, "gene_disruptive_sv_per_genome.by_cancer", sep="_"),
                                             paste(sex, "only", sep="_"), sep="."),
                                       main.title=paste(freq.names[freq], "gene-disruptive SVs / sample"),
                                       barplot.height=barplot.height, barplot.width=barplot.width,
                                       autosomal=TRUE)
    }
  }
}


# Number of single-gene rare/vrare/singleton LoF, CG, IEDs per genome
for(freq in c("rare", "vrare", "singleton")){
  sv.subsets <- list(
    list(paste(freq, "single_gene_LoF_SVs", sep="_"),
         as.character(unlist(sapply(data, function(info){
           rownames(info$bed)[which(sapply(info$bed$PREDICTED_LOF, length) == 1)]
         }))),
         paste(freq.names[freq], "single-gene LoF SVs")),
    list(paste(freq, "single_gene_LoF_DELs", sep="_"),
         as.character(unlist(sapply(data, function(info){
           rownames(info$bed)[which(sapply(info$bed$PREDICTED_LOF, length) == 1
                                    & info$bed$SVTYPE == "DEL")]
         }))),
         paste(freq.names[freq], "single-gene LoF dels.")),
    list(paste(freq, "single_gene_LoF_nonDELs", sep="_"),
         as.character(unlist(sapply(data, function(info){
           rownames(info$bed)[which(sapply(info$bed$PREDICTED_LOF, length) == 1
                                    & info$bed$SVTYPE != "DEL")]
         }))),
         paste(freq.names[freq], "single-gene non-del. LoF SVs")),
    list(paste(freq, "single_gene_CG_SVs", sep="_"),
         as.character(unlist(sapply(data, function(info){
           rownames(info$bed)[which(sapply(info$bed$PREDICTED_COPY_GAIN, length) == 1)]
         }))),
         paste(freq.names[freq], "single-gene CG SVs")),
    list(paste(freq, "single_gene_IED_SVs", sep="_"),
         as.character(unlist(sapply(data, function(info){
           rownames(info$bed)[which(sapply(info$bed$PREDICTED_INTRAGENIC_EXON_DUP, length) == 1)]
         }))),
         paste(freq.names[freq], "single-gene IED SVs"))
  )
  all.stats <- main.burden.wrapper(data, query=paste(freq, "single_gene_disruptive", sep="."),
                                   meta, action="count",
                                   af.fields, ac.fields, sv.subsets=sv.subsets, all.stats,
                                   paste(args$out_prefix, paste(freq, "single_gene_disruptive_sv_per_genome.by_cancer", sep="_"), sep="."),
                                   main.title=paste(freq.names[freq], "single-gene-disruptive SVs"),
                                   barplot.height=barplot.height, barplot.width=barplot.width)
  for(sex in c("MALE", "FEMALE")){
    all.stats <- main.burden.wrapper(data, query=paste(freq, "single_gene_disruptive",
                                                       paste(sex, "only", sep="_"), sep="."),
                                     meta[which(meta$inferred_sex == sex), ], action="count",
                                     af.fields, ac.fields, sv.subsets=sv.subsets, all.stats,
                                     paste(args$out_prefix, paste(freq, "single_gene_disruptive_sv_per_genome.by_cancer",
                                                                  paste(sex, "only", sep="_"), sep="_"), sep="."),
                                     main.title=paste(freq.names[freq], "single-gene-disruptive SVs"),
                                     barplot.height=barplot.height, barplot.width=barplot.width,
                                     autosomal=TRUE)
  }
}


# Number of genes impacted by rare/vrare/singleton LoF per genome
for(freq in c("rare", "vrare", "singleton")){
  ad.vals <- get.ad.values(data, query=paste(freq, "genes_disrupted.lof", sep="."),
                           action="sum", af.fields, ac.fields)
  new.stats <- burden.test(data, paste(freq, "genes_disrupted.lof", sep="."),
                           meta, ad.vals, binomial(), af.fields, ac.fields)
  all.stats <- rbind(all.stats, new.stats)
  pdf(paste(args$out_prefix, freq, "n_genes_disrupted_per_genome.by_cancer.pdf", sep="."),
      height=barplot.height, width=barplot.width)
  barplot.by.phenotype(stats2barplotdf(new.stats, ci.mode="normal"),
                       title=paste("# Genes with", tolower(freq.names[freq]),
                                   "LoF SV per sample"))
  dev.off()
}


# Carrier rates of rare/vrare/singleton LoF/CG/IED in gene lists
if(!is.null(args$gene_lists)){
  for(freq in c("rare", "vrare", "singleton")){
    for(i in 1:length(gene.lists)){

      set.name <- names(gene.lists)[i]
      set.lower <- tolower(gsub(" ", "_", set.name, fixed=T))
      gene.list <- gene.lists[[i]]

      # Get indexes for SVs with predicted effects on any gene in gene set
      lof.idx.list <- lapply(data, function(info){
        lof <- which(sapply(info$bed[, c("PREDICTED_LOF")], function(g){any(gene.list %in% g)}))
        ped <- which(sapply(info$bed[, c("PREDICTED_PARTIAL_EXON_DUP")], function(g){any(gene.list %in% g)}))
        sort(unique(c(lof, ped)))
      })
      cg.idx.list <- lapply(data, function(info){
        sort(unique(which(sapply(info$bed[, c("PREDICTED_COPY_GAIN")], function(g){any(gene.list %in% g)}))))
      })
      ied.idx.list <- lapply(data, function(info){
        sort(unique(which(sapply(info$bed[, c("PREDICTED_INTRAGENIC_EXON_DUP")], function(g){any(gene.list %in% g)}))))
      })
      all.idx.list <- lapply(1:length(data), function(k){
        sort(unique(c(lof.idx.list[[k]], cg.idx.list[[k]], ied.idx.list[[k]])))
      })
      sv.subsets <- list(
        list(paste(freq, set.lower, "LoF", sep="_"),
             unlist(lapply(1:length(data), function(k){rownames(data[[k]]$bed)[lof.idx.list[[k]]]})),
             paste(freq.names[freq], set.name, "LoF")),
        list(paste(freq, set.lower, "CG", sep="_"),
             unlist(lapply(1:length(data), function(k){rownames(data[[k]]$bed)[cg.idx.list[[k]]]})),
             paste(freq.names[freq], set.name, "CG")),
        list(paste(freq, set.lower, "IED", sep="_"),
             unlist(lapply(1:length(data), function(k){rownames(data[[k]]$bed)[ied.idx.list[[k]]]})),
             paste(freq.names[freq], set.name, "IED"))
      )
      all.stats <- main.burden.wrapper(data, query=paste(freq, set.lower, "gene_disruptive", sep="."),
                                       meta, action="any",
                                       af.fields, ac.fields, sv.subsets=sv.subsets, all.stats,
                                       keep.idx.list=all.idx.list,
                                       paste(args$out_prefix, paste(freq, "disruptive_sv_carrier_rate", sep="_"),
                                             set.lower, "by_cancer", sep="."),
                                       main.title=paste(freq.names[freq], set.name, "disruption"),
                                       barplot.height=barplot.height, barplot.width=barplot.width,
                                       barplot.units="percent")
    }
  }
}


# LoF constrained rare SVs after excluding COSMIC/CPG
if(!is.null(args$gene_lists)){
  if(length(intersect(names(gene.lists), c("LoF Constrained", "COSMIC", "CPG"))) == 3){

    # Exclude all SVs with a predicted coding consequence on any COSMIC gene or CPG
    x.genes <- sort(unique(unlist(gene.lists[c("COSMIC", "CPG")])))
    x.csqs <- c("PREDICTED_COPY_GAIN", "PREDICTED_DUP_PARTIAL", "PREDICTED_INTRAGENIC_EXON_DUP",
                "PREDICTED_LOF", "PREDICTED_PARTIAL_EXON_DUP")
    no.xgenes <- lapply(data, function(l){
      which(apply(l$bed[, x.csqs], 1, function(gl){
        length(intersect(unlist(gl), x.genes)) == 0
      }))
    })

    # Set test parameters
    freq <- "rare"
    set.name <- "LoF Constrained"
    set.lower <- tolower(gsub(" ", "_", set.name, fixed=T))
    gene.list <- gene.lists[[set.name]]

    # Get indexes for SVs with predicted effects on any constrained gene
    lof.idx.list <- lapply(data, function(info){
      lof <- which(sapply(info$bed[, c("PREDICTED_LOF")], function(g){any(gene.list %in% g)}))
      ped <- which(sapply(info$bed[, c("PREDICTED_PARTIAL_EXON_DUP")], function(g){any(gene.list %in% g)}))
      sort(unique(c(lof, ped)))
    })
    cg.idx.list <- lapply(data, function(info){
      sort(unique(which(sapply(info$bed[, c("PREDICTED_COPY_GAIN")], function(g){any(gene.list %in% g)}))))
    })
    ied.idx.list <- lapply(data, function(info){
      sort(unique(which(sapply(info$bed[, c("PREDICTED_INTRAGENIC_EXON_DUP")], function(g){any(gene.list %in% g)}))))
    })
    all.idx.list <- lapply(1:length(data), function(k){
      sort(intersect(unique(c(lof.idx.list[[k]], cg.idx.list[[k]], ied.idx.list[[k]])), no.xgenes[[k]]))
    })

    # Run burden test
    all.stats <- main.burden.wrapper(data, query=paste(freq, set.lower, "gene_disruptive", sep="."),
                                     meta, action="any", af.fields, ac.fields,
                                     sv.subsets=NULL, all.stats=all.stats,
                                     keep.idx.list=all.idx.list,
                                     out.prefix=paste(args$out_prefix, paste(freq, "disruptive_sv_carrier_rate", sep="_"),
                                                      set.lower, "by_cancer", "no_COSMIC_no_CPG", sep="."),
                                     main.title=paste(freq.names[freq], set.name, "disruption"),
                                     barplot.height=barplot.height, barplot.width=barplot.width,
                                     barplot.units="percent")
  }
}

# Write all stats to outfile
colnames(all.stats)[1] <- paste("#", colnames(all.stats)[1], sep="")
write.table(all.stats,
            paste(args$out_prefix, "global_burden_tests.tsv", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)
