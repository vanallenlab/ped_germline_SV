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
require(PedSV, quietly=TRUE)
require(Hmisc, quietly=TRUE)
require(argparse, quietly=TRUE)
PedSV::load.constants("all")


##################
# Data functions #
##################
# Wrapper to format AD values as a single vector for a given query
# Note: will return list of data frames, one per cohort, if action=="verbose"
get.ad.values <- function(data, query, action, af.fields, ac.fields,
                          keep.idx.list=NULL, autosomal=TRUE){
  res <- lapply(1:length(data), function(i){
    info <- data[[i]]
    query.bed <- filter.bed(info$bed, query, af.fields[i], ac.fields[i],
                            keep.idx=keep.idx.list[[i]], autosomal=autosomal)
    if("genes_disrupted" %in% unlist(strsplit(query, split=".", fixed=T))){
      weights <- apply(query.bed[, c("PREDICTED_LOF", "PREDICTED_COPY_GAIN",
                                     "PREDICTED_INTRAGENIC_EXON_DUP")],
                       1, function(vals){sum(sapply(vals, length), na.rm=T)})
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
    X <- X[overlapping.ids]
    Y <- Y[overlapping.ids]
    n.carrier.ctrl <- length(which(Y[which(X > 0 & !is.na(X))] == 0))
    n.carrier.case <- length(which(Y[which(X > 0 & !is.na(X))] == 1))
    c(cancer, length(case.ids), n.carrier.case, mean(X[case.ids], na.rm=T),
      sd(X[case.ids], na.rm=T), length(control.ids), n.carrier.ctrl,
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
                                      extra.terms=NULL, autosomal=TRUE){
  ad.df <- get.ad.values(data, query=query, action="verbose", af.fields, ac.fields,
                         keep.idx.list=keep.idx.list, autosomal=autosomal)
  ad.df <- lapply(ad.df, function(df){df[, intersect(colnames(df), rownames(meta))]})
  ad.vals <- unlist(lapply(ad.df, compress.ad.matrix, action=action))
  query.parts <- unlist(strsplit(query, split=".", fixed=T))
  if(length(intersect(query.parts, c("large", "karyotypic"))) > 0
     & length(intersect(names(sv.colors), query.parts)) == 0){
    has.aneu <- intersect(rownames(meta)[which(meta$any_aneuploidy)], names(ad.vals))
    if(length(has.aneu) > 0){
      for(sid in has.aneu){
        if(action == "any"){
          ad.vals[sid] <- 1
        }else if(action %in% c("count", "sum")){
          ad.vals[sid] <- if(is.na(ad.vals[sid])){1}else{ad.vals[sid] + 1}
        }
      }
    }
  }
  new.stats <- burden.test(data, query=query, meta=meta, ad.vals=ad.vals,
                           family=family, extra.terms=extra.terms)
  return(list("ad.df" = ad.df, "new.stats" = new.stats))
}

# Helper function to convert burden test stats to barplot-compliant data.frame
stats2plot <- function(stats, ci.mode="normal"){
  n.phenos <- nrow(stats)
  values <- as.numeric(c(stats$case.mean, stats$control.mean))
  ns <- as.numeric(c(stats$n.case, stats$n.control))
  if(ci.mode == "normal"){
    stdevs <- as.numeric(c(stats$case.stdev, stats$control.stdev))
    ci.margins <- stdevs / sqrt(ns) * qnorm(0.975)
    ci.lowers <- values - ci.margins
    ci.uppers <- values + ci.margins
  }else if(ci.mode == "binomial"){
    require(Hmisc, quietly=TRUE)
    cis <- Hmisc::binconf(x=round(values * ns), n=ns)
    ci.lowers <- cis[, "Lower"]
    ci.uppers <- cis[, "Upper"]
  }
  pvals <- as.numeric(stats$P.value)
  plot.df <- data.frame("case.value" = values[1:n.phenos],
                        "case.ci.lower" = ci.lowers[1:n.phenos],
                        "case.ci.upper" = ci.uppers[1:n.phenos],
                        "control.value" = values[(1:n.phenos) + n.phenos],
                        "control.ci.lower" = ci.lowers[(1:n.phenos) + n.phenos],
                        "control.ci.upper" = ci.uppers[(1:n.phenos) + n.phenos],
                        "p" = pvals, row.names=stats$disease)
  plot.df[intersect(names(cancer.colors[1:4]), rownames(plot.df)), ]
}

# Main wrapper function to run one category of burden test for each SV type
# Note: sv.subsets must be a list of three-element lists, where each inner list
# is of the format (prefix, SV IDs to include, plot title)
main.burden.wrapper <- function(data, query, meta, action, af.fields, ac.fields,
                                sv.subsets, all.stats, out.prefix, keep.idx.list=NULL,
                                extra.terms=NULL, main.title="Values", barplot.height=2.5,
                                barplot.width=3, barplot.units=NULL,
                                custom.hypothesis=NULL, autosomal=TRUE){
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
                                                 autosomal=autosomal)
  new.stats <- supercategory.res[["new.stats"]]
  if(!is.null(custom.hypothesis)){
    new.stats$hypothesis <- custom.hypothesis
  }
  all.stats <- rbind(all.stats, new.stats)
  pdf(paste(out.prefix, "ALL.pdf", sep="."),
      height=barplot.height, width=barplot.width)
  barplot.by.phenotype(stats2plot(new.stats, ci.mode), title=main.title,
                       top.axis.units=barplot.units)
  dev.off()
  ad.dfs <- supercategory.res[["ad.df"]]
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
      barplot.by.phenotype(stats2plot(new.stats, ci.mode), title=subset.info[[3]],
                           top.axis.units=barplot.units)
      dev.off()
    }
  }
  return(all.stats)
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
# args <- list("bed" = c("~/scratch/PedSV.v2.5.2.case_control_cohort.analysis_samples.sites.bed.gz"),
#              "ad" = c("~/scratch/PedSV.v2.5.2.case_control_cohort.analysis_samples.allele_dosages.bed.gz"),
#              "metadata" = "~/scratch/PedSV.v2.5.2.cohort_metadata.w_control_assignments.tsv.gz",
#              "gene_lists" = "~/scratch/PedSV.gene_lists.tsv",
#              "genomic_disorder_hits" = NULL,
#              "subset_samples" = "/Users/ryan/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_v2_callset_generation/v2.2/PedSV.v2.2.case_control_analysis_cohort.samples.list",
#              "exclude_variants" = NULL,
#              "af_field" = "POPMAX_AF",
#              "ac_field" = "AC",
#              "out_prefix" = "~/scratch/PedSV.v2.5.2.case_control.dev")
# args <- list("bed" = c("~/scratch/PedSV.v2.1.trio_cohort.analysis_samples.sites.bed.gz"),
#              "ad" = c("~/scratch/PedSV.v2.1.trio_cohort.analysis_samples.allele_dosages.bed.gz"),
#              "metadata" = "~/scratch/PedSV.v2.1.cohort_metadata.w_control_assignments.tsv.gz",
#              "gene_lists" = "~/scratch/PedSV.gene_lists.tsv",
#              "genomic_disorder_hits" = NULL,
#              "subset_samples" = "/Users/ryan/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_v2_callset_generation/PedSV.v2.1.trio_analysis_cohort.samples.list",
#              "exclude_variants" = NULL,
#              "af_field" = "POPMAX_AF",
#              "ac_field" = "AC",
#              "out_prefix" = "~/scratch/PedSV.v2.1.trio.dev")
# args <- list("bed" = c("~/scratch/PedSV.v2.5.2.full_cohort.analysis_samples.sites.bed.gz"),
#              "ad" = c("~/scratch/PedSV.v2.5.2.full_cohort.analysis_samples.allele_dosages.bed.gz"),
#              "metadata" = "~/scratch/PedSV.v2.5.2.cohort_metadata.w_control_assignments.tsv.gz",
#              "gene_lists" = "~/scratch/PedSV.gene_lists.tsv",
#              "genomic_disorder_hits" = NULL,
#              "subset_samples" = "~/scratch/PedSV.v2.5.2.final_analysis_cohort.samples.list",
#              "exclude_variants" = NULL,
#              "af_field" = "POPMAX_AF",
#              "ac_field" = "AC",
#              "out_prefix" = "~/scratch/PedSV.v2.5.2.full_cohort.dev")

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
meta <- load.sample.metadata(args$metadata, keep.samples=keepers, reassign.parents=FALSE)

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

# Prepare data.frame for collecting test statistics
all.stats <- data.frame("hypothesis"=character(0), "disease"=character(0),
                        "model"=character(0), "n.case"=numeric(0),
                        "case.mean"=numeric(0), "n.control"=numeric(0),
                        "control.mean"=numeric(0), "coefficient"=numeric(0),
                        "std.err"=numeric(0), "test.statistic"=numeric(0),
                        "P.value"=numeric(0), "test"=character(0))

# Set plot dimensions
barplot.height <- 0.5 + (length(unique(meta$disease)) / 4)
barplot.width <- 3

# Absolute sum of nucleotides rearranged per genome by SV type
# TODO: implement this (need to weight AD matrix by query ID)

# Count of large variants per genome by SV type
sv.subsets <- lapply(c("DEL", "DUP", "INV", "CPX"), function(svtype){
  list(svtype,
       as.character(unlist(sapply(data, function(info){
         rownames(info$bed)[which(info$bed$SVTYPE == svtype)]
       }))),
       paste("Samples with", sv.abbreviations[svtype], ">1Mb"))
})
sv.subsets <- c(sv.subsets,
                list(list("CTX",
                          as.character(unlist(sapply(data, function(info){
                            rownames(info$bed)[which(info$bed$SVTYPE == "CTX")]
                          }))),
                          "Samples with Reciprocal Translocations")))
all.stats <- main.burden.wrapper(data, query="large", meta, action="any",
                                 af.fields, ac.fields, sv.subsets=sv.subsets, all.stats,
                                 paste(args$out_prefix, "large_sv_per_genome.by_cancer", sep="."),
                                 main.title="Samples with any SV >1Mb",
                                 barplot.height=barplot.height, barplot.width=barplot.width,
                                 barplot.units="percent", autosomal=FALSE)


# Carrier rate of rare/vrare/singleton large variants per genome by SV type
for(freq in c("rare", "vrare", "singleton")){
  sv.subsets <- lapply(c("DEL", "DUP", "INV", "CPX", "CTX"), function(svtype){
    list(svtype,
         as.character(unlist(sapply(data, function(info){
           rownames(info$bed)[which(info$bed$SVTYPE == svtype)]
         }))),
         paste("Samples with", freq.names[freq], sv.abbreviations[svtype],
               if(svtype != "CTX"){">1Mb"}))
  })
  all.stats <- main.burden.wrapper(data, query=paste("large", freq, sep="."), meta, action="any",
                                   af.fields, ac.fields, sv.subsets=sv.subsets, all.stats,
                                   paste(args$out_prefix, paste(freq, "large_sv_per_genome.by_cancer", sep="_"), sep="."),
                                   main.title=paste("Samples with", freq.names[freq],
                                                    "SV >1Mb"),
                                   barplot.height=barplot.height, barplot.width=barplot.width,
                                   barplot.units="percent", autosomal=FALSE)
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
             paste("Samples w/", freq.names[freq], " GD ", sv.abbreviations[svtype]), sep="")
      })
      keep.idx.list <- lapply(data, function(d){which(rownames(d$bed) %in% gd.ids)})
      all.stats <- main.burden.wrapper(data, query=freq, meta, action="any",
                                       af.fields, ac.fields, sv.subsets=sv.subsets, all.stats,
                                       out.prefix=paste(args$out_prefix, paste(freq, "genomic_disorders.by_cancer", sep="_"), sep="."),
                                       keep.idx.list=keep.idx.list,
                                       main.title=paste("Samples with", freq.names[freq], "Genomic Disorder"),
                                       barplot.height=barplot.height, barplot.width=barplot.width,
                                       barplot.units="percent",
                                       custom.hypothesis=paste("genomic_disorders", freq, sep="."))
    }
  }
}


# Number of rare/vrare/singleton LoF, CG, IEDs per genome
for(freq in c("rare", "vrare", "singleton")){
  sv.subsets <- list(
    list(paste(freq, "LoF_SVs", sep="_"),
         as.character(unlist(sapply(data, function(info){
           rownames(info$bed)[which(sapply(info$bed$PREDICTED_LOF, length) > 0)]
         }))),
         paste(freq.names[freq], "LoF SVs per Sample")),
    list(paste(freq, "LoF_DELs", sep="_"),
         as.character(unlist(sapply(data, function(info){
           rownames(info$bed)[which(sapply(info$bed$PREDICTED_LOF, length) > 0
                                    & info$bed$SVTYPE == "DEL")]
         }))),
         paste(freq.names[freq], "LoF Dels. per Sample")),
    list(paste(freq, "LoF_nonDELs", sep="_"),
         as.character(unlist(sapply(data, function(info){
           rownames(info$bed)[which(sapply(info$bed$PREDICTED_LOF, length) > 0
                                    & info$bed$SVTYPE != "DEL")]
         }))),
         paste(freq.names[freq], "Non-Del. LoF SVs")),
    list(paste(freq, "CG_SVs", sep="_"),
         as.character(unlist(sapply(data, function(info){
           rownames(info$bed)[which(sapply(info$bed$PREDICTED_COPY_GAIN, length) > 0)]
         }))),
         paste(freq.names[freq], "CG SVs per Sample")),
    list(paste(freq, "IED_SVs", sep="_"),
         as.character(unlist(sapply(data, function(info){
           rownames(info$bed)[which(sapply(info$bed$PREDICTED_INTRAGENIC_EXON_DUP, length) > 0)]
         }))),
         paste(freq.names[freq], "IED SVs per Sample")),
    list(paste(freq, "nonLoF_disruptive_SVs", sep="_"),
         as.character(unlist(sapply(data, function(info){
           rownames(info$bed)[which(sapply(info$bed$PREDICTED_LOF, length) == 0)]
         }))),
         paste(freq.names[freq], "Non-LoF Gene-Disruptive SVs"))
  )
  all.stats <- main.burden.wrapper(data, query=paste(freq, "gene_disruptive", sep="."),
                                   meta, action="count",
                                   af.fields, ac.fields, sv.subsets=sv.subsets, all.stats,
                                   paste(args$out_prefix, paste(freq, "gene_disruptive_sv_per_genome.by_cancer", sep="_"), sep="."),
                                   main.title=paste(freq.names[freq], "Gene-Disruptive SVs / Sample"),
                                   barplot.height=barplot.height, barplot.width=barplot.width)
}


# Number of single-gene rare/vrare/singleton LoF, CG, IEDs per genome
for(freq in c("rare", "vrare", "singleton")){
  sv.subsets <- list(
    list(paste(freq, "single_gene_LoF_SVs", sep="_"),
         as.character(unlist(sapply(data, function(info){
           rownames(info$bed)[which(sapply(info$bed$PREDICTED_LOF, length) == 1)]
         }))),
         paste(freq.names[freq], "Single-Gene LoF SVs")),
    list(paste(freq, "single_gene_LoF_DELs", sep="_"),
         as.character(unlist(sapply(data, function(info){
           rownames(info$bed)[which(sapply(info$bed$PREDICTED_LOF, length) == 1
                                    & info$bed$SVTYPE == "DEL")]
         }))),
         paste(freq.names[freq], "Single-Gene LoF Dels.")),
    list(paste(freq, "single_gene_LoF_nonDELs", sep="_"),
         as.character(unlist(sapply(data, function(info){
           rownames(info$bed)[which(sapply(info$bed$PREDICTED_LOF, length) == 1
                                    & info$bed$SVTYPE != "DEL")]
         }))),
         paste(freq.names[freq], "Single-Gene Non-Del. LoF SVs")),
    list(paste(freq, "single_gene_CG_SVs", sep="_"),
         as.character(unlist(sapply(data, function(info){
           rownames(info$bed)[which(sapply(info$bed$PREDICTED_COPY_GAIN, length) == 1)]
         }))),
         paste(freq.names[freq], "Single-Gene CG SVs")),
    list(paste(freq, "single_gene_IED_SVs", sep="_"),
         as.character(unlist(sapply(data, function(info){
           rownames(info$bed)[which(sapply(info$bed$PREDICTED_INTRAGENIC_EXON_DUP, length) == 1)]
         }))),
         paste(freq.names[freq], "Single-Gene IED SVs"))
  )
  all.stats <- main.burden.wrapper(data, query=paste(freq, "single_gene_disruptive", sep="."),
                                   meta, action="count",
                                   af.fields, ac.fields, sv.subsets=sv.subsets, all.stats,
                                   paste(args$out_prefix, paste(freq, "single_gene_disruptive_sv_per_genome.by_cancer", sep="_"), sep="."),
                                   main.title=paste(freq.names[freq], "Single-Gene-Disruptive SVs"),
                                   barplot.height=barplot.height, barplot.width=barplot.width)
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
  barplot.by.phenotype(stats2plot(new.stats, ci.mode="normal"),
                       title=paste("# Genes with", freq.names[freq], "LoF SV per Sample"))
  dev.off()
}


# Carrier rates of rare/vrare/singleton LoF/CG/IED in gene lists
if(!is.null(args$gene_lists)){
  gene.lists <- load.gene.lists(args$gene_lists)
  for(freq in c("rare", "vrare")){
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
                                       main.title=paste(freq.names[freq], set.name, "Disruption"),
                                       barplot.height=barplot.height, barplot.width=barplot.width,
                                       barplot.units="percent")
    }
  }
}


# Write all stats to outfile
colnames(all.stats)[1] <- paste("#", colnames(all.stats)[1], sep="")
write.table(all.stats,
            paste(args$out_prefix, "global_burden_tests.tsv", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)

