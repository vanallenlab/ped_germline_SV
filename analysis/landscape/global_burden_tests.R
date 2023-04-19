#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, and the Van Allen Laboratory
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
# Filter a BED based on various categories of interest
filter.bed <- function(bed, query, af.field="AF", ac.field="AC", autosomal=TRUE){
  query.parts <- unlist(strsplit(query, split=".", fixed=T))
  keep.idx <- 1:nrow(bed)
  if(autosomal){
    keep.idx <- intersect(keep.idx, which(bed$chrom %in% c(1:22, paste("chr", 1:22, sep=""))))
  }
  if("rare" %in% query.parts){
    keep.idx <- intersect(keep.idx, which(bed[, af.field] < 0.01 & bed$SVTYPE != "CNV"))
  }
  if("vrare" %in% query.parts){
    keep.idx <- intersect(keep.idx, which(bed[, af.field] < 0.001 & bed$SVTYPE != "CNV"))
  }
  if("singleton" %in% query.parts){
    keep.idx <- intersect(keep.idx, which(bed[, ac.field] <= 1 & bed$SVTYPE != "CNV"))
  }
  if("large" %in% query.parts){
    keep.idx <- intersect(keep.idx, which(bed$SVLEN > 1000000 | bed$SVTYPE == "CTX"))
  }
  lof.idx <- which(!unlist(sapply(bed$PREDICTED_LOF, is.na)))
  cg.idx <- which(!unlist(sapply(bed$PREDICTED_COPY_GAIN, is.na)))
  ied.idx <- which(!unlist(sapply(bed$PREDICTED_INTRAGENIC_EXON_DUP, is.na)))
  if("gene_disruptive" %in% query.parts){
    keep.idx <- intersect(keep.idx, unique(c(lof.idx, cg.idx, ied.idx)))
  }
  if("lof" %in% query.parts){
    keep.idx <- intersect(keep.idx, lof.idx)
  }
  if("cg" %in% query.parts){
    keep.idx <- intersect(keep.idx, cg.idx)
  }
  if("ied" %in% query.parts){
    keep.idx <- intersect(keep.idx, ied.idx)
  }
  bed[keep.idx, ]
}

# Wrapper to apply filter.bed over script input data format (list)
filter.data <- function(data, query, af.fields, ac.fields){
  lapply(1:length(data), function(i){
    info <- data[[i]]
    list("bed" = filter.bed(info$bed, query, af.fields[i], ac.fields[i]),
         "ad" = info$ad)})
}

# Wrapper to format AD values as a single vector for a given query
# Note: will return list of data frames, one per cohort, if action=="verbose"
get.ad.values <- function(data, query, action, af.fields, ac.fields){
  res <- lapply(1:length(data), function(i){
    info <- data[[i]]
    query.bed <- filter.bed(info$bed, query, af.fields[i], ac.fields[i])
    query.ad.from.sv.bed(info$ad, query.bed, action)
  })
  if(is.data.frame(res[[1]])){
    return(res)
  }else{
    unlist(res)
  }
}

# Wrapper to execute a single burden association test
burden.test <- function(data, query, meta, ad.vals, family,
                        af.fields, ac.field){
  # Run one comparison for each cancer type with at least one sample present in data
  cancers <- setdiff(unique(meta[names(ad.vals), "disease"]), "control")
  cancers <- c("pancan", metadata.cancer.label.map[cancers])
  res <- lapply(cancers, function(cancer){
    # Get list of eligible samples
    elig.sids <- get.eligible.samples(meta, cancer)
    case.ids <- elig.sids[["cases"]]
    control.ids <- elig.sids[["controls"]]

    # Run burden test
    X <- ad.vals[intersect(names(ad.vals), c(case.ids, control.ids))]
    Y <- get.phenotype.vector(case.ids, control.ids)
    overlapping.ids <- intersect(names(X), names(Y))
    X <- X[overlapping.ids]
    Y <- Y[overlapping.ids]
    c(cancer, length(case.ids), mean(X[case.ids], na.rm=T),
      sd(X[case.ids], na.rm=T), length(control.ids),
      mean(X[control.ids], na.rm=T), sd(X[control.ids], na.rm=T),
      pedsv.glm(meta, X, Y, family=family))
  })

  # Format and return data.frame of test statistics
  res.df <- as.data.frame(do.call("rbind", res), row.names=NULL)
  colnames(res.df) <- c("disease", "n.case", "case.mean", "case.stdev",
                        "n.control", "control.mean", "control.stdev",
                        "coefficient", "std.err",
                        "test.statistic", "P.value")
  res.df$hypothesis <- query
  res.df$test <- if(family$family=="binomial"){"logit"}else{"linear"}
  res.df[, c("hypothesis", "disease", "test", "n.case", "case.mean", "case.stdev",
             "n.control", "control.mean", "control.stdev",
             "coefficient", "std.err", "test.statistic", "P.value")]
}

# Supercategory burden test for a category to prepare for rapid re-testing of subsets
supercategory.burden.test <- function(data, query, meta, action, family,
                                      af.fields, ac.fields){
  ad.df <- get.ad.values(data, query=query, action="verbose", af.fields, ac.fields)
  ad.vals <- unlist(lapply(ad.df, compress.ad.matrix, action=action))
  new.stats <- burden.test(data, query=query, meta=meta, ad.vals=ad.vals,
                           family=family)
  return(list("ad.df" = ad.df, "new.stats" = new.stats))
}

# Helper function to convert burden test stats to barplot-compliant data.frame
stats2plot <- function(stats, ci.mode="normal"){
  values <- c(mean(as.numeric(stats$control.mean)), as.numeric(stats$case.mean))
  ns <- c(mean(as.numeric(stats$n.control)), as.numeric(stats$n.case))
  if(ci.mode == "normal"){
    stdevs <- c(mean(as.numeric(stats$control.stdev)), as.numeric(stats$case.stdev))
    ci.margins <- stdevs / sqrt(ns) * qnorm(0.975)
    ci.lowers <- values - ci.margins
    ci.uppers <- values + ci.margins
  }else if(ci.mode == "binomial"){
    require(Hmisc, quietly=TRUE)
    cis <- Hmisc::binconf(x=round(values * ns), n=ns)
    ci.lowers <- cis[, "Lower"]
    ci.uppers <- cis[, "Upper"]
  }
  pvals <- as.numeric(c(NA, stats$P.value))
  plot.df <- data.frame("value" = values, "ci.lower" = ci.lowers, "ci.upper" = ci.uppers,
                        "p" = pvals, row.names=c("control", stats$disease))
  plot.df[c("control", intersect(names(cancer.colors[1:4]), rownames(plot.df))), ]
}

# Main wrapper function to run one category of burden test for each SV type
# Note: sv.subsets must be a list of three-element lists, where each inner list
# is of the format (prefix, SV IDs to include, plot title)
main.burden.wrapper <- function(data, query, meta, action, af.fields, ac.fields,
                                sv.subsets, all.stats, out.prefix,
                                main.title="Values", barplot.height=2.5,
                                barplot.width=3, barplot.units=NULL){
  # Get other parameters dictated by value of action
  if(action == "any"){
    ci.mode <- "binomial"
    family <- binomial()
  }else if(action %in% c("sum", "count")){
    ci.mode <- "normal"
    family <- gaussian()
  }

  # Gather verbose AD matrix for reuse later, and test all SVs before splitting by type
  supercategory.res <- supercategory.burden.test(data, query=query, meta=meta,
                                                 action=action, family=family,
                                                 af.fields=af.fields, ac.fields=ac.fields)
  new.stats <- supercategory.res[["new.stats"]]
  all.stats <- rbind(all.stats, new.stats)
  pdf(paste(out.prefix, "ALL.pdf", sep="."),
      height=barplot.height, width=barplot.width)
  barplot.by.phenotype(stats2plot(new.stats, ci.mode), title=main.title,
                       top.axis.units=barplot.units)
  dev.off()
  ad.dfs <- supercategory.res[["ad.df"]]
  sapply(sv.subsets, function(subset.info){
    ad.vals <- unlist(lapply(ad.dfs, function(df){
      compress.ad.matrix(df[intersect(rownames(df), subset.info[[2]]), ],
                         action=action)
    }))
    new.stats <- burden.test(data, query, meta, ad.vals, family, af.fields, ac.fields)
    new.stats$hypothesis <- paste(new.stats$hypothesis, subset.info[[1]], sep=".")
    all.stats <- rbind(all.stats, new.stats)
    pdf(paste(out.prefix, subset.info[[1]], "pdf", sep="."),
        height=barplot.height, width=barplot.width)
    barplot.by.phenotype(stats2plot(new.stats, ci.mode), title=subset.info[[3]],
                         top.axis.units=barplot.units)
    dev.off()
  })
  return(all.stats)
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Plot SV counts and sizes")
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
parser$add_argument("--subset-samples", metavar=".tsv", type="character",
                    help="list of samples to subset [default: use all samples]")
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
# args <- list("bed" = c("~/scratch/PedSV.v1.1.validation_cohort.analysis_samples.wAFs.bed.gz"),
#              "ad" = c("~/scratch/PedSV.v1.1.validation_cohort.analysis_samples.wAFs.allele_dosages.bed.gz"),
#              "metadata" = "~/scratch/gatk_sv_pediatric_cancers_combined_cohort_metadata_3_31_23.txt",
#              "subset_samples" = "/Users/collins/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/data/ancestry_and_relatedness/PedSV.v1.both_cohorts_final_samples.list",
#              "af_field" = "AF",
#              "ac_field" = "AC",
#              "out_prefix" = "~/scratch/PedSV.dev")

# Load BEDs and pair AD paths with each
beds <- lapply(args$bed, load.sv.bed)
data <- lapply(1:length(beds), function(i){
  list("bed" = beds[[i]], "ad" = args$ad[i])
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
                        "test"=character(0), "n.case"=numeric(0),
                        "case.mean"=numeric(0), "n.control"=numeric(0),
                        "control.mean"=numeric(0), "coefficient"=numeric(0),
                        "std.err"=numeric(0), "test.statistic"=numeric(0),
                        "P.value"=numeric(0))

# Set plot dimensions
barplot.height <- 0.5 + (length(unique(meta$disease)) / 4)
barplot.width <- 3

# Absolute sum of nucleotides rearranged per genome by SV type
# TODO: implement this (need to weight AD matrix by query ID)

# # Count of large variants per genome by SV type
sv.subsets <- lapply(c("DEL", "DUP", "INV", "CPX"), function(svtype){
  list(svtype,
       as.character(unlist(sapply(data, function(info){
         rownames(info$bed)[which(info$bed$SVTYPE == svtype)]
       }))),
       paste(sv.abbreviations[svtype], ">1Mb per Sample"))
})
all.stats <- main.burden.wrapper(data, query="large", meta, action="count",
                                 af.fields, ac.fields, sv.subsets=sv.subsets, all.stats,
                                 paste(args$out_prefix, "large_sv_per_genome.by_cancer", sep="."),
                                 main.title="SVs >1Mb per Sample",
                                 barplot.height, barplot.width)


# Carrier rate of rare, large variants per genome by SV type
sv.subsets <- lapply(c("DEL", "DUP", "INV", "CPX", "CTX"), function(svtype){
  list(svtype,
       as.character(unlist(sapply(data, function(info){
         rownames(info$bed)[which(info$bed$SVTYPE == svtype)]
       }))),
       paste("Samples with Rare", sv.abbreviations[svtype], ">1Mb"))
})
all.stats <- main.burden.wrapper(data, query="large.rare", meta, action="any",
                                 af.fields, ac.fields, sv.subsets=sv.subsets, all.stats,
                                 paste(args$out_prefix, "rare_large_sv_per_genome.by_cancer", sep="."),
                                 main.title="Samples with Rare SV >1Mb",
                                 barplot.height, barplot.width, barplot.units="percent")

# Number of rare LoF, CG, IEDs per genome
sv.subsets <- list(
  list("rare_LoF_SVs",
       as.character(unlist(sapply(data, function(info){
         rownames(info$bed)[which(!unlist(sapply(info$bed$PREDICTED_LOF, is.na)))]
       }))),
       "Rare LoF SVs per Sample"),
  list("rare_CG_SVs",
       as.character(unlist(sapply(data, function(info){
         rownames(info$bed)[which(!unlist(sapply(info$bed$PREDICTED_COPY_GAIN, is.na)))]
       }))),
       "Rare CG SVs per Sample"),
  list("rare_IED_SVs",
       as.character(unlist(sapply(data, function(info){
         rownames(info$bed)[which(!unlist(sapply(info$bed$PREDICTED_INTRAGENIC_EXON_DUP, is.na)))]
       }))),
       "Rare IED SVs per Sample")
)
all.stats <- main.burden.wrapper(data, query="rare.gene_disruptive", meta, action="count",
                                 af.fields, ac.fields, sv.subsets=sv.subsets, all.stats,
                                 paste(args$out_prefix, "rare_gene_disruptive_sv_per_genome.by_cancer", sep="."),
                                 main.title="Rare Gene-Disruptive SVs / Sample",
                                 barplot.height, barplot.width)

# Number of genes impacted by rare LoF, CG, IED per genome
# TODO: add variant weighting option to compress.ad.matrix

# Carrier rate of rare LoF in constrained genes per genome

# Carrier rate of rare LoF in haploinsufficient genes per genome

# Carrier rate of rare CG in triplosensitive genes per genome

# Carrier rate of rare LoF in tumor suppressor genes per genome

# Carrier rate of rare CG in oncogenes per genome

# Write all stats to outfile
colnames(all.stats)[1] <- paste("#", colnames(all.stats)[1], sep="")
write.table(all.stats,
            paste(args$out_prefix, "global_burden_tests.tsv", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)
