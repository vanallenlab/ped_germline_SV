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
filter.bed <- function(bed, query, af.field="AF", ac.field="AC",
                       autosomal=TRUE, biallelic=TRUE, keep.idx=NULL){
  query.parts <- unlist(strsplit(query, split=".", fixed=T))
  if(is.null(keep.idx)){
    keep.idx <- 1:nrow(bed)
  }
  if(autosomal){
    keep.idx <- intersect(keep.idx, which(bed$chrom %in% c(1:22, paste("chr", 1:22, sep=""))))
  }
  if(biallelic){
    keep.idx <- intersect(keep.idx, which(bed$FILTER == "PASS"))
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
  lof.count <- sapply(bed$PREDICTED_LOF, length) + sapply(bed$PREDICTED_PARTIAL_EXON_DUP, length)
  cg.count <- sapply(bed$PREDICTED_COPY_GAIN, length)
  ied.count <- sapply(bed$PREDICTED_INTRAGENIC_EXON_DUP, length)
  if(length(intersect(c("gene_disruptive", "genes_disrupted"), query.parts)) > 0){
    keep.idx <- intersect(keep.idx, unique(c(which(lof.count > 0),
                                             which(cg.count > 0),
                                             which(ied.count > 0))))
  }
  if("single_gene_disruptive" %in% query.parts){
    keep.idx <- intersect(keep.idx, unique(c(which(lof.count == 1),
                                             which(cg.count == 1),
                                             which(ied.count == 1))))
  }
  if("lof" %in% query.parts){
    keep.idx <- intersect(keep.idx, which(lof.count > 0))
  }
  if("cg" %in% query.parts){
    keep.idx <- intersect(keep.idx, which(cg.count > 0))
  }
  if("ied" %in% query.parts){
    keep.idx <- intersect(keep.idx, which(ied.count > 0))
  }
  if("cg_and_ied" %in% query.parts){
    keep.idx <- intersect(keep.idx, which(cg.count > 0 | ied.count > 0))
  }
  bed[keep.idx, ]
}

# Wrapper to apply filter.bed over script input data format (list)
filter.data <- function(data, query, af.fields, ac.fields, keep.idx=NULL){
  lapply(1:length(data), function(i){
    info <- data[[i]]
    list("bed" = filter.bed(info$bed, query, af.fields[i], ac.fields[i],
                            keep.idx=keep.idx),
         "ad" = info$ad)})
}

# Wrapper to format AD values as a single vector for a given query
# Note: will return list of data frames, one per cohort, if action=="verbose"
get.ad.values <- function(data, query, action, af.fields, ac.fields, keep.idx.list=NULL){
  res <- lapply(1:length(data), function(i){
    info <- data[[i]]
    query.bed <- filter.bed(info$bed, query, af.fields[i], ac.fields[i],
                            keep.idx=keep.idx.list[[i]])
    if(query %in% c("rare.genes_disrupted")){
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
      pedsv.glm(meta, X, Y, family=family, extra.terms=extra.terms))
  })

  # Format and return data.frame of test statistics
  res.df <- as.data.frame(do.call("rbind", res), row.names=NULL)
  colnames(res.df) <- c("disease", "n.case", "case.mean", "case.stdev",
                        "n.control", "control.mean", "control.stdev",
                        "coefficient", "std.err",
                        "test.statistic", "P.value")
  res.df$hypothesis <- query
  res.df$model <- if(family$family=="binomial"){"logit"}else{"linear"}
  res.df[, c("hypothesis", "disease", "model", "n.case", "case.mean", "case.stdev",
             "n.control", "control.mean", "control.stdev",
             "coefficient", "std.err", "test.statistic", "P.value")]
}

# Supercategory burden test for a category to prepare for rapid re-testing of subsets
supercategory.burden.test <- function(data, query, meta, action, family,
                                      af.fields, ac.fields, keep.idx.list=NULL,
                                      extra.terms=NULL){
  ad.df <- get.ad.values(data, query=query, action="verbose", af.fields, ac.fields,
                         keep.idx.list=keep.idx.list)
  ad.vals <- unlist(lapply(ad.df, compress.ad.matrix, action=action))
  new.stats <- burden.test(data, query=query, meta=meta, ad.vals=ad.vals,
                           family=family, extra.terms=extra.terms)
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
                                sv.subsets, all.stats, out.prefix, keep.idx.list=NULL,
                                extra.terms=NULL, main.title="Values", barplot.height=2.5,
                                barplot.width=3, barplot.units=NULL){
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
                                                 extra.terms=extra.terms)
  new.stats <- supercategory.res[["new.stats"]]
  all.stats <- rbind(all.stats, new.stats)
  pdf(paste(out.prefix, "ALL.pdf", sep="."),
      height=barplot.height, width=barplot.width)
  barplot.by.phenotype(stats2plot(new.stats, ci.mode), title=main.title,
                       top.axis.units=barplot.units)
  dev.off()
  ad.dfs <- supercategory.res[["ad.df"]]
  for(subset.info in sv.subsets){
    ad.vals <- unlist(lapply(ad.dfs, function(df){
      compress.ad.matrix(df[intersect(rownames(df), subset.info[[2]]), ],
                         action=action)
    }))
    new.stats <- burden.test(data, query, meta, ad.vals, family, af.fields, ac.fields, extra.terms)
    new.stats$hypothesis <- paste(new.stats$hypothesis, subset.info[[1]], sep=".")
    all.stats <- rbind(all.stats, new.stats)
    pdf(paste(out.prefix, subset.info[[1]], "pdf", sep="."),
        height=barplot.height, width=barplot.width)
    barplot.by.phenotype(stats2plot(new.stats, ci.mode), title=subset.info[[3]],
                         top.axis.units=barplot.units)
    dev.off()
  }
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
parser$add_argument("--gene-lists", metavar=".tsv", type="character",
                    help="Two-column .tsv of gene lists (set name and path to gene list)")
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
# args <- list("bed" = c("~/scratch/PedSV.v1.1.validation_cohort.analysis_samples.wAFs.bed.gz"),
#              "ad" = c("~/scratch/PedSV.v1.1.validation_cohort.analysis_samples.wAFs.allele_dosages.bed.gz"),
#              "metadata" = "~/scratch/gatk_sv_pediatric_cancers_combined_cohort_metadata_3_31_23.txt",
#              "gene_lists" = "~/scratch/PedSV.gene_lists.tsv",
#              "subset_samples" = "/Users/collins/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/data/ancestry_and_relatedness/PedSV.v1.validation_cohort_final_samples.list",
#              "exclude_variants" = "~/scratch/PedSV.validation.bad_IDs.dev.txt",
#              "af_field" = "POPMAX_AF",
#              "ac_field" = "AC",
#              "out_prefix" = "~/scratch/PedSV.dev")

# Load BEDs and pair AD paths with each
beds <- lapply(args$bed, load.sv.bed, drop.vids=args$exclude_variants)
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
                        "model"=character(0), "n.case"=numeric(0),
                        "case.mean"=numeric(0), "n.control"=numeric(0),
                        "control.mean"=numeric(0), "coefficient"=numeric(0),
                        "std.err"=numeric(0), "test.statistic"=numeric(0),
                        "P.value"=numeric(0))

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
                                 barplot.units="percent")


# Carrier rate of rare, large variants per genome by SV type
sv.subsets <- lapply(c("DEL", "DUP", "INV", "CPX", "CTX"), function(svtype){
  list(svtype,
       as.character(unlist(sapply(data, function(info){
         rownames(info$bed)[which(info$bed$SVTYPE == svtype)]
       }))),
       paste("Samples with Rare", sv.abbreviations[svtype],
             if(svtype != "CTX"){">1Mb"}))
})
all.stats <- main.burden.wrapper(data, query="large.rare", meta, action="any",
                                 af.fields, ac.fields, sv.subsets=sv.subsets, all.stats,
                                 paste(args$out_prefix, "rare_large_sv_per_genome.by_cancer", sep="."),
                                 main.title="Samples with Rare SV >1Mb",
                                 barplot.height=barplot.height, barplot.width=barplot.width,
                                 barplot.units="percent")


# Number of rare LoF, CG, IEDs per genome
sv.subsets <- list(
  list("rare_LoF_SVs",
       as.character(unlist(sapply(data, function(info){
         rownames(info$bed)[which(sapply(info$bed$PREDICTED_LOF, length) > 0)]
       }))),
       "Rare LoF SVs per Sample"),
  list("rare_LoF_DELs",
       as.character(unlist(sapply(data, function(info){
         rownames(info$bed)[which(sapply(info$bed$PREDICTED_LOF, length) > 0
                                  & info$bed$SVTYPE == "DEL")]
       }))),
       "Rare LoF Dels. per Sample"),
  list("rare_LoF_nonDELs",
       as.character(unlist(sapply(data, function(info){
         rownames(info$bed)[which(sapply(info$bed$PREDICTED_LOF, length) > 0
                                  & info$bed$SVTYPE != "DEL")]
       }))),
       "Rare LoF SVs (No Dels.)"),
  list("rare_CG_SVs",
       as.character(unlist(sapply(data, function(info){
         rownames(info$bed)[which(sapply(info$bed$PREDICTED_COPY_GAIN, length) > 0)]
       }))),
       "Rare CG SVs per Sample"),
  list("rare_IED_SVs",
       as.character(unlist(sapply(data, function(info){
         rownames(info$bed)[which(sapply(info$bed$PREDICTED_INTRAGENIC_EXON_DUP, length) > 0)]
       }))),
       "Rare IED SVs per Sample"),
  list("rare_nonLoF_disruptive_SVs",
       as.character(unlist(sapply(data, function(info){
         rownames(info$bed)[which(sapply(info$bed$PREDICTED_LOF, length) == 0)]
       }))),
       "Rare Non-LoF Gene-Disruptive SVs")
)
all.stats <- main.burden.wrapper(data, query="rare.gene_disruptive", meta, action="count",
                                 af.fields, ac.fields, sv.subsets=sv.subsets, all.stats,
                                 paste(args$out_prefix, "rare_gene_disruptive_sv_per_genome.by_cancer", sep="."),
                                 main.title="Rare Gene-Disruptive SVs / Sample",
                                 barplot.height=barplot.height, barplot.width=barplot.width)


# Number of single-gene rare LoF, CG, IEDs per genome
sv.subsets <- list(
  list("rare_single_gene_LoF_SVs",
       as.character(unlist(sapply(data, function(info){
         rownames(info$bed)[which(sapply(info$bed$PREDICTED_LOF, length) == 1)]
       }))),
       "Rare Single-Gene LoF SVs"),
  list("rare_single_gene_LoF_DELs",
       as.character(unlist(sapply(data, function(info){
         rownames(info$bed)[which(sapply(info$bed$PREDICTED_LOF, length) == 1
                                  & info$bed$SVTYPE == "DEL")]
       }))),
       "Rare Single-Gene LoF Dels."),
  list("rare_single_gene_LoF_nonDELs",
       as.character(unlist(sapply(data, function(info){
         rownames(info$bed)[which(sapply(info$bed$PREDICTED_LOF, length) == 1
                                  & info$bed$SVTYPE != "DEL")]
       }))),
       "Rare Single-Gene LoF SVs (No Dels.)"),
  list("rare_single_gene_CG_SVs",
       as.character(unlist(sapply(data, function(info){
         rownames(info$bed)[which(sapply(info$bed$PREDICTED_COPY_GAIN, length) == 1)]
       }))),
       "Rare Single-Gene CG SVs"),
  list("rare_single_gene_IED_SVs",
       as.character(unlist(sapply(data, function(info){
         rownames(info$bed)[which(sapply(info$bed$PREDICTED_INTRAGENIC_EXON_DUP, length) == 1)]
       }))),
       "Rare Single-Gene IED SVs")
)
all.stats <- main.burden.wrapper(data, query="rare.single_gene_disruptive", meta, action="count",
                                 af.fields, ac.fields, sv.subsets=sv.subsets, all.stats,
                                 paste(args$out_prefix, "rare_single_gene_disruptive_sv_per_genome.by_cancer", sep="."),
                                 main.title="Rare Single-Gene-Disruptive SVs",
                                 barplot.height=barplot.height, barplot.width=barplot.width)


# Number of genes impacted by rare LoF per genome
ad.vals <- get.ad.values(data, query="rare.genes_disrupted.lof",
                         action="sum", af.fields, ac.fields)
new.stats <- burden.test(data, "rare.genes_disrupted.lof",
                         meta, ad.vals, binomial(), af.fields, ac.fields)
all.stats <- rbind(all.stats, new.stats)
pdf(paste(args$out_prefix, "rare.n_genes_disrupted_per_genome.by_cancer.pdf", sep="."),
    height=barplot.height, width=barplot.width)
barplot.by.phenotype(stats2plot(new.stats, ci.mode="normal"),
                     title="# Genes with Rare LoF SV per Sample")
dev.off()

# Carrier rates of LoF/CG/IED in gene lists
if(!is.null(args$gene_lists)){
  gene.lists <- load.gene.lists(args$gene_lists)
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
      list(paste("rare", set.lower, "LoF", sep="_"),
           unlist(lapply(1:length(data), function(k){rownames(data[[k]]$bed)[lof.idx.list[[k]]]})),
           paste("Samples w/", set.name, "LoF")),
      list(paste("rare", set.lower, "CG", sep="_"),
           unlist(lapply(1:length(data), function(k){rownames(data[[k]]$bed)[cg.idx.list[[k]]]})),
           paste("Samples w/", set.name, "CG")),
      list(paste("rare", set.lower, "IED", sep="_"),
           unlist(lapply(1:length(data), function(k){rownames(data[[k]]$bed)[ied.idx.list[[k]]]})),
           paste("Samples w/", set.name, "IED"))
    )
    all.stats <- main.burden.wrapper(data, query=paste("rare", set.lower, "gene_disruptive", sep="."),
                                     meta, action="any",
                                     af.fields, ac.fields, sv.subsets=sv.subsets, all.stats,
                                     keep.idx.list=all.idx.list,
                                     paste(args$out_prefix, "rare_disruptive_sv_carrier_rate",
                                           set.lower, "by_cancer", sep="."),
                                     main.title=paste("Samples w/", set.name, "Disruption"),
                                     barplot.height=barplot.height, barplot.width=barplot.width,
                                     barplot.units="percent")
  }
}

# Carrier rates of LoF/CG/IED in gene lists after adjusting for baseline # of LoF & CG SVs
lof.counts <- get.ad.values(data, query="rare.gene_disrputive.lof",
                            action="count", af.fields, ac.fields)
cg.counts <- get.ad.values(data, query="rare.gene_disrputive.cg",
                            action="count", af.fields, ac.fields)
meta.ext <- meta
meta.ext$n.lof <- lof.counts[rownames(meta.ext)]
meta.ext$n.cg <- cg.counts[rownames(meta.ext)]
if(!is.null(args$gene_lists)){
  gene.lists <- load.gene.lists(args$gene_lists)
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
      list(paste("rare", set.lower, "LoF", sep="_"),
           unlist(lapply(1:length(data), function(k){rownames(data[[k]]$bed)[lof.idx.list[[k]]]})),
           paste("Samples w/", set.name, "LoF")),
      list(paste("rare", set.lower, "CG", sep="_"),
           unlist(lapply(1:length(data), function(k){rownames(data[[k]]$bed)[cg.idx.list[[k]]]})),
           paste("Samples w/", set.name, "CG")),
      list(paste("rare", set.lower, "IED", sep="_"),
           unlist(lapply(1:length(data), function(k){rownames(data[[k]]$bed)[ied.idx.list[[k]]]})),
           paste("Samples w/", set.name, "IED"))
    )
    all.stats <- main.burden.wrapper(data, query=paste("rare", set.lower, "gene_disruptive", sep="."),
                                     meta.ext, action="any",
                                     af.fields, ac.fields, sv.subsets=sv.subsets, all.stats,
                                     keep.idx.list=all.idx.list,
                                     paste(args$out_prefix, "rare_disruptive_sv_carrier_rate",
                                           set.lower, "baseline_adjusted.by_cancer", sep="."),
                                     extra.terms=c("n.lof", "n.cg"),
                                     main.title=paste("Samples w/", set.name, "Disruption"),
                                     barplot.height=barplot.height, barplot.width=barplot.width,
                                     barplot.units="percent")
  }
}

# Write all stats to outfile
colnames(all.stats)[1] <- paste("#", colnames(all.stats)[1], sep="")
write.table(all.stats,
            paste(args$out_prefix, "global_burden_tests.tsv", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)
