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
  if(length(intersect(query.parts, c("notsmall", "large", "karyotypic"))) > 0
     & length(intersect(names(sv.colors), query.parts)) == 0
     & !("balanced" %in% query.parts)){
    if(autosomal){
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
# args <- list("bed" = c("~/scratch/PedSV.v2.5.3.case_control_cohort.analysis_samples.sites.bed.gz"),
#              "ad" = c("~/scratch/PedSV.v2.5.3.case_control_cohort.analysis_samples.allele_dosages.bed.gz"),
#              "metadata" = "~/scratch/PedSV.v2.5.3.cohort_metadata.w_control_assignments.tsv.gz",
#              "gene_lists" = NULL,
#              "genomic_disorder_hits" = NULL,
#              "subset_samples" = "~/scratch/PedSV.v2.5.3.case_control_analysis_cohort.samples.list",
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
#              "exclude_variants" = NULL,
#              "af_field" = "POPMAX_AF",
#              "ac_field" = "AC",
#              "out_prefix" = "~/scratch/PedSV.v2.5.3.trio.dev")
# args <- list("bed" = c("~/scratch/PedSV.v2.5.3.full_cohort.analysis_samples.sites.bed.gz"),
#              "ad" = c("~/scratch/PedSV.v2.5.3.full_cohort.analysis_samples.allele_dosages.bed.gz"),
#              "metadata" = "~/scratch/PedSV.v2.5.3.cohort_metadata.w_control_assignments.tsv.gz",
#              "gene_lists" = "~/scratch/PedSV.gene_lists.tsv",
#              "genomic_disorder_hits" = NULL,
#              "subset_samples" = "~/scratch/PedSV.v2.5.3.final_analysis_cohort.samples.list",
#              "exclude_variants" = NULL,
#              "af_field" = "POPMAX_AF",
#              "ac_field" = "AC",
#              "out_prefix" = "~/scratch/PedSV.v2.5.3.full_cohort.dev")

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


# Count of large autosomal variants per genome by SV type
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


# Carrier rate of rare/vrare/singleton large (1Mb) & notsmall (100kb) variants per genome by SV type
for(freq in c("rare", "vrare", "singleton")){
  for(size in c("large", "notsmall")){
    sv.subsets <- lapply(c("DEL", "DUP", "INV", "CPX"), function(svtype){
      list(svtype,
           as.character(unlist(sapply(data, function(info){
             rownames(info$bed)[which(info$bed$SVTYPE == svtype)]
           }))),
           paste("Samples with", freq.names[freq], sv.abbreviations[svtype],
                 if(svtype != "CTX"){if(size == "large"){">1Mb"}else{">100kb"}}))
    })
    all.stats <- main.burden.wrapper(data, query=paste(size, freq, sep="."), meta, action="any",
                                     af.fields, ac.fields, sv.subsets=sv.subsets, all.stats,
                                     paste(args$out_prefix, paste(freq, size, "sv_per_genome.by_cancer", sep="_"), sep="."),
                                     main.title=paste("Samples with", freq.names[freq], "SV",
                                                      if(size == "large"){">1Mb"}else{">100kb"}),
                                     barplot.height=barplot.height, barplot.width=barplot.width,
                                     barplot.units="percent", autosomal=FALSE)
    # Add one separate test for all large *unbalanced* SVs
    # (DEL + DUP + aneuploidy + qualifying CPX)
    all.stats <- main.burden.wrapper(data, query=paste(size, freq, "unbalanced", sep="."), meta, action="any",
                                     af.fields, ac.fields, sv.subsets=NULL, all.stats,
                                     paste(args$out_prefix, paste(freq, size, "unbalanced_sv_per_genome.by_cancer", sep="_"), sep="."),
                                     main.title=paste("Pct. w/", freq.names[freq],
                                                      "Unbalanced SV",
                                                      if(size == "large"){">1Mb"}else{">100kb"}),
                                     barplot.height=barplot.height, barplot.width=barplot.width,
                                     barplot.units="percent", autosomal=FALSE)
    # Add one separate test for all large *unbalanced* *autosomal* SVs
    # (DEL + DUP + aneuploidy + qualifying CPX)
    all.stats <- main.burden.wrapper(data, query=paste(size, freq, "unbalanced", sep="."), meta, action="any",
                                     af.fields, ac.fields, sv.subsets=NULL, all.stats,
                                     paste(args$out_prefix, paste(freq, size, "unbalanced_sv_per_genome.by_cancer.autosomal_only", sep="_"), sep="."),
                                     main.title=paste("Pct. w/", freq.names[freq],
                                                      "Unbalanced SV",
                                                      if(size == "large"){">1Mb"}else{">100kb"}),
                                     barplot.height=barplot.height, barplot.width=barplot.width,
                                     barplot.units="percent", autosomal=TRUE,
                                     custom.hypothesis=paste(size, freq, "unbalanced", "autosomal_only", sep="."))
    # Add one separate test for all large *unbalanced* *allosomal* SVs
    # (DEL + DUP + aneuploidy + qualifying CPX)
    all.stats <- main.burden.wrapper(data, query=paste(size, freq, "unbalanced", sep="."), meta, action="any",
                                     af.fields, ac.fields, sv.subsets=NULL, all.stats,
                                     keep.idx.list=lapply(data, function(d){which(d$bed$chrom %in% c("chrX", "chrY"))}),
                                     paste(args$out_prefix, paste(freq, size, "unbalanced_sv_per_genome.by_cancer.allosomal_only", sep="_"), sep="."),
                                     main.title=paste("Pct. w/", freq.names[freq],
                                                      "Unbalanced SV",
                                                      if(size == "large"){">1Mb"}else{">100kb"}),
                                     barplot.height=barplot.height, barplot.width=barplot.width,
                                     barplot.units="percent", autosomal=FALSE,
                                     custom.hypothesis=paste(size, freq, "unbalanced", "allosomal_only", sep="."))
    # Add one separate test for all large *balanced* SVs
    # (INS + INV + qualifying CPX)
    all.stats <- main.burden.wrapper(data, query=paste(size, freq, "balanced", sep="."), meta, action="any",
                                     af.fields, ac.fields, sv.subsets=NULL, all.stats,
                                     paste(args$out_prefix, paste(freq, size, "balanced_sv_per_genome.by_cancer", sep="_"), sep="."),
                                     main.title=paste("Pct. w/", freq.names[freq],
                                                      "Balanced SV",
                                                      if(size == "large"){">1Mb"}else{">100kb"}),
                                     barplot.height=barplot.height, barplot.width=barplot.width,
                                     barplot.units="percent", autosomal=FALSE)
  }
  # Add one separate test for all huge (>5Mb) *unbalanced* SVs
  # (DEL + DUP + aneuploidy + qualifying CPX)
  all.stats <- main.burden.wrapper(data, query=paste("karyotypic", freq, "unbalanced", sep="."), meta, action="any",
                                   af.fields, ac.fields, sv.subsets=NULL, all.stats,
                                   paste(args$out_prefix, paste(freq, "karyotypic_unbalanced_sv_per_genome.by_cancer", sep="_"), sep="."),
                                   main.title=paste("Pct. w/", freq.names[freq],
                                                    "Unbalanced SV >5Mb"),
                                   barplot.height=barplot.height, barplot.width=barplot.width,
                                   barplot.units="percent", autosomal=FALSE)
}



# Focused secondary analyses of largest rare unbalanced SV per genome
for(freq in c("rare", "vrare", "singleton")){
  # Collect all SV data
  largest.sv.data <-
    supercategory.burden.test(data=data,
                              query=paste(freq, "unbalanced.notsmall.genomic_imbalance", sep="."),
                              meta=meta, family=binomial(), action="max",
                              af.fields=af.fields, ac.fields=ac.fields)
  all.stats <- rbind(all.stats, largest.sv.data$new.stats)

  # K-M style visualization of largest unbalanced SV vs. case status
  largest.sv <- apply(largest.sv.data$ad.df[[1]], 2, max, na.rm=T)
  cancers.km.layer.ordered <- c("control", "pancan")
  surv.models <- lapply(cancers.km.layer.ordered, function(cancer){
    elig.sids <- intersect(get.eligible.samples(meta, cancer)$cases, names(largest.sv))
    survfit(Surv(log10(largest.sv[elig.sids]), rep(1, length(elig.sids))) ~ 1)
  })
  ylims <- c(0, max(sapply(surv.models, function(ss){
    max(ss$surv[which(ss$time>log10(50000))], na.rm=T)
  }), na.rm=T) + 0.025)
  pdf(paste(args$out_prefix, freq, "genomic_imbalance_km", "pdf", sep="."),
      height=2, width=2.5)
  svlen.line.plot(x.svlen=lapply(surv.models, function(sm){sm$time}),
                  y.value=lapply(surv.models, function(sm){sm$surv}),
                  ci.lower=lapply(surv.models, function(sm){sm$lower}),
                  ci.upper=lapply(surv.models, function(sm){sm$upper}),
                  colors=cancer.colors[cancers.km.layer.ordered], ci.alpha=0,
                xlab=paste("Largest", freq.names[freq], "SV"),
                ylab="Samples (%)", xlim=log10(c(50000, 5000000)),
                lwds=c(3, 3), x.axis.labels=c("50kb", "500kb", "5Mb"),
                ylims=ylims, y.axis.units="percent",
                x.axis.labels.at=log10(c(50000, 500000, 5000000)),
                parmar=c(2, 3, 0.25, 1))
  rect(xleft=log10(500000), xright=log10(5000000), ybottom=0, ytop=0.05,
       col=NA, lty=2, xpd=T)
  dev.off()

  # Small inset plot starting at 500kb
  pdf(paste(args$out_prefix, freq, "genomic_imbalance_km.inset", "pdf", sep="."),
      height=0.6*(12/5), width=0.65*(12/5))
  svlen.line.plot(x.svlen=lapply(surv.models, function(sm){sm$time}),
                y.value=lapply(surv.models, function(sm){sm$surv}),
                ci.lower=lapply(surv.models, function(sm){sm$lower}),
                ci.upper=lapply(surv.models, function(sm){sm$upper}),
                colors=cancer.colors[cancers.km.layer.ordered], ci.alpha=0,
                xlab=NA, ylab=NA, xlim=log10(c(500000, 5000000)),
                ylims=c(0, 0.05), y.axis.units="percent",
                lwds=c(rep(2, length(cancers.km.layer.ordered)-2), 3, 3),
                x.axis.labels=c("1Mb", "5Mb"),
                x.axis.labels.at=log10(c(1000000, 5000000)), x.tck=-0.025,
                parmar=c(1.15, 1.75, 0.4, 0.75))
  dev.off()

  # Burden test series every log-step
  ad.df <- get.ad.values(data, query=paste(freq, "notsmall", "unbalanced", sep="."),
                         action="verbose", af.fields=af.fields, ac.fields=ac.fields,
                         autosomal=FALSE)
  unbal.sv.size <- unlist(lapply(data, function(l){v <- calc.genomic.imbalance(l$bed); names(v) <- rownames(l$bed); return(v)}))
  size.burden.stats <- empty.stats.df
  size.or.x <- 10^seq(4, log10(5000000), 0.05)
  for(size in size.or.x){
    qual.sv.ids <- names(unbal.sv.size)[which(unbal.sv.size >= size)]
    ad.vals <- unlist(lapply(ad.df, compress.ad.matrix, action="any", keep.vids=qual.sv.ids))
    ad.vals[intersect(rownames(meta)[which(meta$any_aneuploidy)], names(ad.vals))] <- 1
    new.burden.stats <- burden.test(data, query=paste(freq, ".", size/1000, "kb", sep=""),
                                    meta=meta, ad.vals=ad.vals, family=binomial(),
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
  pdf(paste(args$out_prefix, freq, "genomic_imbalance_effect_size", "pdf", sep="."),
      height=2, width=2)
  svlen.line.plot(x.svlen=list(log10(size.or.x)),
                  y.value=list(size.burden.y.smooth),
                  ci.lower=list(size.burden.ci.lower.smooth),
                  ci.upper=list(size.burden.ci.upper.smooth),
                  step=FALSE, colors=cancer.colors["pancan"], ci.alpha=0.25,
                  xlab=paste("Largest", freq.names[freq], "SV"),
                  ylab="Odds Ratio", y.title.line=0.1, xlim=log10(c(50000, 5000000)),
                  lwds=c(3, 3), x.axis.labels=c("50kb", "500kb", "5Mb"),
                  ylims=c(min(c(0.8, min(size.burden.y.smooth))),
                          max(c(1.2, 1.1*size.burden.y.smooth))),
                  x.axis.labels.at=log10(c(50000, 500000, 5000000)),
                  parmar=c(2, 2, 0.25, 1))
  dev.off()
}


# TODO: need to implement this in a more memory-efficient way
# # Sum total of genomic imbalance per genome by SV type
# for(freq in c("rare", "vrare", "singleton")){
#   sv.subsets <- lapply(c("DEL", "DUP"), function(svtype){
#     list(svtype,
#          as.character(unlist(sapply(data, function(info){
#            rownames(info$bed)[which(info$bed$SVTYPE == svtype)]
#          }))),
#          paste(freq.names[freq], sv.abbreviations[svtype], "Alleles per Sample"))
#   })
#   all.stats <- main.burden.wrapper(data, query=paste(freq, "genomic_imbalance"), meta, action="sum",
#                                    af.fields, ac.fields, sv.subsets=sv.subsets, all.stats,
#                                    paste(args$out_prefix, paste(freq, "total_genomic_imbalance.by_cancer", sep="_"), sep="."),
#                                    main.title=paste(freq.names[freq], "Dosage Imbalance (log10)"),
#                                    barplot.height=barplot.height, barplot.width=barplot.width)
# }


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


# TODO: need to implement this in a more memory-efficient way
# # Number of rare/vrare/singleton variants (of any size) per genome by SV type
# for(freq in c("rare", "vrare", "singleton")){
#   sv.subsets <- lapply(c("DEL", "DUP", "INV", "CPX", "CTX"), function(svtype){
#     list(svtype,
#          as.character(unlist(sapply(data, function(info){
#            rownames(info$bed)[which(info$bed$SVTYPE == svtype)]
#          }))),
#          paste(freq.names[freq], sv.abbreviations[svtype], "Alleles per Sample"))
#   })
#   all.stats <- main.burden.wrapper(data, query=freq, meta, action="sum",
#                                    af.fields, ac.fields, sv.subsets=sv.subsets, all.stats,
#                                    paste(args$out_prefix, paste(freq, "sv_per_genome.by_cancer", sep="_"), sep="."),
#                                    main.title=paste(freq.names[freq], "SV Alleles per Sample"),
#                                    barplot.height=barplot.height, barplot.width=barplot.width)
# }


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
  barplot.by.phenotype(stats2barplotdf(new.stats, ci.mode="normal"),
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
