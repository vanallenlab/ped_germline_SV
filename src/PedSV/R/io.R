#!/usr/bin/env R

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Data I/O handlers


#' Load a principal components matrix
#'
#' Load a .tsv of genetic principal components (columns) for a set of samples (rows)
#'
#' @param tsv.in Path to input PC matrix .tsv
#' @param keep.n Keep top N PCs \[default: keep all PCs\]
#'
#' @returns data.frame
#'
#' @export load.pc.matrix
#' @export
load.pc.matrix <- function(tsv.in, keep.n=10e10){
  pc <- read.table(tsv.in, header=T, comment.char="", sep="\t", check.names=F,
                   stringsAsFactors=F)
  rownames(pc) <- pc[, 1]
  pc[, 2:ncol(pc)] <- apply(pc[, 2:ncol(pc)], 2, as.numeric)
  pc[, 1] <- NULL
  pc[, 1:keep.n]
}


#' Load relatedness metrics
#'
#' Load a .tsv of PC-Relate output metrics for pairs of samples
#'
#' @param tsv.in Path to input .tsv
#'
#' @returns data.frame
#'
#' @export load.kinship.metrics
#' @export
load.kinship.metrics <- function(tsv.in){
  kdf <- read.table(tsv.in, header=T, comment.char="", sep="\t", check.names=F,
                    stringsAsFactors=F)
  rownames(kdf) <- apply(kdf[, 1:2], 1,
                         function(pair){paste(sort(pair), collapse="|")})
  kdf[, 3:ncol(kdf)] <- apply(kdf[, 3:ncol(kdf)], 2, as.numeric)
  kdf[, -c(1:2)]
}


#" Load sample metadata
#'
#' Load a .tsv of sample metadata
#'
#' @param tsv.in Path to input .tsv
#' @param keep.samples Either vector of sample IDs to retain or path to
#' single-column flat text file of sample IDs to retain [default: keep all samples]
#' @param reassign.parents Assign all parents to have `control` disease labels [default: TRUE]
#'
#' @returns data.frame
#'
#' @details Row names in output dataframe are sample IDs
#'
#' @export load.sample.metadata
#' @export
load.sample.metadata <- function(tsv.in, keep.samples=NULL, reassign.parents=TRUE){
  # Load data
  df <- read.table(tsv.in, header=T, comment.char="", sep="\t", check.names=F,
                   stringsAsFactors=F)

  # Set row names as sample IDs
  rownames(df) <- df[, 1]
  df[, 1] <- NULL
  if(!is.null(keep.samples)){
    if(length(keep.samples) == 1){
      if(file.exists(keep.samples)){
        keep.samples <- unique(read.table(keep.samples, header=F)[, 1])
      }
    }
    df <- df[keep.samples, ]
  }

  # Rename columns
  for(cnames in list(c("ancestry_short_variant_inferred_or_reported", "reported_ancestry"),
                     c("melt_insert_size", "insert_size"),
                     c("melt_read_length", "read_length"),
                     c("sex_inferred_by_ploidy", "inferred_sex"),
                     c("ancestry_inferred_by_SVs", "inferred_ancestry"),
                     c("neuroblastoma_control", "NBL_control"),
                     c("ewing_control", "EWS_control"),
                     c("osteosarcoma_control", "OS_control"))) {
    if(cnames[1] %in% colnames(df)){
      colnames(df)[which(colnames(df) == cnames[1])] <- cnames[2]
    }
  }

  # Convert types as necessary
  df$proband <- c("Yes" = TRUE, "No" = FALSE)[df$proband]
  PedSV::load.constants("colors")
  for(cancer in names(cancer.colors)){
    cname <- paste(cancer, "control", sep="_")
    if(cname %in% colnames(df)){
      df[, cname] <- c("True" = TRUE, "False" = FALSE)[df[, cname]]
    }
  }

  # Reassign parents as controls unless disabled
  if(reassign.parents){
    parent.idxs <- which(!is.na(df$family_id) & !df$proband)
    df[parent.idxs, "disease"] <- "control"
  }

  # Infer whether each sample carries an aneuploidy
  auto.aneu <- apply(df[, paste("chr", 1:22, "_CopyNumber", sep="")], 1,
                     function(ploidies){any(ploidies > 2.9 | ploidies < 1.1)})
  sex.aneu <- abs(df$chrX_CopyNumber + df$chrY_CopyNumber - 2) > 0.9
  df$autosomal_aneuploidy <- auto.aneu
  df$sex_aneuploidy <- sex.aneu
  df$any_aneuploidy <- (auto.aneu | sex.aneu)

  # Reorder columns and sort on sample ID before returning
  df[sort(rownames(df)),
     c("study_phase", "batch", "study", "disease", "proband", "family_id",
       colnames(df)[grep("_control$", colnames(df))],
       "reported_ancestry", "inferred_ancestry", "inferred_sex",
       "insert_size", "median_coverage", "wgd_score",
       colnames(df)[grep("^PC", colnames(df))],
       colnames(df)[grep("_CopyNumber$", colnames(df))],
       colnames(df)[grep("_aneuploidy$", colnames(df))])]
}


#" Load SV callset BED file
#'
#' Load a .bed for an SV callset
#'
#' @param bed.in Path to input .bed
#' @param keep.coordinates Should coordinates be retained? \[Default: TRUE\]
#' @param pass.only Should only PASS variants be included? \[Default: TRUE\]
#' @param split.coding Should coding consequence annotations be split from
#' comma-delimited strings to character vectors? \[Default: TRUE\]
#' @param split.noncoding Should noncoding consequence annotations be split from
#' comma-delimited strings to character vectors? \[Default: FALSE\]
#' @param drop.vids Either a vector or a path to a file containing variant IDs
#' to be excluded drop the input file. \[Default: Keep all IDs\]
#' @param drop.cohort.frequencies A character vector of one or more prefixes
#' corresponding to cohort names for which cohort-specific frequencies should be
#' dropped. Useful for reducing in-memory size. Examples include "case_control"
#' and "trio". \[default: Keep all frequencies\]
#' @param keep.all.pop.frequencies Optional boolean indicator to retain population-
#' specific frequency info. \[default: FALSE\]
#' @param keep.all.sex.frequencies Optional boolean indicator to retain sex-specific
#' frequency info. \[default: FALSE\]
#'
#' @returns data.frame
#'
#' @details Row names in output dataframe are variant IDs
#'
#' @export load.sv.bed
#' @export
load.sv.bed <- function(bed.in, keep.coords=TRUE, pass.only=TRUE,
                        split.coding=TRUE, split.noncoding=FALSE,
                        drop.vids=NULL, drop.cohort.frequencies=c(),
                        keep.all.pop.frequencies=FALSE,
                        keep.all.sex.frequencies=FALSE){
  # Load data
  df <- read.table(bed.in, header=T, comment.char="", sep="\t", check.names=F,
                   stringsAsFactors=F)
  colnames(df)[1] <- gsub("#", "", colnames(df)[1])

  # Drop variant IDs, if specified
  if(!is.null(drop.vids)){
    if(length(drop.vids) == 1){
      if(file.exists(drop.vids)){
        drop.vids <- unique(read.table(drop.vids, header=F)[, 1])
      }
    }
    df <- df[which(!df$name %in% drop.vids), ]
  }

  # Restrict to PASS-only, if optioned
  if(pass.only){
    df <- df[which(df$FILTER %in% c("PASS", "MULTIALLELIC")), ]
  }

  # Set row names as variant IDs
  rownames(df) <- df$name
  df$name <- NULL

  # Drop unnecessary columns
  drop.cols <- c("svtype", "STRANDS", "MINSL", "NCN",
                 "GENOTYPE_CONCORDANCE", "NON_REF_GENOTYPE_CONCORDANCE",
                 "PCRMINUS_NCR",
                 colnames(df)[grep("^TRUTH_", colnames(df))],
                 colnames(df)[grep("^VAR_", colnames(df))])
  if(!keep.coords){
    drop.cols <- c(drop.cols, "chrom", "start", "end", "CHR2", "END", "END2",
                   "CPX_INTERVALS", "SOURCE")
  }
  if(!keep.all.pop.frequencies){
    pop.prefixes <- c("AFR", "AMI", "AMR", "ASJ", "EAS", "EUR", "FIN", "MID", "NFE", "SAS", "OTH")
    for(pp in pop.prefixes){
      drop.cols <- c(drop.cols, colnames(df)[grep(paste(pp, "_", sep=""), colnames(df))])
    }
  }
  if(!keep.all.sex.frequencies){
    drop.cols <- c(drop.cols, colnames(df)[grep("MALE_", colnames(df))])
  }
  for(prefix in drop.cohort.frequencies){
    prefix.cols <- grep(paste("^", prefix, "_", sep=""), colnames(df))
    if(length(prefix.cols) > 0){
      drop.cols <- c(drop.cols, colnames(df)[prefix.cols])
    }else{
      cat(paste("Warning: could not find any columns with the prefix '",
                prefix, "_' in ", bed.in, sep=""))
    }
  }
  df[, intersect(unique(drop.cols), colnames(df))] <- NULL

  # Backfill missing PREDICTED_NEAREST_TSS using PREDICTED_PROMOTER
  if(length(intersect(c("PREDICTED_NEAREST_TSS", "PREDICTED_PROMOTER"), colnames(df))) == 2){
    ntss_missing <- is.na(df$PREDICTED_NEAREST_TSS)
    if(any(ntss_missing)){
      df$PREDICTED_NEAREST_TSS[which(ntss_missing)] <- df$PREDICTED_PROMOTER[which(ntss_missing)]
    }
  }

  # Fill gnomAD annotation info for variants with no match in gnomAD
  gnomad.idxs <- intersect(grep("^gnom", colnames(df)), grep("_AF$", colnames(df)))
  if(length(gnomad.idxs) > 0){
    df[, gnomad.idxs] <- apply(df[, gnomad.idxs], 2, function(vals){
      vals[which(is.na(vals))] <- 0
      return(vals)
    })
  }


  # Parse list-style columns, if optioned
  list.cols <- setdiff(colnames(df)[grep("^PREDICTED_", colnames(df))], "PREDICTED_INTERGENIC")
  if(!split.coding){
    list.cols <- list.cols[grep("NONCODING", list.cols)]
  }
  if(!split.noncoding){
    list.cols <- list.cols[grep("NONCODING", list.cols, invert=TRUE)]
  }
  df[, list.cols] <- apply(df[, list.cols], 2, function(col.vals){
    sapply(as.character(col.vals), function(str){
      vals <- strsplit(str, split=",", fixed=T)
      unlist(vals[which(!is.na(vals))])
    })
  })

  # Ensure numeric frequency columns
  freq.suffixes <- c("AC", "AN", "AF", "N_BI_GENOS", "N_HOMREF", "N_HET", "N_HOMALT",
                     "FREQ_HOMREF", "FREQ_HET", "FREQ_HOMALT", "CN_NUMBER",
                     "CN_NONDIPLOID_COUNT", "CN_NONDIPLOID_FREQ")
  for(suf in freq.suffixes){
    hits <- c(which(colnames(df) == suf),
              grep(paste("_", suf, "$", sep=""), colnames(df)))
    if(length(hits) > 0){
      df[, hits] <- apply(df[, hits], 2, as.numeric)
    }
  }

  # Special parsing for multiallelic frequency columns
  mcnv.freq.suffixes <- c("CN_COUNT", "CN_FREQ")
  for(suf in mcnv.freq.suffixes){
    hits <- c(which(colnames(df) == suf),
              grep(paste("_", suf, "$", sep=""), colnames(df)))
    if(length(hits) > 0){
      df[, hits] <- apply(df[, hits], 2, function(col.vals){
        sapply(as.character(col.vals), function(x){
          as.numeric(unlist(strsplit(as.character(x), split=",")))
        })
      })
    }
  }

  # Return
  return(df)
}


#' Load gene list(s)
#'
#' Load one or more list(s) of gene symbols
#'
#' @param input.tsv Path to two-column .tsv of gene list prefix(es) and path to gene list(s)
#'
#' @returns list of gene symbol character vectors
#'
#' @export load.gene.lists
#' @export
load.gene.lists <- function(input.tsv){
  info.df <- read.table(input.tsv, header=F, sep="\t")
  set.names <- info.df[, 1]
  sets <- lapply(info.df[, 2], function(path){sort(unique(read.table(path, header=F)[, 1]))})
  names(sets) <- set.names
  return(sets)
}
