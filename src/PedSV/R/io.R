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
  pc <- read.table(tsv.in, header=T, comment.char="", sep="\t", check.names=F)
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
  kdf <- read.table(tsv.in, header=T, comment.char="", sep="\t", check.names=F)
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
#' @param keep.samples Vector of sample IDs to retain [default: keep all samples]
#'
#' @returns data.frame
#'
#' @details Row names in output dataframe are sample IDs
#'
#' @export load.sample.metadata
#' @export
load.sample.metadata <- function(tsv.in, keep.samples=NULL){
  # Load data
  df <- read.table(tsv.in, header=T, comment.char="", sep="\t", check.names=F)

  # Set row names as sample IDs
  rownames(df) <- df[, 1]
  df[, 1] <- NULL
  if(!is.null(keep.samples)){
    df <- df[keep.samples, ]
  }

  # Rename columns
  for( cnames in list(c("ancestry_short_variant_inferred_or_reported", "reported_ancestry"),
                      c("melt_insert_size", "insert_size"),
                      c("melt_read_length", "read_length"),
                      c("sex_inferred_by_ploidy", "inferred_sex"),
                      c("ancestry_inferred_by_SVs", "inferred_ancestry"))) {
    colnames(df)[which(colnames(df) == cnames[1])] <- cnames[2]
  }

  # Convert types as necessary
  df$proband <- c("Yes" = TRUE, "No" = FALSE)[df$proband]

  # Reorder columns and sort on sample ID before returning
  df[sort(rownames(df)),
     c("study_phase", "batch", "study", "disease", "proband", "family_id",
       "reported_ancestry", "inferred_ancestry", "inferred_sex", "read_length",
       "insert_size", "median_coverage", "melt_coverage", "wgd_score",
       colnames(df)[grep("^PC", colnames(df))],
       colnames(df)[grep("_CopyNumber$", colnames(df))])]
}


#" Load SV callset BED file
#'
#' Load a .bed for an SV callset
#'
#' @param bed.in Path to input .bed
#' @param keep.coordinates Should coordinates be retained? [Default: TRUE]
#' @param pass.only Should only PASS variants be included? [Default: TRUE]
#'
#' @returns data.frame
#'
#' @details Row names in output dataframe are variant IDs
#'
#' @export load.sv.bed
#' @export
load.sv.bed <- function(bed.in, keep.coords=TRUE, pass.only=TRUE){
  # Load data
  df <- read.table(bed.in, header=T, comment.char="", sep="\t", check.names=F)
  colnames(df)[1] <- gsub("#", "", colnames(df)[1])

  # Restrict to PASS-only, if optioned
  if(pass.only){
    df <- df[which(df$FILTER %in% c("PASS", "MULTIALLELIC")), ]
  }

  # Set row names as variant IDs
  rownames(df) <- df$name
  df$name <- NULL

  # Drop coordinates, if optioned
  if(!keep.coordinates){
    df <- df[, -c(1:3)]
  }

  # Drop unnecessary columns
  drop.cols <- c("svtype", "STRANDS")
  if(!keep.coordinates){
    drop.cols <- c(drop.cols, "chrom", "start", "end", "CHR2", "END", "END2",
                   "CPX_INTERVALS", "SOURCE")
  }
  df[, drop.cols] <- NULL

  # Parse list-style columns

  # Ensure numeric frequency columns
  freq.suffixes <- c("AC", "AN", "AF", "N_BI_GENOS", "N_HOMREF", "N_HET", "N_HOMALT",
                     "FREQ_HOMREF", "FREQ_HET", "FREQ_HOMALT", "CN_NUMBER", "CN_COUNT",
                     "CN_FREQ", "CN_NONDIPLOID_COUNT", "CN_NONDIPLOID_FREQ")
  for(suf in freq.suffixes){
    hits <- c(which(colnames(df) == suf),
              grep(paste("_", suf, "$", sep=""), colnames(df)))
    if(length(hits) > 0){
      df[, hits] <- apply(df[, hits], 2, as.numeric)
    }
  }

}


#' Query allele dosage matrix
#'
#' Subsets and optionally compresses an allele dosage matrix
#'
#' @param ad Either a path to an allele dosage .bed.gz matrix or a data.frame
#' from a previous call of [query.ad.matrix]
#' @param query.regions List of tripartite coordinate vectors to query; see `Details`
#' @param query.ids Vector of variant IDs to query
#' @param action Action to apply to query rows; see `Details`
#'
#' @details The value for `query.regions` is expected to be a list of vectors,
#' where each vector must have three elements of the format (chromosome, start, end).
#'
#' Recognized values for `action` include:
#' * `"verbose"` : return the full query matrix \[default\]
#' * `"any"` : return numeric indicator if any of the query rows are non-zero
#' for each sample
#' * `"all"` : return numeric indicator if all of the query rows are non-zero
#' for each sample
#' * `"sum"` : return the sum of allele dosages for all query rows per sample
#'
#' Note that
#'
#' @return numeric vector or data.frame, depending on `action`
#'
#' @importFrom bedr tabix
#'
#' @export query.ad.matrix
#' @export
query.ad.matrix <- function(ad, query.regions=NULL, query.ids=NULL,
                            action="verbose"){

  # Load region(s) of interest or whole matrix if query.regions is NULL
  if(inherits(ad, "character")){
    if(is.null(query.regions)){
      ad <- read.table(ad, sep="\t", comment.char="", check.names=F)
    }else{
      require(bedr, quietly=TRUE)
      do.call("rbind", lapply(query.regions, function(coords){
        bedr::tabix(paste(coords[1], ":", coords[2], "-", coords[3], sep=""),
                    ad, verbose=FALSE)
        }))
    }
    rownames(ad) <- ad[, 4]
    ad <- ad[, -c(1:4)]
  }

  # Subset to variants of interest and return directly if optioned
  sub.df <- ad[which(rownames(ad) %in% vids), ]
  if(action == "verbose"){
    return(sub.df)
  }

  # Otherwise, apply various compression strategies prior to returning
  col.all.na <- apply(sub.df, 2, function(vals){all(is.na(vals))})
  if(action == "any"){
    query.res <- as.numeric(apply(sub.df, 2, function(vals){any(as.logical(vals), na.rm=T)}))
  }else if(action == "all"){
    query.res <- as.numeric(apply(sub.df, 2, function(vals){all(as.logical(vals), na.rm=T)}))
  }else if(action == "sum"){
    query.res <- apply(sub.df, 2, sum, na.rm=T)
  }
  query.res[col.all.na] <- NA
  names(query.res) <- colnames(sub.df)
  return(query.res)
}
