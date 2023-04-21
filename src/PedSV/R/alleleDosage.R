#!/usr/bin/env R

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Functions for manipulating allele dosage matrixes


#' Query allele dosage matrix
#'
#' Subsets and optionally compresses an allele dosage matrix
#'
#' @param ad Either a path to an allele dosage .bed.gz matrix or a data.frame
#' from a previous call of [query.ad.matrix]
#' @param query.regions List of tripartite coordinate vectors to query; see `Details`
#' @param query.ids Vector of variant IDs to query
#' @param action Action to apply to query rows; see `Details`
#' @param weights Optional named vector of weights for each variant. Variants
#' not present in this vector will be assigned weight = 1. \[default: no weighting\]
#'
#' @details The value for `query.regions` is expected to be a list of vectors,
#' where each vector must either be a single chromosome or have exactly
#' three elements in the format (chromosome, start, end).
#'
#' Recognized values for `action` include:
#' * `"verbose"` : return the full query matrix \[default\]
#' * `"any"` : return numeric indicator if any of the query rows are non-zero
#' for each sample
#' * `"all"` : return numeric indicator if all of the query rows are non-zero
#' for each sample
#' * `"count"` : return the count of non-reference variants per sample
#' irrespective of genotype
#' * `"sum"` : return the sum of allele dosages for all query rows per sample
#'
#' @return numeric vector or data.frame, depending on `action`
#'
#' @importFrom bedr tabix
#'
#' @seealso [compress.ad.matrix]
#'
#' @export query.ad.matrix
#' @export
query.ad.matrix <- function(ad, query.regions=NULL, query.ids=NULL,
                            action="verbose", weights=NULL){

  # Load region(s) of interest or whole matrix if query.regions is NULL
  if(inherits(ad, "character")){
    if(is.null(query.regions)){
      ad <- read.table(ad, sep="\t", comment.char="", check.names=F,
                       stringsAsFactors=F)
    }else{
      require(bedr, quietly=TRUE)
      tabix.query <- sapply(query.regions, function(coords){
        if(length(coords) == 1){
          coords <- c(coords, 0, 400000000)
        }
        coords <- as.character(coords)
        paste(coords[1], ":", coords[2], "-", coords[3], sep="")
      })
      tabix.query <- bedr.merge.region(bedr.sort.region(tabix.query, verbose=FALSE),
                                       verbose=FALSE)
      ad <- bedr::tabix(region=tabix.query, file.name=ad)
      if(!is.null(query.ids)){
        ad <- ad[which(ad[, 4] %in% query.ids), ]
      }
    }
    ad <- as.data.frame(ad)
    ad <- ad[!duplicated(ad[, 1:4]), ]
    ad[, -c(1:4)] <- apply(ad[, -c(1:4)], 2, as.numeric)
    rownames(ad) <- ad[, 4]
    ad <- ad[, -c(1:4)]
  }

  # Second subsetting step based on query IDs
  # This is necessary for use cases when ad is passed as a data.frame
  if(!is.null(query.ids)){
    ad <- ad[which(rownames(ad) %in% query.ids), ]
  }

  # Otherwise, apply various compression strategies prior to returning
  compress.ad.matrix(ad, action, weights)
}


#' Compress allele dosage matrix
#'
#' Compresses information in an allele dosage matrix
#'
#' @param ad.df Data.frame of allele dosages loaded by [query.ad.matrix]
#' @param action Action to apply to query rows; see `Details` in [query.ad.matrix]
#' @param weights Optional named vector of weights for each variant. Variants
#' not present in this vector will be assigned weight = 1. \[default: no weighting\]
#'
#' @return numeric vector
#'
#' @seealso [query.ad.matrix]
#'
#' @export compress.ad.matrix
#' @export
compress.ad.matrix <- function(ad.df, action, weights=NULL){
  col.all.na <- apply(ad.df, 2, function(vals){
    all(is.na(as.numeric(vals)))
  })

  # Reorder weights, fill missing weights, and apply weights to all ACs
  ordered.weights <- rep(1, times=nrow(ad.df))
  if(!is.null(weights)){
    for(svid in names(weights)){
      hits <- which(rownames(ad.df) == svid)
      if(length(hits) > 0){
        ordered.weights[hits] <- weights[svid]
      }
    }
  }
  ad.df <- as.data.frame(apply(ad.df, 2, function(vals){vals * ordered.weights}))

  if(action == "verbose"){
    query.res <- ad.df
  }
  if(action == "any"){
    query.res <- as.numeric(apply(ad.df, 2, function(vals){
      any(as.logical(as.numeric(vals)), na.rm=T)
    }))
  }else if(action == "all"){
    query.res <- as.numeric(apply(ad.df, 2, function(vals){
      all(as.logical(as.numeric(vals)), na.rm=T)
    }))
  }else if(action == "count"){
    query.res <- apply(ad.df, 2, function(vals){
      sum(as.numeric(vals != 0), na.rm=T)
    })
  }else if(action == "sum"){
    query.res <- apply(ad.df, 2, function(vals){
      sum(as.numeric(vals), na.rm=T)
    })
  }
  query.res[col.all.na] <- NA
  names(query.res) <- colnames(ad.df)
  return(query.res)
}


#' Query allele dosage matrix using an SV BED file
#'
#' Helper function to format and quickly execute multiple allele dosage matrix
#' queries corresponding to the variants present in an input SV BED file
#'
#' @param ad Allele dosage matrix. See [query.ad.matrix].
#' @param bed Data frame of SVs to query as loaded by [load.sv.bed].
#' @param action Action to apply to query. See [query.ad.matrix].
#' @param weights Optional named vector of weights for each variant. Variants
#' not present in this vector will be assigned weight = 1. \[default: no weighting\]
#'
#' @return Data.frame or vector depending on value of `action`
#'
#' @seealso [query.ad.matrix], [load.sv.bed]
#'
#' @export query.ad.from.sv.bed
#' @export
query.ad.from.sv.bed <- function(ad, bed, action="verbose", weights=NULL){
  # Make query intervals list from BED
  query.regions <- lapply(1:nrow(bed), function(i){
    c(bed[i, "chrom"], bed[i, "start"] - 5, bed[i, "start"] + 5)
  })
  query.ids <- rownames(bed)

  # Query AD matrix
  query.ad.matrix(ad, query.regions, query.ids, action, weights)
}
