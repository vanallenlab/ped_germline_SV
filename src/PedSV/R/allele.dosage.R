#!/usr/bin/env R

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Functions for manipulating allele dosage matrixes


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
#' @param padding Number of bp to pad each breakpoint for query \[default: 5\]
#' @param keep.samples Optional. Either character vector of sample IDs to retain
#' or path to one-column flat text file of samples to retain. \[default: keep all samples\]
#' See [query.ad.matrix].
#'
#' @return Data.frame or vector depending on value of `action`
#'
#' @seealso [query.ad.matrix], [load.sv.bed]
#'
#' @export query.ad.from.sv.bed
#' @export
query.ad.from.sv.bed <- function(ad, bed, action="verbose", weights=NULL, padding=5,
                                 keep.samples=NULL){
  # Make query intervals list from BED & query AD matrix
  if(nrow(bed) > 0){
    query.regions <- lapply(1:nrow(bed), function(i){
      c(bed[i, "chrom"], bed[i, "start"] - padding, bed[i, "start"] + padding)
    })
    query.ids <- rownames(bed)
    query.ad.matrix(ad, query.regions, query.ids=query.ids, action=action,
                    weights=weights, keep.samples=keep.samples)
  }else{
    data.frame()
  }
}
