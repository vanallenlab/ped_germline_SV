#!/usr/bin/env R

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>


# Miscellaneous helper functions


#' Get trios from metadata
#'
#' Infer complete trios from sample metadata
#'
#' @param meta Sample metadata loaded by [load.sample.metadata]
#'
#' @details The function will attempt to label mothers and fathers based on chrX
#' ploidy. This may not always match expectations for trios with data derived from
#' cell lines (e.g., 1000 Genomes LCLs) given the potential for sex chromosome
#' aneuploidies arising in culture.
#'
#' @returns data.frame
#'
#' @export get.complete.trios
#' @export
get.complete.trios <- function(meta){
  # Get list of family IDs with exactly three members with one proband and two non-probands
  fam.ids <- names(which(table(meta$family_id) == 3))
  keep.fams <- sapply(fam.ids, function(fid){
    (length(which(meta$family_id == fid & meta$proband) == 1)
    & length(which(meta$family_id == fid & !meta$proband) == 2))
  })

  # Process each family ID separately
  trios <- data.frame(do.call("rbind", lapply(fam.ids[keep.fams], function(fid){
    pro.id <- rownames(meta)[which(meta$proband & meta$family_id == fid)]
    parent.ids <- rownames(meta)[which(!meta$proband & meta$family_id == fid)]

    # Assign mom/dad based on chrX ploidy rather than inferred sex labels
    # This is necessary due to cell culture artifacts (see details above)
    parent.chrX <- meta[parent.ids, "chrX_CopyNumber"]
    mom.id <- parent.ids[head(which(parent.chrX == max(parent.chrX)), 1)]
    dad.id <- setdiff(parent.ids, mom.id)
    c(fid, pro.id, mom.id, dad.id)
  })))
  colnames(trios) <- c("family", "proband", "mother", "father")
  rownames(trios) <- trios$family
  trios$family <- NULL

  return(trios)
}


#' Stretch Vector
#'
#' Expand the length of a vector by repeating (or "stuttering") values
#'
#' @param values Vector of values to be stretched
#' @param k Number of times to duplicate each element in `values`
#'
#' @examples
#' stretch.vector(values=c(1, 2, 3), k=4)
#'
#' @export stretch.vector
#' @export
stretch.vector <- function(values, k){
  as.vector(sapply(values, rep, times=k))
}

