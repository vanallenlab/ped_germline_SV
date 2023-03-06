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
#' @param tsv_in Path to input PC matrix .tsv
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

