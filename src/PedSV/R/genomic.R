#!/usr/bin/env R

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>


# Functions for manipulating genomic-style data loaded by this package


#' Filter an SV BED
#'
#' Filter an SV BED based on various categories of interest
#'
#' @param bed SV BED file loaded by [load.sv.bed]
#' @param query Terms to indicate what kind of filtering should be performed. See `details`.
#' @param af.field BED column to use for allele frequencies \[default: POPMAX_AF\]
#' @param ac.field BED column to use for allele counts \[default: AC\]
#' @param use.gnomad.freqs Should gnomAD SV AFs also be applied for
#' frequency filtering? \[default: TRUE\]
#' @param gnomad.max.freq Maximum frequency allowed in gnomAD \[default: 0.01\]
#' @param gnomad.column Specify the column header to use for gnomAD frequency
#' filtering. Only used if `use.gnomad.freqs == TRUE`. \[default: "gnomad_v3.1_sv_POPMAX_AF"\]
#' @param autosomal Keep only autosomal variants \[default: FALSE\]
#' @param pass.only Keep only variants where FILTER is PASS \[default: TRUE\]
#' @param keep.idx Optional vector to specify which row numbers from `bed` should
#' be retained. Can be useful for applying more complex or custom filtering in addition
#' to the basic functions in `query`.
#' @param return.idxs Should row indexes be returned instead of a filtered version
#' of `bed`? \[default: FALSE\]
#'
#' @details The `query` argument can be passed as either a character vector or a
#' single character string. In both cases, the argument should indicate what filter(s)
#' should be applied to `bed`. Recognized filters currently include:
#' * `DEL`|`DUP`|`CNV`|`INS`|`INV`|`CPX`|`CTX` : retain only this SV type
#' * `rare` : AF < 1%
#' * `vrare` : AF < 0.1%
#' * `singleton` : AC = 1
#' * `large` : SVLEN > 1,000,000 or SVTYPE = CTX
#' * `karyotypic` : SVLEN > 5,000,000 or SVTYPE = CTX
#' * `gene_disruptive` or `genes_disrupted` : any SV with predicted LoF, PED, CG, or IED
#' * `single_gene_disruptive` : as above, but further restricting to SVs impacting just one gene
#' * `lof` : predicted LoF and/or PED consequences
#' * `cg` : predicted CG consequence
#' * `ied` : predicted IED consequence
#' * `cg_and_ied` : predicted CG and/or IED
#'
#' If `query` is provided as a vector, it should contain any of the above terms.
#' If it is provided as a single character string, it must be period-delimited.
#'
#' If `use.gnomad.freqs` is `TRUE` and any of `rare`, `vrare`, or `singleton`
#' are included in `query`, then variants will be further filtered
#' to be < `gnomad.max.freq` according to `gnomad.column`
#'
#' @return data.frame
#'
#' @export filter.bed
#' @export
filter.bed <- function(bed, query, af.field="POPMAX_AF", ac.field="AC",
                       use.gnomad.freqs=TRUE, gnomad.max.freq=0.01,
                       gnomad.column="gnomad_v3.1_sv_POPMAX_AF",
                       autosomal=FALSE, pass.only=TRUE,
                       keep.idx=NULL, return.idxs=FALSE){
  if(length(query) == 1 & length(grep(".", query, fixed=T)) > 0){
    query.parts <- unlist(strsplit(query, split=".", fixed=T))
  }else{
    query.parts <- query
  }
  has.gnomad <- gnomad.column %in% colnames(bed)
  if(is.null(keep.idx)){
    keep.idx <- 1:nrow(bed)
  }
  if(autosomal){
    keep.idx <- intersect(keep.idx, which(bed$chrom %in% c(1:22, paste("chr", 1:22, sep=""))))
  }
  if(pass.only){
    keep.idx <- intersect(keep.idx, which(bed$FILTER == "PASS"))
  }
  for(svtype in c("DEL", "DUP", "CNV", "INS", "INV", "CPX", "CTX")){
    if(svtype %in% query.parts){
      keep.idx <- intersect(keep.idx, which(bed$SVTYPE == svtype))
    }
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
  if(length(intersect(query.parts, c("rare", "vrare", "singleton"))) > 0
     & use.gnomad.freqs & has.gnomad){
      keep.idx <- intersect(keep.idx, which(bed[, gnomad.column] < 0.01))
  }
  if("large" %in% query.parts){
    keep.idx <- intersect(keep.idx, which(bed$SVLEN > 1000000 | bed$SVTYPE == "CTX"))
  }
  if("karyotypic" %in% query.parts){
    keep.idx <- intersect(keep.idx, which(bed$SVLEN > 5000000 | bed$SVTYPE == "CTX"))
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
  if(return.idxs){
    keep.idx
  }else{
    bed[keep.idx, ]
  }
}

