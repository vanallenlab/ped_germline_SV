#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Conduct gene-based collapsing burden tests for rare SVs in all cancer types


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(PedSV, quietly=TRUE)
require(argparse, quietly=TRUE)
PedSV::load.constants("all")
cancers <- setdiff(names(cancer.colors), "control")
csq.sets <- list("LOF" = c("PREDICTED_LOF", "PREDICTED_PARTIAL_EXON_DUP"),
              "CG" = c("PREDICTED_COPY_GAIN"),
              "GeneDisruptive" = c("PREDICTED_LOF", "PREDICTED_PARTIAL_EXON_DUP",
                                   "PREDICTED_INTRAGENIC_EXON_DUP",
                                   "PREDICTED_COPY_GAIN"))


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Gene-based rare SV association tests")
parser$add_argument("--bed", metavar=".tsv", type="character", required=TRUE,
                    help="SV sites .bed. Required.")
parser$add_argument("--ad", metavar=".tsv", type="character", required=TRUE,
                    help="Allele dosage .bed. Required.")
parser$add_argument("--genes", metavar=".bed", type="character", required=TRUE,
                    help=paste("BED of gene coordinates and gene symbols to",
                               "be tested. Required."))
parser$add_argument("--metadata", metavar=".tsv", type="character",
                    help="sample metadata .tsv", required=TRUE)
parser$add_argument("--subset-samples", metavar=".tsv", type="character",
                    help="list of samples to subset [default: use all samples]")
parser$add_argument("--out-bed", metavar=".tsv", type="character", required=TRUE,
                    help="Path to output .tsv of summary statistics. Required.")
args <- parser$parse_args()

# # DEV:
# args <- list("bed" = c("~/scratch/PedSV.v2.1.case_control_cohort.analysis_samples.sites.bed.gz"),
#              "ad" = c("~/scratch/PedSV.v2.1.case_control_cohort.analysis_samples.allele_dosages.bed.gz"),
#              "genes" = "~/scratch/test.genes.bed",
#              "metadata" = "~/scratch/PedSV.v2.1.cohort_metadata.w_control_assignments.tsv.gz",
#              "subset_samples" = "/Users/ryan/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_v2_callset_generation/PedSV.v2.1.case_control_analysis_cohort.samples.list",
#              "out_bed" = "~/scratch/PedSV.v2.1.case_control.dev.rvas_test.bed")


# Load metadata
keepers <- NULL
if(!is.null(args$subset_samples)){
  keepers <- read.table(args$subset_samples, header=F)[, 1]
}
meta <- load.sample.metadata(args$metadata, keep.samples=keepers,
                             reassign.parents=FALSE)

# Load BED of eligible genes
genes.bed <- read.table(args$genes, header=F, sep="\t")
colnames(genes.bed) <- c("chr", "start", "end", "gene")
genes <- sort(unique(genes.bed$gene))

# Load BED and subset to rare variants impacting at least one of the eligible genes
bed <- load.sv.bed(args$bed)
bed <- filter.bed(bed, query="rare.gene_disruptive")
hit.idxs <- which(apply(bed[, grep("PREDICTED_", colnames(bed))], 1, function(l){
  any(sapply(l, function(g){any(g %in% genes)}))
}))
bed <- bed[hit.idxs, ]

# Load allele dosage matrix
ad <- query.ad.from.sv.bed(args$ad, bed, keep.samples=rownames(meta))

# Apply over genes
rvas.res <- as.data.frame(do.call("rbind", lapply(genes, function(gene){

  # Apply over consequences
  do.call("rbind", lapply(1:length(csq.sets), function(csq.idx){
    csq.name <- names(csq.sets)[csq.idx]
    csqs <- unlist(csq.sets[csq.idx])

    # Apply over diseases
    c(gene, csq.name, as.vector(sapply(cancers, function(cancer){

      # Run logit glm to test association
      esamps <- get.eligible.samples(meta, cancer)
      Y <- get.phenotype.vector(esamps$cases, esamps$controls)
      hit.idxs <- which(apply(bed[csqs], 1, function(l){
        any(unlist(l) == gene)
      }))
      if(length(hit.idxs) == 0){
        glm.res <- c(NA, NA, NA, NA, "glm")
        n.sv <- n.case <- case.nonref <- case.mean <- n.control <- control.nonref <- control.mean <- 0
      }else{
        vids <- rownames(bed)[hit.idxs]
        X <- compress.ad.matrix(ad[vids, ], action="any")
        keepers <- intersect(names(X), names(Y))
        X <- X[keepers]
        Y <- Y[keepers]
        ct <- table(data.frame(X, Y))
        if(nrow(ct) > 1 & ncol(ct) > 1){
          glm.res <- pedsv.glm(meta, X=X, Y=Y, family=binomial(),
                               extra.terms="cohort", firth.fallback=TRUE)
          glm.res[4] <- -log10(as.numeric(glm.res[4]))
        }else{
          glm.res <- c(NA, NA, NA, NA, "glm")
        }

        # Compute summary metrics
        n.sv <- length(vids)
        n.case <- length(esamps$cases)
        case.nonref <- length(which(X[esamps$cases] != 0))
        case.mean <- mean(X[esamps$cases], na.rm=T)
        n.control <- length(esamps$controls)
        control.nonref <- length(which(X[esamps$controls] != 0))
        control.mean <- mean(X[esamps$controls], na.rm=T)
      }

      c(n.sv, n.case, case.nonref, case.mean, n.control, control.nonref, control.mean, glm.res)
    })))
  }))
})))
cancer.cols <- as.character(sapply(cancers, function(cancer){
  paste(cancer, c("n_sv", "n_case", "case_nonref", "case_mean",
                  "n_control", "control_nonref", "control_mean",
                  "beta", "beta_se", "test_stat", "neglog10_p", "model"), sep=".")
}))
colnames(rvas.res) <- c("gene", "consequence", cancer.cols)

# Merge gene coordinates with summary statistics and sort by coordinate
out.res <- merge(genes.bed, rvas.res, by="gene", all=T, sort=F)
out.res <- out.res[with(out.res, order(chr, start, end, gene, consequence)),
                   c("chr", "start", "end", colnames(rvas.res))]

# Write results to output BED
write.table(out.res, args$out_bed, col.names=T, row.names=F, sep="\t", quote=F)

