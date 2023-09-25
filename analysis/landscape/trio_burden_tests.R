#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Conduct genome-wide trio-based burden tests for all cancer types in one cohort


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
# Map sample AD/AC values onto trios and compute midparent
get.trio.values <- function(trios, ad.vals){
  vals <- as.data.frame(t(apply(trios, 1, function(ids){ad.vals[ids]})))
  colnames(vals) <- colnames(trios)
  vals$midparent <- apply(vals[, 2:3], 1, mean, na.rm=T)
  vals$maxparent <- apply(vals[, 2:3], 1, max, na.rm=T)
  vals$proband.delta <- vals$proband - vals$midparent
  return(vals)
}



###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Trio-based midparent burden tests")
parser$add_argument("--bed", metavar=".tsv", type="character",
                    help="SV sites .bed", required=TRUE)
parser$add_argument("--ad", metavar=".tsv", type="character",
                    help="Allele dosage .bed", required=TRUE)
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
parser$add_argument("--af-field", metavar="string", type="character", default="parent_AF",
                    help=paste("Column header to use for AF-related analyses.",
                               "[default: parent_AF]"))
parser$add_argument("--ac-field", metavar="string", type="character", default="AC",
                    help=paste("Column header to use for AC-related analyses.",
                               "[default: AC]"))
parser$add_argument("--out-prefix", metavar="path", type="character", required=TRUE,
                    help="path/prefix for all output files")
args <- parser$parse_args()

# # DEV:
# args <- list("bed" = "~/scratch/PedSV.v1.1.trio_cohort.analysis_samples.wAFs.bed.gz",
#              "ad" = "~/scratch/PedSV.v1.1.trio_cohort.analysis_samples.wAFs.allele_dosages.bed.gz",
#              "metadata" = "~/scratch/gatk_sv_pediatric_cancers_combined_cohort_metadata_3_31_23.txt",
#              "gene_lists" = "~/scratch/PedSV.gene_lists.tsv",
#              "subset_samples" = "/Users/collins/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/data/ancestry_and_relatedness/PedSV.v1.trio_cohort_final_samples.list",
#              "af_field" = "POPMAX_parent_AF",
#              "ac_field" = "parent_AC",
#              "out_prefix" = "~/scratch/PedSV.trio.dev")

# Load BED
bed <- load.sv.bed(args$bed, drop.vids=args$exclude_variants)

# Load metadata, subset to samples of interest, and get a list of complete trios
keepers <- NULL
if(!is.null(args$subset_samples)){
  keepers <- read.table(args$subset_samples, header=F)[, 1]
}
meta <- load.sample.metadata(args$metadata, keep.samples=keepers, reassign.parents=FALSE)
trios <- get.complete.trios(meta)

# # Prepare data.frame for collecting test statistics
# all.stats <- data.frame("hypothesis"=character(0), "disease"=character(0),
#                         "model"=character(0), "n.case"=numeric(0),
#                         "case.mean"=numeric(0), "n.control"=numeric(0),
#                         "control.mean"=numeric(0), "coefficient"=numeric(0),
#                         "std.err"=numeric(0), "test.statistic"=numeric(0),
#                         "P.value"=numeric(0))
#
# # Set plot dimensions
# barplot.height <- 0.5 + (length(unique(meta$disease)) / 4)
# barplot.width <- 3

# DEV:
bed <- bed[which(bed$ALGORITHMS != "depth" & !(bed$EVIDENCE %in% c("RD", "BAF,RD"))), ]
rare.lof.bed <- filter.bed(bed, "rare.lof", af.field=args$af_field, ac.field=args$ac_field, autosomal=TRUE)
variant.weights <- sapply(bed$PREDICTED_LOF, length) + sapply(bed$PREDICTED_PARTIAL_EXON_DUP, length)
names(variant.weights) <- rownames(bed)
lof.ad <- query.ad.from.sv.bed(args$ad, rare.lof.bed, action="count")
lof.ac <- query.ad.from.sv.bed(args$ad, rare.lof.bed, action="sum")
lof.genes <- query.ad.from.sv.bed(args$ad, rare.lof.bed, action="sum", weights=variant.weights)
ad.vals <- get.trio.values(trios, lof.ad)
ac.vals <- get.trio.values(trios, lof.ac)
lof.gene.count <- get.trio.values(trios, lof.genes)

nbl.trios <- unique(meta$family_id[which(meta$disease == "neuroblastoma")])
ews.trios <- unique(meta$family_id[which(meta$disease == "ewing")])
control.trios <- unique(meta$family_id[which(meta$disease == "control")])

ac.vals.case <- ac.vals[!ac.vals$control, ]
ac.vals.control <- ac.vals[ac.vals$control, ]

t.test(ac.vals.case$proband, c(ac.vals.case$mother, ac.vals.case$father))
t.test(ac.vals.control$proband, c(ac.vals.control$mother, ac.vals.control$father))
wilcox.test(ac.vals.case$proband, c(ac.vals.case$mother, ac.vals.case$father))
wilcox.test(ac.vals.control$proband, c(ac.vals.control$mother, ac.vals.control$father))
poisson.test(x=c(sum(ac.vals.case$proband), sum(c(ac.vals.case$father, ac.vals.case$mother))), T=c(nrow(ac.vals.case), 2*nrow(ac.vals.case)))
poisson.test(x=c(sum(ac.vals.control$proband), sum(c(ac.vals.control$father, ac.vals.control$mother))), T=c(nrow(ac.vals.control), 2*nrow(ac.vals.control)))

require(rCNV2)
dens.scatter(ac.vals[!ac.vals$control, "midparent"], ac.vals[!ac.vals$control, "proband"], pt.cex=3, parmar=c(2.5, 2.5, 0.5, 0.5))
clean.axis(1, title="Midparent LoF Alleles")
clean.axis(2, title="Proband LoF Alleles")
abline(0, 1, col="gray", lty=5)
dens.scatter(ac.vals[ac.vals$control, "midparent"], ac.vals[ac.vals$control, "proband"], pt.cex=3, parmar=c(2.5, 2.5, 0.5, 0.5))
clean.axis(1, title="Midparent LoF Alleles")
clean.axis(2, title="Proband LoF Alleles")
abline(0, 1, col="gray", lty=5)


# Check constrained genes
constr <- read.table("~/scratch/PedSV_gene_lists/gnomad.v2.1.1.LoF_constrained.genes.list",
                     header=F, sep="\t")[, 1]
constr.lof <- which(sapply(bed$PREDICTED_LOF, function(g){any(g %in% constr)}))
constr.rare.lof.bed <- filter.bed(bed, "rare.lof", keep.idx=constr.lof, autosomal=TRUE)
lof.ac <- query.ad.from.sv.bed(args$ad, constr.rare.lof.bed, action="sum")
ac.vals <- get.trio.values(trios, lof.ac)
ac.vals$control <- rownames(ac.vals) %in% control.trios
boxplot(proband ~ control, data=ac.vals)
fdf <- matrix(c(length(which(ac.vals$proband == 0 & ac.vals$control)),
         length(which(ac.vals$proband == 0 & !ac.vals$control)),
         length(which(ac.vals$proband > 0 & ac.vals$control)),
         length(which(ac.vals$proband > 0 & !ac.vals$control))),
       byrow=T, nrow=2)

# # Write all stats to outfile
# colnames(all.stats)[1] <- paste("#", colnames(all.stats)[1], sep="")
# write.table(all.stats,
#             paste(args$out_prefix, "global_burden_tests.tsv", sep="."),
#             col.names=T, row.names=F, sep="\t", quote=F)
