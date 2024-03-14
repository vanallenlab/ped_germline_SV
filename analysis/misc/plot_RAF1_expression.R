#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2024-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Plot RAF1 tumor expression to highlight increased expression for duplication carriers


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(PedSV, quietly=TRUE)
require(beeswarm, quietly=TRUE)
PedSV::load.constants("all")
dup.carriers <- c("PT_5CZHK9CR", "PT_ZT2NW6WA", "SJ071781")


##################
# Data functions #
##################
# Merge TPMs for a single gene from two or more TPM matrixes
query.tpms <- function(tpm.list, gene){
  unlist(lapply(tpm.list, function(m){
    v <- m[gene, ]
    names(v) <- colnames(m)
    return(v)
  }))
}


######################
# Plotting functions #
######################


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Plot RAF1 RNA expression by tumor type")
parser$add_argument("--nbl-tpm", metavar=".tsv", type="character",
                    help="Matrix of NBL log2(TPM+1); row names are gene symbols",
                    required=TRUE)
parser$add_argument("--os-tpm", metavar=".tsv", type="character",
                    help="Matrix of OS log2(TPM+1); row names are gene symbols",
                    required=TRUE)
parser$add_argument("--metadata", metavar=".tsv", type="character",
                    help="Sample metadata .tsv", required=TRUE)
parser$add_argument("--outfile", metavar="path", type="character", required=TRUE,
                    help="Path to output .pdf")
args <- parser$parse_args()

# DEV
args <- list("nbl_tpm" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/data/matched_tumor_rnaseq/GMKF.NBL.all_samples.tpm_plus1_log2.tsv",
             "os_tpm" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/data/matched_tumor_rnaseq/StJude.OS.dx_samples.tpm_plus1_log2.tsv",
             "metadata" = "~/scratch/PedSV.v2.5.3.cohort_metadata.w_control_assignments.tsv.gz",
             "outfile" = "~/scratch/RAF1.expression.pdf")


# Load TPMs
nbl.tpm <- read.table(args$nbl_tpm, sep="\t", header=T)
os.tpm <- read.table(args$os_tpm, sep="\t", header=T)

# Load sample metadata
meta <- load.sample.metadata(args$metadata)

# Subset TPMs and metadata to overlapping samples
nbl.keepers <- sort(intersect(names(nbl.tpm), rownames(meta)))
nbl.tpm <- nbl.tpm[, nbl.keepers]
os.keepers <- sort(intersect(names(os.tpm), rownames(meta)))
os.tpm <- os.tpm[, os.keepers]
meta <- meta[c(nbl.keepers, os.keepers), ]
dup.carriers <- intersect(dup.carriers, rownames(meta))

# Get plot data
raf1 <- query.tpms(list(nbl.tpm, os.tpm), "RAF1")
tmem40 <- query.tpms(list(nbl.tpm, os.tpm), "TMEM40")
pt.cex <- 0.3
shade <- "light2"
highlight.cex <- 1
highlight.shade <- "main"

# Plot
pdf(args$outfile, height=2.6, width=2.8)
par(mfrow=c(1, 2))
sdf <- swarmplot.by.phenotype(raf1, meta, title=bquote(italic(RAF1)),
                              y.axis.title=NULL, add.sample.size=TRUE,
                              shorten.cancer.names=TRUE, pt.cex=pt.cex,
                              shade=shade, return.swarm.df=TRUE,
                              parmar=c(3, 1, 1.5, 0.6))
points(sdf[dup.carriers, ], pch=23, cex=highlight.cex,
       bg=sapply(dup.carriers, function(sid){cancer.palettes[[metadata.cancer.label.map[meta[sid, "disease"]]]][highlight.shade]}))
sdf <- swarmplot.by.phenotype(tmem40, meta, title=bquote(italic(TMEM40)),
                              y.axis.title=NULL, add.sample.size=TRUE,
                              shorten.cancer.names=TRUE, pt.cex=pt.cex,
                              shade=shade, return.swarm.df=TRUE,
                              parmar=c(3, 1.5, 1.5, 0.4))
points(sdf[dup.carriers, ], pch=23, cex=highlight.cex,
       bg=sapply(dup.carriers, function(sid){cancer.palettes[[metadata.cancer.label.map[meta[sid, "disease"]]]][highlight.shade]}))
dev.off()

