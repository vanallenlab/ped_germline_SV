#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Plot principal components for manuscript


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(PedSV, quietly=TRUE)
require(argparse, quietly=TRUE)
PedSV::load.constants("all")


##################
# Data functions #
##################
# Compile table of sample counts per ancestry per cohort
tabulate.ancestries <- function(meta){
  res <- lapply(c("ALL", sort(unique(meta$study_phase))), function(phase){
    subres <- lapply(c("ALL", sort(unique(meta$disease))), function(disease){
      do.call("rbind",
              lapply(c("ALL", unique(meta$inferred_ancestry)), function(pop){
                keep.idx <- 1:nrow(meta)
                if(phase != "ALL"){
                  keep.idx <- intersect(keep.idx, which(meta$study_phase == phase))
                }
                if(disease != "ALL"){
                  keep.idx <- intersect(keep.idx, which(meta$disease == disease))
                  disease <- metadata.cancer.label.map[disease]
                }
                if(pop != "ALL"){
                  keep.idx <- intersect(keep.idx, which(meta$inferred_ancestry == pop))
                }
                c(phase, disease, pop, length(keep.idx))
              }))
    })
    do.call("rbind", subres)
  })
  df <- as.data.frame(do.call("rbind", res))
  colnames(df) <- c("#study_phase", "disease", "ancestry", "N")
  return(df[which(as.numeric(df$N) > 0), ])
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Plot sample PCs")
parser$add_argument("metadata", metavar=".tsv", type="character",
                    help="sample metadata .tsv")
parser$add_argument("--subset-samples", metavar=".tsv", type="character",
                    help="list of samples to subset [default: use all samples]")
parser$add_argument("--out-prefix", metavar="path", type="character", required=TRUE,
                    help="path/prefix for all output files")
args <- parser$parse_args()

# # DEV:
# args <- list("metadata" = "~/scratch/PedSV.v2.5.3.cohort_metadata.w_control_assignments.tsv.gz",
#              "subset_samples" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_v2_callset_generation/v2.5.3/PedSV.v2.5.3.final_analysis_cohort.samples.list",
#              "out_prefix" = "~/scratch/PedSV.v2.5.3.dev")
# args <- list("metadata" = "/Users/ryan/Desktop/Collins/VanAllen/jackie_younglung/younglung_metadata/YL.SV.v1.1.analysis_metadata.tsv.gz",
#              "subset_samples" = "/Users/ryan/Desktop/Collins/VanAllen/jackie_younglung/YL_analysis/YL.analysis_samples.list",
#              "out_prefix" = "~/scratch/YL.SV.v1.1")

# Load PCs
keepers <- NULL
if(!is.null(args$subset_samples)){
  keepers <- read.table(args$subset_samples, header=F)[, 1]
}
meta <- PedSV::load.sample.metadata(args$metadata, keep.samples=keepers)

# Suffle order of samples to avoid control overplotting, and put controls first
set.seed(2023)
meta <- meta[sample(1:nrow(meta), nrow(meta), replace=F), ]
meta <- meta[c(which(meta$disease == "control"),
               which(meta$disease != "control")), ]

# Plot PCs
pdf.dim <- 3.25
axis.title.line <- 1.25
parmar <- c(3.25, 3.25, 0.5, 0.5)
pop.legend.labels <- sapply(names(pop.colors), function(pop){
  paste(pop.names.short[pop], " (",
        round(100 * length(which(meta$inferred_ancestry == pop)) / nrow(meta), 0),
        "%)", sep="")
})
disease.legend.labels <- sapply(names(cancer.colors), function(disease){
  paste(cancer.names.short[disease], " (",
        round(100 * length(which(metadata.cancer.label.map[meta$disease] == disease)) / nrow(meta), 0),
        "%)", sep="")
})
sapply(list(c(1, 2), c(3, 4), c(5, 6), c(7, 8)), function(pc.idxs){

  # Check for presence of both SNV- and SV-derived PCs
  # If present, plot & label them separately
  if(length(intersect(c("SV.PC1", "PC1"), colnames(meta))) == 2){
    pc.col.prefixes <- c("PC", "SV.PC")
    pc.label.prefixes <- c("SNV principal component", "SV principal component")
    plot.suffixes <- c(".SNV_PCs", "")
  }else{
    pc.col.prefixes <- c("PC")
    pc.label.prefixes <- "SV principal component"
    plot.suffixes <- ""
  }

  # Plot PCs
  for(pctype.idx in 1:length(pc.col.prefixes)){
    # Colored by ancestry
    pdf(paste(args$out_prefix, ".pc", pc.idxs[1], "_vs_pc",
              pc.idxs[2], ".by_ancestry", plot.suffixes[pctype.idx],
              ".pdf", sep=""),
        height=pdf.dim, width=pdf.dim)
    pc.scatterplot(meta,
                   paste(pc.col.prefixes[pctype.idx], pc.idxs[1], sep=""),
                   paste(pc.col.prefixes[pctype.idx], pc.idxs[2], sep=""),
                   colors=pop.colors[meta$inferred_ancestry],
                   x.title=paste(pc.label.prefixes[pctype.idx], pc.idxs[1]),
                   x.title.line=axis.title.line-0.5,
                   x.label.line=-0.7,
                   y.title=paste(pc.label.prefixes[pctype.idx], pc.idxs[2]),
                   y.title.line=axis.title.line,
                   y.label.line=-0.6,
                   legend.vals=pop.colors[which(names(pop.colors) != "OTH")],
                   legend.labels=pop.legend.labels[which(names(pop.colors) != "OTH")],
                   parmar=parmar)
    dev.off()

    # Colored by cancer type
    pdf(paste(args$out_prefix, ".pc", pc.idxs[1], "_vs_pc",
              pc.idxs[2], ".by_disease", plot.suffixes[pctype.idx],
              ".pdf", sep=""),
        height=pdf.dim, width=pdf.dim)
    pc.scatterplot(meta,
                   paste(pc.col.prefixes[pctype.idx], pc.idxs[1], sep=""),
                   paste(pc.col.prefixes[pctype.idx], pc.idxs[2], sep=""),
                   colors=cancer.colors[metadata.cancer.label.map[meta$disease]],
                   x.title=paste(pc.label.prefixes[pctype.idx], pc.idxs[1]),
                   x.title.line=axis.title.line-0.5,
                   x.label.line=-0.7,
                   y.title=paste(pc.label.prefixes[pctype.idx], pc.idxs[2]),
                   y.title.line=axis.title.line,
                   y.label.line=-0.6,
                   legend.vals=cancer.colors[which(names(cancer.colors) %in% unique(metadata.cancer.label.map[meta$disease]))],
                   legend.labels=disease.legend.labels[which(names(cancer.colors) %in% unique(metadata.cancer.label.map[meta$disease]))],
                   parmar=parmar)
    dev.off()
  }
})

# Compile table of sample counts per ancestry
write.table(tabulate.ancestries(meta),
            paste(args$out_prefix, "sample_counts_by_ancestry.tsv", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)
