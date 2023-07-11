#!/usr/bin/env R

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Define outlier samples in an ancestry-specific manner

# Note that this analysis code has not been polished or productionized
# It is intended for local, quick-and-dirty use
# Contact Ryan <Ryan_Collins@dfci.harvard.edu> with questions


require(beeswarm)

# Load data
x <- read.table("~/scratch/PedSV_v2_minGQ_postGQR_July_2023.rare.counts.tsv", header=T, sep="\t", comment.char="")
y <- read.table("~/scratch/PedSV_v2_minGQ_postGQR_July_2023.common.counts.tsv", header=T, sep="\t", comment.char="")
m <- read.table("/Users/ryan/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_v2_callset_generation/v2_mega_batching/gatk_sv_pediatric_cancers_combined_sample_data_model_updated_links_6_2_23_with_annotations_for_batching.w_batch_assignments.tsv.gz",
                header=T, sep="\t", comment.char="")
rownames(m) <- m$entity.sample_id

# Update ancestry labels with v1 labels
om <- PedSV::load.sample.metadata("~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/data/gatk_sv_pediatric_cancers_combined_cohort_metadata_3_31_23.txt")
overlaps <- m$entity.sample_id[which(m$entity.sample_id %in% rownames(om))]
m[overlaps, "reported_or_SNV_ancestry"] <- om[overlaps, "inferred_ancestry"]

# Merge count & metadata
x <- merge(x, m, by.x="sample", by.y="entity.sample_id", sort=F, all.y=F, all.x=T)
RASMod::load.constants("colors")
pop.colors["OTH"] <- "gray70"
# x$reported_or_SNV_ancestry[which(is.na(x$reported_or_SNV_ancestry))] <- "OTH"
y <- merge(y, m, by.x="sample", by.y="entity.sample_id", sort=F, all.y=F, all.x=T)
# y$reported_or_SNV_ancestry[which(is.na(y$reported_or_SNV_ancestry))] <- "OTH"
outlier.pch <- c("FALSE" = 19, "TRUE" = 4)
outlier.cex <- c("FALSE" = 0.15, "TRUE" = 0.6)

# Define outliers
x$OUTLIER <- FALSE
for(svtype in c("DUP", "INS", "DEL", "CPX", "BND")){
  for(pop in unique(x$reported_or_SNV_ancestry)){
    idx <- which(x$svtype == svtype & x$reported_or_SNV_ancestry == pop)
    vals <- x[idx, "count"]
    q3 <- quantile(vals, prob=0.75)
    MAD <- mad(vals)
    cutoff <- q3 + (6*MAD)
    x$OUTLIER[idx[which(vals > cutoff)]] <- TRUE
  }
}

# Plot SVTYPE x POP
png("~/scratch/PedSV.v2.rare_counts.by_pop.png", height=2.4*300, width=8*300, res=300)
par(mfrow=c(1, 6), mar=c(0.1, 2, 1, 1))
sapply(setdiff(unique(x$svtype), c("CTX")), function(svtype){
  idxs <- which(x$svtype == svtype)
  boxplot(count ~ reported_or_SNV_ancestry, data=x[idxs, ], main=svtype, xlab="", ylab="",
          xaxt="n", names=NA, outline=F, lty=1, col="gray85")
  beeswarm(count ~ reported_or_SNV_ancestry, data=x[idxs, ], add=T,
           xlab="", ylab="", xaxt="n", yaxt="n",
           pwpch=outlier.pch[as.character(x[idxs, "OUTLIER"])],
           pwcex=outlier.cex[as.character(x[idxs, "OUTLIER"])],
           corral="wrap", corralWidth=0.7, lwd=0.75,
           pwcol=pop.colors[x$reported_or_SNV_ancestry[idxs]])
})
dev.off()

# Plot COHORT x POP
png("~/scratch/PedSV.v2.rare_counts.by_cohort.png", height=3*300, width=10*300, res=300)
par(mfrow=c(1, 6), mar=c(4, 2, 1, 1))
sapply(setdiff(unique(x$svtype), c("CTX")), function(svtype){
  idxs <- which(x$svtype == svtype)
  boxplot(count ~ study, data=x[idxs, ], main=svtype, xlab="", ylab="",
          outline=F, lty=1, col="gray85", las=2)
  beeswarm(count ~ study, data=x[idxs, ], add=T,
           xlab="", ylab="", xaxt="n", yaxt="n",
           pch=19, cex=0.15, corral="wrap", corralWidth=0.7, lwd=0.75,
           pwcol=pop.colors[x$reported_or_SNV_ancestry[idxs]])
})
dev.off()

# Plot common by ancestry, for comparison
png("~/scratch/PedSV.v2.common_counts.by_pop.png", height=2.4*300, width=10*300, res=300)
par(mfrow=c(1, 7), mar=c(0.1, 2, 1, 1))
sapply(unique(y$svtype), function(svtype){
  idxs <- which(y$svtype == svtype)
  boxplot(count ~ reported_or_SNV_ancestry, data=y[idxs, ], main=svtype, xlab="", ylab="",
          xaxt="n", names=NA, outline=F, lty=1, col="gray85")
  beeswarm(count ~ reported_or_SNV_ancestry, data=y[idxs, ], add=T,
           xlab="", ylab="", xaxt="n", yaxt="n",
           pch=19, cex=0.15, corral="wrap", corralWidth=0.7, lwd=0.75,
           pwcol=pop.colors[y$reported_or_SNV_ancestry[idxs]])
})
dev.off()

# Plot common by ancestry, for comparison
png("~/scratch/PedSV.v2.common_counts.by_study.png", height=3*300, width=10*300, res=300)
par(mfrow=c(1, 7), mar=c(4, 2, 1, 1))
sapply(unique(y$svtype), function(svtype){
  idxs <- which(y$svtype == svtype)
  boxplot(count ~ study, data=y[idxs, ], main=svtype, xlab="", ylab="",
          outline=F, lty=1, col="gray85", las=2)
  beeswarm(count ~ study, data=y[idxs, ], add=T,
           xlab="", ylab="", xaxt="n", yaxt="n",
           pch=19, cex=0.15, corral="wrap", corralWidth=0.7, lwd=0.75,
           pwcol=pop.colors[y$reported_or_SNV_ancestry[idxs]])
})
dev.off()

# Write list of samples
bad.ids <- unique(x$sample[which(x$OUTLIER)])
write.table(bad.ids, "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_v2_callset_generation/PedSV.v2.rare_sv_outliers.list",
            col.names=F, row.names=F, quote=F)
