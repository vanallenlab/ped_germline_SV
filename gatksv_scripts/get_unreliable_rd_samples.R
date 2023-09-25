#!/usr/bin/env R

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Define samples with unreliable RD GTs based on rates of rare RD-only CNVs in an ancestry-specific manner

# Note that this analysis code has not been polished or productionized
# It is intended for local, quick-and-dirty use
# Contact Ryan <Ryan_Collins@dfci.harvard.edu> with questions


require(beeswarm)

# Load data
x <- read.table("~/scratch/PedSV.v2.1.minGQ_postGQR.noOutliers.bfx_marked.rare_RD_only.counts.tsv", header=T, sep="\t", comment.char="")
x$svtype <- paste("RD-Only", x$svtype)
elig.samples <- read.table("~/scratch/PedSV.v2.1.minGQ_postGQR.noOutliers.bfx_marked.rare_RD_only.all_samples.list", header=F, comment.char="")[, 1]
m <- read.table("/Users/ryan/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_v2_callset_generation/v2_mega_batching/gatk_sv_pediatric_cancers_combined_sample_data_model_updated_links_6_2_23_with_annotations_for_batching.w_batch_assignments.tsv.gz",
                header=T, sep="\t", comment.char="")
rownames(m) <- m$entity.sample_id
m <- m[elig.samples, ]

# Update ancestry labels with v1 labels
om <- PedSV::load.sample.metadata("~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/data/gatk_sv_pediatric_cancers_combined_cohort_metadata_3_31_23.txt")
overlaps <- m$entity.sample_id[which(m$entity.sample_id %in% rownames(om))]
m[overlaps, "reported_or_SNV_ancestry"] <- om[overlaps, "inferred_ancestry"]
m <- rbind(cbind(m, "svtype"="RD-Only DEL"),
           cbind(m, "svtype"="RD-Only DUP"))

# Merge count & metadata
x <- merge(x, m, by.x=c("sample", "svtype"), by.y=c("entity.sample_id", "svtype"), sort=F, all.y=T, all.x=T)
x$count[which(is.na(x$count))] <- 0
RASMod::load.constants("colors")
x$reported_or_SNV_ancestry[which(is.na(x$reported_or_SNV_ancestry))] <- "EUR"
# y <- merge(y, m, by.x="sample", by.y="entity.sample_id", sort=F, all.y=F, all.x=T)
# y$reported_or_SNV_ancestry[which(is.na(y$reported_or_SNV_ancestry))] <- "OTH"
outlier.pch <- c("FALSE" = 19, "TRUE" = 4)
outlier.cex <- c("FALSE" = 0.15, "TRUE" = 0.6)

# Define outliers
x$OUTLIER <- FALSE
for(svtype in unique(x$svtype)){
  for(pop in unique(x$reported_or_SNV_ancestry)){
    idx <- which(x$svtype == svtype & x$reported_or_SNV_ancestry == pop)
    vals <- x[idx, "count"]
    q3 <- quantile(vals, prob=0.75)
    MAD <- mad(vals)
    cutoff <- q3 + (6*MAD)
    x$OUTLIER[idx[which(vals > cutoff)]] <- TRUE
  }
}

x$count <- log10(x$count + 1)

# Plot SVTYPE x POP
png("~/scratch/PedSV.v2.rare_RD_only_counts.by_pop.png", height=2.5*300, width=5*300, res=300)
par(mfrow=c(1, 2), mar=c(0.1, 4, 1, 1))
sapply(setdiff(unique(x$svtype), c("CTX")), function(svtype){
  idxs <- which(x$svtype == svtype)
  boxplot(count ~ reported_or_SNV_ancestry, data=x[idxs, ], main=svtype, xlab="",
          xaxt="n", names=NA, outline=F, lty=1, col="gray85", ylim=c(0, max(x$count[idxs])),
          ylab=bquote(log[10]("SV Count" + 1)))
  beeswarm(count ~ reported_or_SNV_ancestry, data=x[idxs, ], add=T,
           xlab="", ylab="", xaxt="n", yaxt="n",
           pwpch=outlier.pch[as.character(x[idxs, "OUTLIER"])],
           pwcex=outlier.cex[as.character(x[idxs, "OUTLIER"])],
           corral="wrap", corralWidth=0.7, lwd=0.75,
           pwcol=pop.colors[x$reported_or_SNV_ancestry[idxs]])
})
dev.off()

# Plot COHORT x POP
png("~/scratch/PedSV.v2.rare_RD_only_counts.by_cohort.png", height=3*300, width=6*300, res=300)
par(mfrow=c(1, 2), mar=c(4, 4, 1, 1))
sapply(setdiff(unique(x$svtype), c("CTX")), function(svtype){
  idxs <- which(x$svtype == svtype)
  boxplot(count ~ study, data=x[idxs, ], main=svtype, xlab="",
          outline=F, lty=1, col="gray85", las=2, ylim=c(0, max(x$count[idxs])),
          ylab=bquote(log[10]("SV Count" + 1)))
  beeswarm(count ~ study, data=x[idxs, ], add=T,
           xlab="", ylab="", xaxt="n", yaxt="n",
           corral="wrap", corralWidth=0.7, lwd=0.75,
           pwcol=pop.colors[x$reported_or_SNV_ancestry[idxs]],
           pwpch=outlier.pch[as.character(x[idxs, "OUTLIER"])],
           pwcex=outlier.cex[as.character(x[idxs, "OUTLIER"])],)
})
dev.off()

# Write list of samples
bad.ids <- unique(x$sample[which(x$OUTLIER)])
write.table(bad.ids, "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_v2_callset_generation/PedSV.v2.1.rare_RD_only_outliers.samples.list",
            col.names=F, row.names=F, quote=F)
