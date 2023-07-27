#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Optimize list of samples to prune from PedSV callset basd on relatedness

options(stringsAsFactors=FALSE)

# Two simple positional arguments
args <- commandArgs(trailingOnly=T)

# # DEV:
# args <- c("~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/data/ancestry_and_relatedness/PedSV.v2.polished.kinship.tsv.gz",
#           "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_v2_callset_generation/v2_mega_batching/gatk_sv_pediatric_cancers_combined_sample_data_model_updated_links_6_2_23_with_annotations_for_batching.w_batch_assignments.tsv.gz")

# Load data
kin <- read.table(args[1], header=T, sep="\t", comment.char="")
meta <- read.table(args[2], header=T, sep="\t", comment.char="")

# Hard-exclude all 1000G samples and cancer trio parents
all.1kg <- meta[which(meta$study == "1000G"), 1]
cancer.parents <- meta[which(meta$study_phase == "trio" & meta$study != "1000G" & meta$proband == "No"), 1]
pruned <- c(all.1kg, cancer.parents)
kin <- kin[apply(kin, 1, function(v){all(!(v %in% pruned))}), ]

# Only bother with third degree relatives or closer
kin <- kin[which(kin$kin >= 1/16), ]
sids <- names(sort(table(c(kin$X.sample_1, kin$sample_2)), decreasing=T))

# Prioritize retaining cases over controls
sids <- sids[c(which(sids %in% meta[which(meta$disease == "control"), 1]),
             which(!(sids %in% meta[which(meta$disease == "control"), 1])))]

# Iterative pruning
for(sid in sids){
  # print(paste("Remaining pairs:", nrow(kin)))
  if(sid %in% c(kin$X.sample_1, kin$sample_2)){
    pruned <- c(pruned, sid)
    kin <- kin[-which(kin$X.sample_1 == sid | kin$sample_2 == sid), ]
  }
  # }else{
    # print(paste("Not found in remaining pairs:", sid))
  # }
  if(nrow(kin) == 0){
    break
  }
}

# Write to stdout
write(pruned, stdout())
