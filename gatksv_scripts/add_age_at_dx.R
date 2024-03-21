#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2024-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Add age-at-diagnosis information to main sample metadata manifest


# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)

# Two simple positional arguments: input metadata, age-at-dx .tsv, and output metadata
args <- commandArgs(trailingOnly=TRUE)

# Read input metadata as simple df
meta <- read.table(args[1], sep="\t", check.names=F, header=T, comment.char="")
rownames(meta) <- meta[, 1]
orig.order <- meta[, 1]

# Read ages
ages <- read.table(args[2], sep="\t", header=T)
colnames(ages) <- c("entity:sample_id", "age_at_dx_years")

# Deduplicate age data, taking the earliest diagnosis age in the case of relapses
ages <- ages[order(ages[, 2]), ]
ages <- ages[which(!duplicated(ages[, 1])), ]

# Append age information to metadata
meta <- merge(meta, ages, by="entity:sample_id", sort=F, all.x=T, all.y=F)
rownames(meta) <- meta[, 1]

# Write updated metadata to disk
write.table(meta[orig.order, ], args[3], row.names=F, col.names=T, sep="\t", quote=F)
