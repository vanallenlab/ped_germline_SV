#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2024-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Adjust noncoding CWAS Z-scores, P-values and effect sizes based on
# saddlepoint approximation of the null for inference


# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(PedSV, quietly=TRUE)
require(EQL, quietly=TRUE)


# Two simple positional arguments: in.tsv and out.tsv
args <- commandArgs(trailingOnly=TRUE)

# # DEV inputs:
# args <- c("~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/data/ncCWAS_sumstats/ewing_all_noncoding_cwas_concatenated_glm_results_12_26_23.txt",
#           "~/scratch/nccwas.correction.test.tsv")

# Read input
ss <- read.table(args[1], header=F, sep="\t")
colnames(ss) <- c("point_estimate", "std_error", "z_score", "p_value",
                  "SV_counts_cases", "SV_counts_cases_max",
                  "number_of_cases_with_zero_SVs", "total_cases",
                  "SV_counts_controls", "SV_counts_controls_max",
                  "number_of_controls_with_zero_SVs", "total_controls",
                  "number_of_unique_SVs", "category_name")

# Saddlepoint adjustment requires all tests to use the same test statistic
# Thus, we force everything to Z-scores
ss$z_score <- ss$point_estimate / ss$std_error

# Saddlepoint adjust Z-scores and P-values
spa.res <- saddlepoint.adj(ss$z_score)

# Transform effect sizes
spa.res$new_lnor <- spa.res$zscores * ss$std_error

# Update original ss data.frame
ss$point_estimate <- spa.res$new_lnor
ss$z_score <- spa.res$zscores
ss$p_value <- spa.res$pvalues

# Write to outfile
write.table(ss, args[2], col.names=F, row.names=F, sep="\t", quote=F)
