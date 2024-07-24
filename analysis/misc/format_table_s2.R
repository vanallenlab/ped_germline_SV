#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2024-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Format Table S2 for PedSV manuscript


#########
# Setup #
#########
# Load necessary libraries and constants
require(PedSV, quietly=TRUE)
options(scipen=1000, stringsAsFactors=F)

# Simple command-line arguments: in.tsv, out.tsv
args <- commandArgs(trailingOnly=TRUE)

# Read burden summary stats
ss <- read.table(args[1], header=T, sep="\t", comment.char="", check.names=F)

# Subset to categories of interest
hyp.keep <- c("large", "large.rare", "large.rare.unbalanced", "large.rare.unbalanced.MALE_only",
              "large.rare.unbalanced.FEMALE_only", "large.rare.unbalanced.autosomal_only",
              "large.rare.unbalanced.autosomal_only.MALE_only",
              "large.rare.unbalanced.autosomal_only.FEMALE_only",
              "large.rare.unbalanced.allosomal_only", "large.rare.balanced",
              "large.singleton", "large.singleton.unbalanced",
              "large.singleton.unbalanced.MALE_only", "large.singleton.unbalanced.FEMALE_only",
              "large.singleton.unbalanced.autosomal_only",
              "large.singleton.unbalanced.autosomal_only.MALE_only",
              "large.singleton.unbalanced.autosomal_only.FEMALE_only",
              "large.singleton.unbalanced.allosomal_only", "large.singleton.balanced",
              "large.rare.unbalanced.autosomal_only.EUR_only",
              "large.rare.unbalanced.autosomal_only.EUR_only.MALE_only",
              "large.rare.unbalanced.autosomal_only.EUR_only.FEMALE_only",
              "large.rare.unbalanced.autosomal_only.no_COSMIC_no_CPG",
              "large.rare.unbalanced.autosomal_only.no_COSMIC_no_CPG.MALE_only",
              "large.rare.unbalanced.autosomal_only.no_COSMIC_no_CPG.FEMALE_only",
              "all_genes_disrupted", "rare.gene_disruptive",
              "rare.gene_disruptive.rare_LoF_SVs", "rare.gene_disruptive.rare_CG_SVs",
              "rare.gene_disruptive.MALE_only", "rare.gene_disruptive.FEMALE_only",
              "rare.gene_disruptive.no_large_unbalanced", "singleton.gene_disruptive",
              "singleton.gene_disruptive.singleton_LoF_SVs",
              "singleton.gene_disruptive.singleton_CG_SVs",
              "singleton.gene_disruptive.singleton_IED_SVs", "singleton.gene_disruptive.MALE_only",
              "singleton.gene_disruptive.FEMALE_only",
              "singleton.gene_disruptive.no_large_unbalanced", "rare.single_gene_disruptive",
              "rare.single_gene_disruptive.rare_single_gene_LoF_SVs",
              "rare.single_gene_disruptive.rare_single_gene_CG_SVs",
              "rare.single_gene_disruptive.rare_single_gene_IED_SVs",
              "singleton.single_gene_disruptive",
              "singleton.single_gene_disruptive.singleton_single_gene_LoF_SVs",
              "singleton.single_gene_disruptive.singleton_single_gene_CG_SVs",
              "rare.cosmic.gene_disruptive", "rare.cosmic.gene_disruptive.rare_cosmic_LoF",
              "rare.cosmic.gene_disruptive.rare_cosmic_CG", "rare.cpg.gene_disruptive",
              "rare.cpg.gene_disruptive.rare_cpg_LoF", "rare.cpg.gene_disruptive.rare_cpg_CG",
              "singleton.cosmic.gene_disruptive",
              "singleton.cosmic.gene_disruptive.singleton_cosmic_LoF",
              "singleton.cosmic.gene_disruptive.singleton_cosmic_CG",
              "singleton.cpg.gene_disruptive", "singleton.cpg.gene_disruptive.singleton_cpg_LoF",
              "singleton.cpg.gene_disruptive.singleton_cpg_CG")
keep.cols <- c("#hypothesis", "disease", "n.case", "n.case.carrier", "case.mean",
               "case.stdev", "n.control", "n.control.carrier", "control.mean",
               "control.stdev", "coefficient", "std.err", "test.statistic", "P.value")
ss <- ss[which(ss$`#hypothesis` %in% hyp.keep), keep.cols]

# Write to out.tsv
write.table(ss, args[2], col.names=T, row.names=F, quote=F, sep="\t")

