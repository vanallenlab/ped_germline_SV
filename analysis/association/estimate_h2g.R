#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Estimate heritability from rare SVs for a single phenotype


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(glmnet, quietly=TRUE)
require(caret, quietly=TRUE)
require(drcarlate, quietly=TRUE)
require(PedSV, quietly=TRUE)
require(argparse, quietly=TRUE)
PedSV::load.constants("all")


##################
# Data functions #
##################
# Load feature matrixes for coding and noncoding CWAS and subset to relevant samples
load.cwas.features <- function(coding.in, noncoding.in, meta){
  # Read data
  coding.df <- as.data.frame(t(read.table(coding.in, header=T, sep="\t", row.names=1)))
  noncoding.df <- as.data.frame(t(read.table(noncoding.in, header=T, sep="\t", row.names=1)))

  # Merge data & subset to samples present in `meta`
  m.df <- merge(coding.df, noncoding.df, by="row.names", sort=F)
  rownames(m.df) <- m.df$Row.names
  m.df$Row.names <- NULL
  cat(paste("Loaded", prettyNum(ncol(m.df), big.mark=","), "CWAS categories\n"))
  m.df[rownames(meta), ]
}

# Reduce dimensions of a feature matrix with LASSO regression
select.features <- function(cwas.features.all, Y){
  Y <- Y[rownames(cwas.features.all)]
  X <- makeX(cwas.features.all, na.impute=TRUE)

  # Optimize LASSO penalty following example from:
  # https://www.statology.org/lasso-regression-in-r/
  cv_model <- cv.glmnet(X, Y, family="binomial", alpha=1)
  best_lambda <- cv_model$lambda.min
  best_model <- glmnet(X, Y, family="binomial", alpha=1, lambda=best_lambda)

  # Extract feature names of non-zero coefficients
  nonzero.features <- colnames(X)[summary(best_model$beta)$i]
  cat(paste("Retained", prettyNum(length(nonzero.features), big.mark=","),
            "CWAS categories after LASSO regression\n"))
  X[, nonzero.features]
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Estimate heritability attributable to rare SVs")
parser$add_argument("--bed", metavar=".tsv", type="character", required=TRUE,
                    help="SV sites .bed.")
parser$add_argument("--ad", metavar=".tsv", type="character", required=TRUE,
                    help="Allele dosage .bed.")
parser$add_argument("--metadata", metavar=".tsv", type="character",
                    help="sample metadata .tsv", required=TRUE)
parser$add_argument("--cancer", metavar="string", type="character", required=TRUE,
                    help="Cancer type to model [EWS, NBL, or OS]")
parser$add_argument("--coding-features", metavar=".tsv", type="character", required=TRUE,
                    help="Matrix of categories by samples from coding CWAS")
parser$add_argument("--noncoding-features", metavar=".tsv", type="character", required=TRUE,
                    help="Matrix of categories by samples from noncoding CWAS")
parser$add_argument("--prevalence", metavar="float", type="numeric", required=TRUE,
                    help="Estimated population prevalence of --cancer")
parser$add_argument("--subset-samples", metavar=".tsv", type="character",
                    help="list of samples to subset [default: use all samples]")
parser$add_argument("--out-tsv", metavar="path", type="character", required=TRUE,
                    help="path to output .tsv with final feature coefficients")
args <- parser$parse_args()

# # DEV:
# args <- list("bed" = "~/scratch/PedSV.v2.5.4.full_cohort.analysis_samples.sites.bed.gz",
#              "ad" = "~/scratch/PedSV.v2.5.4.full_cohort.analysis_samples.allele_dosages.bed.gz",
#              "metadata" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_v2_callset_generation/v2.5.4/PedSV.v2.5.4.cohort_metadata.w_control_assignments.tsv.gz",
#              "cancer" = "NBL",
#              "coding_features" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/data/heritability_inputs/neuroblastoma_all_coding_CWAS_concatenated_SV_count_results_formatted_nominally_significant_categories_8_9_24.txt",
#              "noncoding_features" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/data/heritability_inputs/neuroblastoma_all_noncoding_CWAS_concatenated_SV_count_results_formatted_nominally_significant_categories_8_9_24.txt",
#              "prevalence" = 0.0001428571,
#              "subset_samples" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_v2_callset_generation/v2.5.4/PedSV.v2.5.4.final_analysis_cohort.samples.list",
#              "out_tsv" = "~/scratch/PedSV.v2.5.4.full_cohort.h2g_fit.coeffs.tsv")

# Print startup to log
cat(paste("\n\nHeritability estimation for ",
          tolower(cancer.names.long[args$cancer]),
          ":\n", sep=""))

# Load metadata
keepers <- NULL
if(!is.null(args$subset_samples)){
  keepers <- read.table(args$subset_samples, header=F)[, 1]
}
meta <- load.sample.metadata(args$metadata, keep.samples=keepers,
                             reassign.parents=FALSE)

# Subset to cases & controls appropriate for cancer type of interest
elig.ids <- get.eligible.samples(meta, args$cancer)
meta <- meta[intersect(rownames(meta), unlist(elig.ids)), ]
cat(paste("Retained", prettyNum(length(elig.ids$cases), big.mark=","),
          args$cancer, "cases and",
          prettyNum(length(elig.ids$controls), big.mark=","),
          "ancestry-matched controls for analysis\n"))

# Fit baseline heritability model based on demographic data alone
Y <- get.phenotype.vector(elig.ids$cases, elig.ids$controls)
meta$Y <- Y[rownames(meta)]
test.df <- meta[, c("Y", "inferred_sex", "inferred_ancestry")]
baseline.fit <- summary(train(Y ~ ., method="lm", data=test.df, trControl=trainControl(method="cv", number=10)))
cat(paste("Baseline model of Y ~ sex + ancestry:\n  Multiple R2 =",
          baseline.fit$r.squared, "\n  Adjusted R2 =",
          baseline.fit$adj.r.squared, "\n"))

# Load CWAS data and reorient according to regression model
cwas.features.all <- load.cwas.features(args$coding_features, args$noncoding_features, meta)

# Feature selection on CWAS data with LASSO penalty
cwas.features.subset <- select.features(cwas.features.all, Y)

# Get large, rare, unbalanced SV counts
lru.bed <- filter.bed(load.sv.bed(args$bed), query="large.rare.unbalanced")
lru.k <- query.ad.from.sv.bed(args$ad, lru.bed, action="any")
test.df$lru <- lru.k[rownames(test.df)]
test.df$lruXsex <- test.df$lru * as.integer(round(meta$chrX_CopyNumber) <= 1)

# Merge selected CWAS features with other features
test.df <- as.data.frame(cbind(test.df[names(Y), ],
                               cwas.features.subset[names(Y), ]))

# Fit heritability model with SVs
test.df <- impute.missing.values(test.df)
h2g.fit <- summary(train(Y ~ ., method="lm", data=test.df, trControl=trainControl(method="cv", number=10)))
cat(paste("Extended model of Y ~ SVs + sex + ancestry:\n  Multiple R2 =",
          h2g.fit$r.squared, "\n  Adjusted R2 =",
          h2g.fit$adj.r.squared, "\n"))
if(!is.null(args$out_tsv)){
  out.df <- as.data.frame(h2g.fit$coefficients)
  out.df$feature <- rownames(out.df)
  write.table(out.df[, c("feature", setdiff(colnames(out.df), "feature"))],
              args$out_tsv, col.names=T, row.names=F, sep="\t", quote=F)
}

# Transform R2 to liability scale per Equation 23 in Lee et al., AJHG, 2011
K <- args$prevalence
P <- sum(Y) / length(Y)
h2l.multi <- h2g.fit$r.squared * ( (K * (1-K)) / (dnorm(norminv(K))^2) ) * ( (K * (1-K)) / (P * (1-P)) )
h2l.adj <- h2g.fit$adj.r.squared * ( (K * (1-K)) / (dnorm(norminv(K))^2) ) * ( (K * (1-K)) / (P * (1-P)) )
cat(paste("Liability-scale heritability explained by rare SVs: ",
          round(100 * h2l.adj, 1), "-", round(100 * h2l.multi, 1), "%\n", sep=""))

