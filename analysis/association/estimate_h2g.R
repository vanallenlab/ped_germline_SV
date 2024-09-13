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
  impute.missing.values(m.df[rownames(meta), ])
}

# Optimize elastic net penalty (lambda) with 10-fold CV following example from:
# https://www.statology.org/lasso-regression-in-r/
optimize.penalty <- function(X, Y, alpha=0.5, k=10){
  set.seed(1234)
  cv_model <- cv.glmnet(X, Y, alpha=alpha, nfolds=k)
  cv_model$lambda.min
}

# Reduce dimensions of a feature matrix with LASSO regression or simple feature correlation
select.features <- function(cwas.features.all, Y=NULL, r2.cutoff=0.8, k=10){
  # If Y is not specified, no phenotype-guided feature selection will be used
  if(is.null(Y)){
    # Compute feature correlation matrix
    feature.r2 <- cor(cwas.features.all, use="complete.obs")^2
    cor.pairs <- which(feature.r2 >= r2.cutoff, arr.ind=TRUE)
    rownames(cor.pairs) <- NULL
    cor.pairs <- t(apply(cor.pairs, 1, function(idxs){sort(c(rownames(feature.r2)[idxs[1]],
                                      colnames(feature.r2)[idxs[2]]))}))
    cor.pairs <- cor.pairs[which(cor.pairs[, 1] != cor.pairs[, 2]), ]
    cor.pairs <- cor.pairs[!duplicated(cor.pairs), ]
    drop.candidates <- unique(as.vector(unlist(cor.pairs)))
    f.var <- apply(cwas.features.all, 2, var, na.rm=T)
    f.mean <- apply(cwas.features.all, 2, mean, na.rm=T)

    # Iteratively prune correlated features, until none remain with r2 > r2.cutoff
    # Preferentially retain features with greater variance & mean to break ties
    features.dropped <- c()
    while(length(cor.pairs) > 0){
      if(is.vector(cor.pairs) & length(cor.pairs) == 2){
        cor.pairs <- matrix(cor.pairs, ncol=2, byrow=T)
      }
      # Recompute drop priority
      n.cor <- sapply(drop.candidates, function(k){length(which(unlist(cor.pairs) == k))})
      drop.priority.df <- data.frame("n.cor"=n.cor,
                                     "f.var"=f.var[drop.candidates],
                                     "f.mean"=f.mean[drop.candidates],
                                     row.names=drop.candidates)
      drop.priority <- with(drop.priority.df, order(-n.cor, f.var, f.mean))

      # Drop one feature
      drop.feature <- rownames(drop.priority.df)[drop.priority][1]
      features.dropped <- c(features.dropped, drop.feature)
      drop.candidates <- setdiff(drop.candidates, drop.feature)
      cor.pairs <- cor.pairs[apply(cor.pairs, 1, function(p){!(drop.feature %in% p)}), ]
    }

    # Return decorrelated features
    keep.features <- setdiff(colnames(cwas.features.all), features.dropped)
    cat(paste("Retained", prettyNum(length(keep.features), big.mark=","),
              "CWAS categories after simple decorrelation\n"))
    cwas.features.all[, keep.features]

  }else{
    Y <- Y[rownames(cwas.features.all)]
    X <- makeX(cwas.features.all, na.impute=TRUE)

    # Optimize LASSO penalty and fit best model
    best.lambda <- optimize.penalty(X, Y, alpha=1, k=k)
    best.model <- glmnet(X, Y, alpha=1, lambda=best.lambda)

    # Extract feature names of non-zero coefficients
    nonzero.features <- colnames(X)[summary(best.model$beta)$i]
    cat(paste("Retained", prettyNum(length(nonzero.features), big.mark=","),
              "CWAS categories after LASSO regression\n"))
    X[, nonzero.features]
  }
}

# Estimate heritability with ridge regression and CV
h2.ridge <- function(test.df, k=10, alpha=0){
  Y <- test.df$Y
  X <- makeX(test.df[, -which(colnames(test.df) == "Y")], na.impute=TRUE)

  # Next, estimate heritability using 10-fold CV to protect against overfitting
  cv.fits <- list()
  set.seed(2024)
  test.idxs <- sample(1:k, size=nrow(X), replace=T)
  for(i in 1:k){
    # Fit model on 90% of data
    X.train <- X[which(test.idxs != i), ]
    Y.train <- Y[which(test.idxs != i)]
    X.test <- X[which(test.idxs == i), ]
    Y.test <- Y[which(test.idxs == i)]
    best.lambda <- optimize.penalty(X.train, Y.train, alpha=alpha, k=k)
    set.seed(i)
    fit <- glmnet(X.train, Y.train, alpha=alpha, lambda=best.lambda)

    # Evaluate (RMSE) on 10% of data
    Y.pred <- predict(fit, newx=X.test)
    rmse <- RMSE(Y.pred, Y.test)
    r2 <- cor(Y.pred, Y.test)^2

    # Add model fit and RMSE to collector
    cv.fits[[i]] <- list("fit"=fit, "rmse"=rmse, "r.squared"=r2)
  }

  # Select best model as lowest test RMSE
  best.idx <- which.min(as.numeric(unlist(lapply(cv.fits, function(l){l$rmse}))))
  best.model <- cv.fits[[best.idx]]$fit

  # Print range of test R2 values for comparison
  test.r2s <- sapply(cv.fits, function(l){l$r.squared})

  cat(paste("  Test R2 range =", paste(sort(round(test.r2s, 5)), collapse=", "), "\n"))
  cat(paste("  CV-optimal model R2 =", test.r2s[best.idx], "\n"))

  # Make predictions on full dataset and return R2 and betas
  Y.pred <- predict(best.model, newx=X)
  r2 <- as.numeric(cor(Y, Y.pred)^2)
  betas <- summary(best.model$beta)
  rownames(betas) <- colnames(X)[betas$i]
  betas[, c("i", "j")] <- NULL
  colnames(betas) <- "beta"
  return(list("r.squared" = r2, "betas" = betas, "test.r.squared" = test.r2s))
}

# Transform R2 to liability scale per Equation 23 in Lee et al., AJHG, 2011
lee.transform <- function(r2, K, P){
  r2 * ( (K * (1-K)) / (dnorm(norminv(K))^2) ) * ( (K * (1-K)) / (P * (1-P)) )
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
parser$add_argument("--feature-selection", default=FALSE, action="store_true",
                    help=paste("should CWAS dimensionality be reduced with LASSO",
                               "vs. case|control status prior to h2 estimation?",
                               "[default: decorrelate features agnostic to",
                               "case|control status]"))
parser$add_argument("--unpenalized-final-estimate", default=FALSE, action="store_true",
                    help=paste("do not penalize final h2 estimation [default:",
                               "penalize final estimate with ridge regression]"))
parser$add_argument("--out-tsv", metavar="path", type="character", required=TRUE,
                    help="path to output .tsv with final feature coefficients")
args <- parser$parse_args()

# # DEV:
# args <- list("bed" = "~/scratch/PedSV.v2.5.4.full_cohort.analysis_samples.sites.bed.gz",
#              "ad" = "~/scratch/PedSV.v2.5.4.full_cohort.analysis_samples.allele_dosages.bed.gz",
#              "metadata" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_v2_callset_generation/v2.5.4/PedSV.v2.5.4.cohort_metadata.w_control_assignments.tsv.gz",
#              "cancer" = "EWS",
#              "coding_features" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/data/heritability_inputs/ewing_all_coding_CWAS_concatenated_SV_count_results_formatted_8_9_24.txt",
#              "noncoding_features" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/data/heritability_inputs/ewing_all_noncoding_CWAS_concatenated_SV_count_results_formatted_8_9_24.txt",
#              "prevalence" = 0.0001428571,
#              "subset_samples" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/PedSV_v2_callset_generation/v2.5.4/PedSV.v2.5.4.final_analysis_cohort.samples.list",
#              "feature_selection" = FALSE,
#              "unpenalized_final_estimate" = FALSE,
#              "out_tsv" = "~/scratch/PedSV.v2.5.4.full_cohort.EWS.h2g_fit.coeffs.tsv")

# Print startup to log
cat(paste("\n\nHeritability estimation for ",
          tolower(cancer.names.long[args$cancer]),
          ":\n", sep=""))

# Load metadata and restrict to European samples
keepers <- NULL
if(!is.null(args$subset_samples)){
  keepers <- read.table(args$subset_samples, header=F)[, 1]
}
meta <- load.sample.metadata(args$metadata, keep.samples=keepers,
                             reassign.parents=FALSE)
meta <- meta[which(meta$inferred_ancestry == "EUR"), ]

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
# test.df <- meta[, c("Y", "inferred_sex", "inferred_ancestry")]
test.df <- meta[, c("Y", "inferred_sex")]
baseline.fit <- summary(train(Y ~ ., method="lm", data=test.df,
                              trControl=trainControl(method="cv", number=5)))
cat(paste("Baseline model of Y ~ sex:\n  Multiple R2 =",
          baseline.fit$r.squared, "\n  Adjusted R2 =",
          baseline.fit$adj.r.squared, "\n"))

# Load CWAS data and reorient according to regression model
cwas.features.all <- load.cwas.features(args$coding_features, args$noncoding_features, meta)

# Feature selection on CWAS data with LASSO penalty (if optioned) or simple decorrelation
if(args$feature_selection){
  cwas.features.subset <- select.features(cwas.features.all, Y, k=5)
}else{
  cwas.features.subset <- select.features(cwas.features.all, Y=NULL)
}

# Get large, rare, unbalanced SV counts
lru.bed <- filter.bed(load.sv.bed(args$bed), query="large.rare.unbalanced")
lru.k <- query.ad.from.sv.bed(args$ad, lru.bed, action="any")
test.df$lru <- lru.k[rownames(test.df)]
test.df$lruXsex <- test.df$lru * as.integer(round(meta$chrX_CopyNumber) <= 1)

# Merge selected CWAS features with other features
test.df <- as.data.frame(cbind(test.df[names(Y), ],
                               cwas.features.subset[names(Y), ]))
test.df <- impute.missing.values(test.df)

# Fit heritability model with SVs
cat("Extended model of Y ~ SVs + sex:\n")
if(args$unpenalized_final_estimate){
  h2g.fit <- summary(train(Y ~ ., method="lm", data=test.df,
                           trControl=trainControl(method="cv", number=5)))
  adj.r2 <- h2g.fit$adj.r.squared
  out.df <- as.data.frame(h2g.fit$coefficients)
  out.df$feature <- rownames(out.df)
  out.df <- out.df[, c("feature", setdiff(colnames(out.df), "feature"))]
}else{
  h2g.fit <- h2.ridge(test.df, k=5)
  adj.r2 <- NULL
  out.df <- h2g.fit$betas
}
multi.r2 <- h2g.fit$r.squared

# Report raw heritability estimate
cat(paste("  Multiple R2 =", multi.r2, "\n  Adjusted R2 =", adj.r2, "\n"))

# Save fitted coefficients to file, if optioned
if(!is.null(args$out_tsv)){
  write.table(out.df, args$out_tsv, col.names=T, row.names=F, sep="\t", quote=F)
}

# Report heritability estimates on the liability scale
K <- args$prevalence
P <- sum(Y) / length(Y)
cat("Liability-scale heritability explained by rare SVs:\n")
cat(paste("  Full cohort in-sample, multi R2:",
          round(lee.transform(h2g.fit$r.squared, K, P), 5), "\n"))
if(!is.null(h2g.fit$adj.r.squared)){
  cat(paste("  Full cohort in-sample, adjusted R2:",
            round(lee.transform(h2g.fit$adj.r.squared, K, P), 5), "\n"))
}
if(!is.null(h2g.fit$test.r.squared)){
  cat("  Held-out test fold R2: ")
  for(r2 in sort(h2g.fit$test.r.squared)){
    cat(paste(round(lee.transform(r2, K, P), 5), ","))
  }
  cat("\n")
}

