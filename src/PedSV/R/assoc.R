#!/usr/bin/env R

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Statistical and helper functions for association testing


#' Construct association dataframe
#'
#' Prepare a dataframe with dependent and independent variables for association
#' testing with additional covariates as optioned
#'
#' @param meta Sample metadata loaded with [PedSV::load.sample.metadata]
#' @param X Vector of values for primary independent variable. See `Details`.
#' @param Y Vector of values for dependent variable. See `Details`.
#' @param use.N.pcs Specify how many principal components should be adjusted in
#' model \[default: 3\]
#' @param extra.terms Specify if any extra terms should be added to the model.
#' Named options include:  "batch", "coverage, "insert.size", and "wgd".
#' Custom terms can be passed using their exact column names in `meta`.
#'
#' @details There are several options for providing `X` and `Y` values:
#'  * As an unnamed vector. In this case, the values are assumed to be in the
#'  same order as the samples in `meta`.
#'  * As a named vector. In this case, the vector names are assumed to be sample
#'  IDs, and any sample ID failing to match in `meta` will be dropped.
#'  * As a data frame. In this case, the row names are assumed to be sample IDs,
#'  and any sample ID failing to match in `meta` will be dropped.
#'
#' @returns data.frame
#'
#' @seealso [pedsv.glm]
#'
#' @export prep.glm.matrix
#' @export
prep.glm.matrix <- function(meta, X, Y, use.N.pcs=3, extra.terms=NULL){
  # Standard covariates
  df <- data.frame("is.female" = as.numeric(meta$inferred_sex == "FEMALE"
                                            | meta$chrX_CopyNumber > 1.5),
                   "trio.phase" = as.numeric(meta$study_phase == "trio"),
                   row.names=rownames(meta))
  if(use.N.pcs > 0){
    df <- cbind(df, apply(meta[paste("PC", 1:use.N.pcs, sep="")], 2, scale))
  }
  if(!is.null(extra.terms) & length(extra.terms) > 0){
    if("batch" %in% extra.terms){
      df$batch = as.factor(meta$batch)
    }
    if("coverage" %in% extra.terms){
      df$coverage = scale(as.numeric(meta$median_coverage))
    }
    if("insert.size" %in% extra.terms){
      df$insert.size <- scale(as.numeric(meta$insert_size))
    }
    if("wgd" %in% extra.terms){
      df$abs.wgd = scale(abs(as.numeric(meta$wgd_score)))
    }
    other.terms <- setdiff(extra.terms, c("cohort", "batch", "coverage",
                                          "insert.size", "wgd"))
    for(term in other.terms){
      if(term %in% colnames(meta)){
        df[, term] <- scale(meta[, term])
      }else{
        warning(message=paste("Other term '", term, "' not found in sample metadata. ",
                              "Will be ignored.", sep=""))
      }
    }
  }

  # Add X and Y values
  for(col.to.add in list(list("X", X), list("Y", Y))){
    new.col.name <- as.character(col.to.add[[1]])
    vals <- col.to.add[[2]]
    if(is.null(names(vals))){
      df[, new.col.name] <- as.numeric(vals)
    }else{
      if(!is.data.frame(vals)){
        vals <- as.data.frame(vals, row.names=names(vals))
      }
      colnames(vals) <- new.col.name
      df <- merge(df, vals, by=0, sort=F, all=F)
      rownames(df) <- df$Row.names; df$Row.names <- NULL
    }
  }

  # Drop redundant features
  drop.idx <- which(apply(df, 2, function(v){length(unique(v)) == 1}))
  if(length(drop.idx) > 0){
    df <- df[, -drop.idx]
  }

  # Standard normalize all covariates
  skip.standardize.variables <- c("X", "Y", "batch", "is.female", "batch", "trio.phase")
  standardize.variables <- setdiff(colnames(df), skip.standardize.variables)
  if(length(standardize.variables) > 0){
    df[, standardize.variables] <- apply(df[, standardize.variables], 2, scale)
  }

  # Ensure X and Y have at least two distinct values
  if(!all(c("X", "Y") %in% colnames(df))){
    warning("prep.glm.matrix: must supply non-separable X and Y values")
  }

  return(df[complete.cases(df), ])
}


#' Get eligible samples
#'
#' Determine the sample IDs eligible for a specified case:control comparison
#'
#' @param meta Sample metadata as loaded by [load.sample.metadata]
#' @param cancer Cancer code; one of "pancan", "OS", "NBL", or "EWS"
#'
#' @return List of character vectors of sample IDs for cases and controls
#'
#' @details If `meta` contains columns for cancer-specific case and control designations,
#' those will be used; otherwise all eligible controls will be returned
#'
#' @export get.eligible.samples
#' @export
get.eligible.samples <- function(meta, cancer){
  control.idx <- which(meta$disease == "control"
                       & (!meta$proband & !is.na(meta$family_id)
                          | is.na(meta$proband)))
  if(cancer == "pancan"){
    case.idx <- which(meta$disease != "control"
                      & (meta$proband | is.na(meta$proband)))
  }else{
    case.idx <- which(metadata.cancer.label.map[meta$disease] == cancer
                      & (meta$proband | is.na(meta$proband)))
  }
  case.cname <- paste(cancer, "case", sep="_")
  if(case.cname %in% colnames(meta)){
    case.idx <- intersect(case.idx, which(meta[, case.cname]))
  }
  control.cname <- paste(cancer, "control", sep="_")
  if(control.cname %in% colnames(meta)){
    control.idx <- intersect(control.idx, which(meta[, control.cname]))
  }else{
    control.idx <- which(meta$disease == "control"
                         & (meta$proband | is.na(meta$proband)))
  }
  list("cases" = rownames(meta)[case.idx],
       "controls" = rownames(meta)[control.idx])
}


#' Build a phenotype indicator vector
#'
#' Construct a one-hot encoding for disease status for a set of cases and controls
#'
#' @param case.ids Vector of case IDs
#' @param control.ids Vector of control IDs
#'
#' @returns Numeric vector
#'
#' @export get.phenotype.vector
#' @export
get.phenotype.vector <- function(case.ids, control.ids){
  Y <- as.numeric(c(case.ids, control.ids) %in% case.ids)
  names(Y) <- c(case.ids, control.ids)
  return(Y)
}


#' PedSV generic association test
#'
#' Fit a GLM for user-defined dependent & independent variables while adjusting
#' for standard sample covariates
#'
#'
#' @param meta Sample metadata loaded with [PedSV::load.sample.metadata]
#' @param X Vector of values for primary independent variable. See [PedSV::prep.glm.matrix].
#' @param Y Vector of values for dependent variable. See [PedSV::prep.glm.matrix].
#' @param use.N.pcs Specify how many principal components should be adjusted in
#' model \[default: 3\]
#' @param family `family` parameter passed to [glm]
#' @param extra.terms Extra covariate terms to include in model. See [PedSV::prep.glm.matrix].
#' @param firth.fallback Attempt to use Firth bias-reduced logistic regression when
#' traditional logistic regression fails to converge or dataset is quasi-separable
#' \[default: TRUE\]
#' @param strict.fallback Implement Firth regression if a standard logit model returns
#' any errors or warnings. Setting this to `FALSE` will only default to Firth
#' regression if logit returns any errors or if the standard error of the genotype
#' coefficient exceeds `se.tolerance` \[default: TRUE\]
#' @param se.tolerance If `firth.fallback` is `TRUE` but `strict.fallback` is
#' `FALSE`, only use Firth regression if a standard logit model produces a
#' genotype effect standard error exceeding this value \[default: 10\]
#' @param firth.always Always use Firth regression \[default: FALSE\]
#' @param return.fit.summary Return the full summary of the fitted model. Only recommended
#' for debugging purposes. Not recommended for analysis. \[default: FALSE\]
#'
#' @return Named vector of test statsitics corresponding to independent variable
#'
#' @seealso [PedSV::prep.glm.matrix]
#'
#' @export pedsv.glm
#' @export
pedsv.glm <- function(meta, X, Y, use.N.pcs=3, family=gaussian(), extra.terms=NULL,
                      firth.fallback=TRUE, strict.fallback=TRUE,
                      se.tolerance=10, firth.always=FALSE,
                      return.all.coefficients=FALSE){
  # Ensure Firth package is loaded
  require(logistf, quietly=TRUE)

  # Build dataframe of covariates
  test.df <- prep.glm.matrix(meta, X, Y, use.N.pcs, extra.terms)

  # Check to make sure dataset is not separable
  if(!all(c("X", "Y") %in% colnames(test.df))){
    return(c(NA, NA, NA, NA, NA))
  }
  ct <- table(test.df[, c("X", "Y")])
  if(nrow(ct) <= 1 | ncol(ct) <= 1){
    return(c(NA, NA, NA, NA, NA))
  }

  # Fit GLM
  logit.regression <- function(data){
    glm(Y ~ ., data=data, family=family)
  }
  firth.regression <- function(data){
    tryCatch(logistf(Y ~ ., data=data, control=logistf.control(maxit=100, maxstep=-1), flic=TRUE),
             error=function(e){
               tryCatch(logit.regression(data),
                        error=function(e){
                          c(NA, NA, NA, NA, "flic")
                        })
             })
  }
  if(firth.always & family$family == "binomial"){
    fit <- tryCatch(firth.regression(test.df),
                    error=function(e){logit.regression(test.df)})
  }else{
    if(firth.fallback & family$family == "binomial"){
      if(strict.fallback){
        fit <- tryCatch(logit.regression(test.df),
                        warning=function(w){firth.regression(test.df)},
                        error=function(e){firth.regression(test.df)})
      }else{
        fit <- tryCatch(logit.regression(test.df),
                        error=function(e){firth.regression(test.df)})
      }
      if(fit$method == "glm.fit"){
        if(summary(fit)$coefficients["X", 2] > se.tolerance){
          fit <- firth.regression(test.df)
        }
      }
    }else{
      fit <- logit.regression(test.df)
    }
  }
  if(length(fit) == 5 & is.na(fit[1])){
    return(fit)
  }else{
    firth <- if(fit$method != "glm.fit"){TRUE}else{FALSE}
  }

  # Extract coefficient corresponding to independent variable
  # Point estimate, stderr, Z-score, P-value
  if(!return.all.coefficients){
    if(firth){
      beta <- fit$coefficients["X"]
      c(as.numeric(c(beta,
                     sqrt(diag(vcov(fit)))["X"],
                     (if(beta >= 0){1}else{-1}) * sqrt(qchisq(1-fit$prob["X"], df=1)),
                     fit$prob["X"])), "flic")
    }else{
      c(summary(fit)[["coefficients"]]["X", ], "glm")
    }
  }else{
    summary(fit)
  }
}


#' Survival model of SV size
#'
#' Constructs a Cox proportional hazards model of the largest SV per sample
#' versus case-control status
#'
#' @param sv.bed SV bed that was imported with [load.sv.bed]
#' @param meta Sample metadata that was imported with [load.sample.metadata]
#' @param ad Either a path to an allele dosage .bed.gz matrix or a data.frame
#' from a previous call of [query.ad.matrix]
#' @param cancer Cancer type of interest \[default: "pancan"\]
#' @param use.N.pcs Specify how many principal components should be adjusted in
#' model \[default: 3\]
#' @param extra.terms Specify if any extra terms should be added to the model.
#' Named options include:  "batch", "coverage, "insert.size", and "wgd".
#' Custom terms can be passed using their exact column names in `meta`.
#'
#' @export sv.size.survfit
#' @export
sv.size.survfit <- function(sv.bed, meta, ad, cancer="pancan", use.N.pcs=3,
                            extra.terms=NULL){
  require(survival, quietly=TRUE)
  # Subset data to cases & controls meeting inclusion criteria
  sids <- get.eligible.samples(meta, cancer)
  case.ids <- sids$cases
  control.ids <- sids$controls
  all.sids <- unlist(sids)
  Y <- get.phenotype.vector(case.ids, control.ids)

  # Get verbose AD matrix for all SVs in sv.bed and convert to binary 0|1 indicator
  ad.df <- query.ad.from.sv.bed(ad, sv.bed, action="verbose",
                                keep.samples=all.sids)
  ad.df <- as.data.frame(apply(ad.df, 2, function(vals){
    vals[which(vals > 1)] <- 1; return(vals)
  }))

  # Assign X value as largest log10(SVLEN) per sample
  sv.bed$log10.SVLEN <- sapply(sv.bed$SVLEN, function(v){log10(max(c(v, 1)))})
  X <- sapply(colnames(ad.df), function(sid){
    max(ad.df[, sid] * sv.bed[rownames(ad.df), "log10.SVLEN"], na.rm=T)
  })

  # Make association df and customize to be appropriate for survival
  test.df <- prep.glm.matrix(meta, X, Y, use.N.pcs=use.N.pcs,
                             extra.terms=extra.terms)
  test.df$event <- 1

  # Fit survival models for cases and controls separately (for visualization)
  surv.fits <- list("control" = summary(survfit(with(test.df[which(test.df$Y == 0), ],
                                                     Surv(X, event)) ~ 1)),
                    "case" = summary(survfit(with(test.df[which(test.df$Y == 1), ],
                                                  Surv(X, event)) ~ 1)))

  # Fit survival model of case|control ~ largest SV length
  cox.res <- summary(coxph(Surv(X, event) ~ ., data=test.df))

  return(list("surv.models" = surv.fits,
              "cox.stats" = cox.res))
}

