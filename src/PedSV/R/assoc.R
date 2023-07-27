#!/usr/bin/env R

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, and the Van Allen Laboratory
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
#' Named options include: "study", "cohort", "coverage, "insert.size",
#' and "wgd". Custom terms can be passed using their exact column names in `meta`.
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
                   row.names=rownames(meta))
  if(use.N.pcs > 0){
    df <- cbind(df, apply(meta[paste("PC", 1:use.N.pcs, sep="")], 2, scale))
  }
  if(!is.null(extra.terms)){
    if("cohort" %in% extra.terms){
      df$trio.phase = as.numeric(meta$study_phase == "trio")
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
    other.terms <- setdiff(extra.terms, c("cohort", "coverage", "insert.size", "wgd"))
    for(term in other.terms){
      df[, term] <- scale(meta[, term])
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
  df[, setdiff(colnames(df), c("X", "Y"))] <- apply(df[, setdiff(colnames(df), c("X", "Y"))], 2, scale)

  # Ensure X and Y have at least two distinct values
  if(!all(c("X", "Y") %in% colnames(df))){
    stop("Error in prep.glm.matrix: must supply non-separable X and Y values")
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
#' @details If `meta` contains columns for cancer-specific control designation,
#' those will be used; otherwise all eligible controls will be returned
#'
#' @export get.eligible.samples
#' @export
get.eligible.samples <- function(meta, cancer){
  control.idx <- which(meta$disease == "control"
                       & (!meta$proband & !is.na(meta$family_id))
                       | is.na(meta$proband))
  if(cancer == "pancan"){
    case.idx <- which(meta$disease != "control"
                      & (meta$proband | is.na(meta$proband)))
  }else{
    case.idx <- which(metadata.cancer.label.map[meta$disease] == cancer
                      & (meta$proband | is.na(meta$proband)))
  }
  control.cname <- paste(cancer, "control", sep="_")
  if(control.cname %in% colnames(meta)){
    control.idx <- which(meta[, control.cname])
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
#' model \[default: 10\]
#' @param family `family` parameter passed to [glm]
#' @param extra.terms Extra covariate terms to include in model. See [PedSV::prep.glm.matrix].
#' @param firth.fallback Attempt to use Firth bias-reduced logistic regression when
#' traditional logistic regression fails to converge or dataset is quasi-separable
#' \[default: TRUE\]
#' @param strict.fallback Implement Firth regression if a standard logit model returns
#' any errors or warnings. Setting this to `FALSE` will only default to Firth
#' regression if logit returns any errors or if the standard error of the genotype
#' coefficient exceeds `nonstrict.se.tolerance` \[default: TRUE\]
#' @param nonstrict.se.tolerance If `strict.fallback` is `FALSE`, only use Firth
#' regression if a standard logit model produces a genotype effect standard error
#' exceeding this value \[default: 10\]
#' @param firth.always Always use Firth regression \[default: FALSE\]
#'
#' @return Named vector of test statsitics corresponding to independent variable
#'
#' @seealso [PedSV::prep.glm.matrix]
#'
#' @export pedsv.glm
#' @export
pedsv.glm <- function(meta, X, Y, use.N.pcs=10, family=gaussian(), extra.terms=NULL,
                      firth.fallback=TRUE, strict.fallback=TRUE,
                      nonstrict.se.tolerance=10, firth.always=FALSE){
  # Ensure Firth package is loaded
  require(logistf, quietly=TRUE)

  # Build dataframe of covariates
  test.df <- prep.glm.matrix(meta, X, Y, use.N.pcs, extra.terms)

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
        if(fit$method == "glm.fit"){
          if(summary(fit)$coefficients["X", 2] > nonstrict.se.tolerance){
            fit <- firth.regression(test.df)
          }
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
  # Point estimate, stderr, test statistic, P-value
  if(firth){
    c(as.numeric(c(fit$coefficients["X"],
                   sqrt(diag(vcov(fit)))["X"],
                   qchisq(1-fit$prob, df=1)["X"],
                   fit$prob["X"])), "flic")
  }else{
    c(summary(fit)[["coefficients"]]["X", ], "glm")
  }
}
