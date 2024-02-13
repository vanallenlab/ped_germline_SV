#!/usr/bin/env R

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2024-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>


# Formal & secondary analysis helper functions


#' Convert landscape test stats to barplot-compliant data.frame
#'
#' Transform the association summary statistics from a single landscape
#' case:control burden test into a format expected by [barplot.by.phenotype()]
#'
#' @param stats data.frame with test statistics as output by [pedsv.glm()]
#' @param ci.mode mode for confidence interval calculation \[choices: `binomial`, `normal`; default: `normal`\]
#'
#' @seealso [barplot.by.phenotype()], [pedsv.glm()]
#'
#' @return data.frame
#'
#' @export stats2barplotdf
#' @export
stats2barplotdf <- function(stats, ci.mode="normal"){
  n.phenos <- nrow(stats)
  values <- as.numeric(c(stats$case.mean, stats$control.mean))
  ns <- as.numeric(c(stats$n.case, stats$n.control))
  if(ci.mode == "normal"){
    stdevs <- as.numeric(c(stats$case.stdev, stats$control.stdev))
    ci.margins <- stdevs / sqrt(ns) * qnorm(0.975)
    ci.lowers <- values - ci.margins
    ci.uppers <- values + ci.margins
  }else if(ci.mode == "binomial"){
    require(Hmisc, quietly=TRUE)
    cis <- Hmisc::binconf(x=round(values * ns), n=ns)
    ci.lowers <- cis[, "Lower"]
    ci.uppers <- cis[, "Upper"]
  }
  pvals <- as.numeric(stats$P.value)
  plot.df <- data.frame("case.value" = values[1:n.phenos],
                        "case.ci.lower" = ci.lowers[1:n.phenos],
                        "case.ci.upper" = ci.uppers[1:n.phenos],
                        "control.value" = values[(1:n.phenos) + n.phenos],
                        "control.ci.lower" = ci.lowers[(1:n.phenos) + n.phenos],
                        "control.ci.upper" = ci.uppers[(1:n.phenos) + n.phenos],
                        "p" = pvals, row.names=stats$disease)
  plot.df[intersect(names(cancer.colors[1:4]), rownames(plot.df)), ]
}

