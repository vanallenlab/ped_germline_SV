#!/usr/bin/env R

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>


# Plotting helper functions


#' Add phenotype barplots
#'
#' Add stacked barplots by phenotype to an existing plot
#'
#' @param plot.df Plotting dataframe. See [barplot.by.phenotype()].
#' @param bar.mids Midpoints for all bars (controls first)
#' @param bar.hex Bar width \[default: 0.5\]
#' @param control.sep Separation constant between cases & controls \[default: 0.1875\]
#' @param horiz Should the bars be plotted horizontally? \[default: TRUE\]
#' @param color.by.sig Should the case bars be shaded by significance? \[default: TRUE\]
#' @param legend.on.bars Should a legend be written on the bars? \[default: FALSE\]
#' @param case.label Label for case bars if `legend.on.bars` is `TRUE` \[default: "Case"\]
#' @param control.label Label for control bars if `legend.on.bars` is `TRUE` \[default: "Control"\]
#' @param add.pvals Should P values be annotated on the opposite margin?
#' \[default: TRUE\]
#' @param cancer.types.override Specify order of cancer types in `plot.df`
#' \[default: infer from `rownames(plot.df)` \]
#'
#' @seealso [barplot.by.phenotype()]
#'
#' @export add.pheno.bars
#' @export
add.pheno.bars <- function(plot.df, bar.mids, bar.hex=0.5, control.sep=0.1875,
                           horiz=TRUE, color.by.sig=TRUE, legend.on.bars=FALSE,
                           case.label="Case", control.label="Control",
                           add.pvals=FALSE, cancer.types.override=NULL){
  # Infer positional parameters
  n.pheno <- nrow(plot.df)
  # if(!horiz){plot.df <- plot.df[n.pheno:1, ]}
  control.idx <- 1:n.pheno
  case.idx <- control.idx + n.pheno
  staggered.bar.mids <- c(bar.mids+control.sep, bar.mids-control.sep)
  bar.width.min <- staggered.bar.mids - (bar.hex/2)
  bar.width.max <- staggered.bar.mids + (bar.hex/2)
  bar.val.min <- rep(0, 2*n.pheno)
  bar.val.max <- as.numeric(c(plot.df[, 4], plot.df[, 1]))
  ci.min <- as.numeric(c(plot.df[, 5], plot.df[, 2]))
  ci.max <- as.numeric(c(plot.df[, 6], plot.df[, 3]))
  longest.control <- head(which(plot.df[, 4] == max(plot.df[, 4], na.rm=T)), 1)
  longest.case <- head(which(plot.df[, 1] == max(plot.df[, 1], na.rm=T)), 1) + n.pheno

  # Set X and Y values according to value of horiz
  if(horiz){
    bar.x.left <- bar.val.min
    bar.x.right <- bar.val.max
    bar.y.bottom <- bar.width.min
    bar.y.top <- bar.width.max
    ci.x0 <- ci.min
    ci.x1 <- ci.max
    ci.y0 <- staggered.bar.mids
    ci.y1 <- staggered.bar.mids
    legend.text.x <- rep(0-(0.05*diff(par("usr")[1:2])), 2)
    legend.text.y <- c(staggered.bar.mids[longest.control]+(bar.hex/5),
                       staggered.bar.mids[longest.case]+(bar.hex/10))
    pval.axis <- 4
  }else{
    bar.x.left <- bar.width.min
    bar.x.right <- bar.width.max
    bar.y.bottom <- bar.val.min
    bar.y.top <- bar.val.max
    ci.x0 <- staggered.bar.mids
    ci.x1 <- staggered.bar.mids
    ci.y0 <- ci.min
    ci.y1 <- ci.max
    legend.text.x <- c(staggered.bar.mids[longest.control]+(bar.hex/5),
                       staggered.bar.mids[longest.case]+(bar.hex/10))
    legend.text.y <- rep(0-(0.05*diff(par("usr")[3:4])), 2)
    pval.axis <- 3
  }

  # Set coloring for case bars based on cancer type and value of color.by.sig
  sig.idx <- which(plot.df[, 7] < 0.05)
  ctypes <- if(!is.null(cancer.types.override)){cancer.types.override}else{rownames(plot.df)}
  bar.pals <- lapply(ctypes, function(pheno){
    if(pheno %in% names(cancer.palettes)){
      cancer.palettes[[pheno]]
    }else{
      cancer.palettes[["pancan"]]
    }
  })
  bar.cols <- sapply(bar.pals, function(pal){pal["light1"]})
  ci.cols <- sapply(bar.pals, function(pal){pal["main"]})
  if(color.by.sig){
    for(i in sig.idx){
      bar.cols[i] <- bar.pals[[i]]["main"]
      ci.cols[i] <- bar.pals[[i]]["dark1"]
    }
  }

  # Add control bars
  rect(xleft=bar.x.left[control.idx], xright=bar.x.right[control.idx],
       ybottom=bar.y.bottom[control.idx], ytop=bar.y.top[control.idx],
       col=control.colors[["main"]], border=NA, bty="n")
  segments(x0=ci.x0[control.idx], x1=ci.x1[control.idx],
           y0=ci.y0[control.idx], y1=ci.y1[control.idx],
           lwd=2, lend="butt", col=cancer.palettes[["control"]]["dark1"])
  if(legend.on.bars){
    text(x=legend.text.x[1], y=legend.text.y[1],
         pos=4, label=control.label, cex=4/6, col=control.colors[["dark2"]])
  }
  rect(xleft=bar.x.left[control.idx], xright=bar.x.right[control.idx],
       ybottom=bar.y.bottom[control.idx], ytop=bar.y.top[control.idx],
       col=NA, xpd=T)

  # Add case bars
  rect(xleft=bar.x.left[case.idx], xright=bar.x.right[case.idx],
       ybottom=bar.y.bottom[case.idx], ytop=bar.y.top[case.idx],
       col=bar.cols, border=NA, bty="n")
  segments(x0=ci.x0[case.idx], x1=ci.x1[case.idx],
           y0=ci.y0[case.idx], y1=ci.y1[case.idx],
           lwd=2, lend="butt", col=ci.cols)
  if(legend.on.bars){
    text(x=legend.text.x[2], y=legend.text.y[2], pos=4, label=case.label,
         cex=4/6, col=bar.pals[[longest.case-n.pheno]][["dark2"]])
  }
  rect(xleft=bar.x.left[case.idx], xright=bar.x.right[case.idx],
       ybottom=bar.y.bottom[case.idx], ytop=bar.y.top[case.idx],
       col=NA, xpd=T)

  # Add P values, if optioned
  if(add.pvals){
    sapply(1:nrow(plot.df), function(i){
      pval <- as.numeric(plot.df[i, 7])
      if(color.by.sig){
        p.col <- if(!is.na(pval) & pval < 0.05){"black"}else{control.colors[["main"]]}
      }else{
        p.col <- "black"
      }
      p.label <- if(is.na(pval)){"NA"}else{PedSV::format.pval(pval, nsmall=0)}
      axis(pval.axis, at=bar.mids[i], tick=F, line=-0.9, las=2, cex.axis=5/6,
           labels=p.label, col.axis=p.col)
    })
  }
}


#' Convert landscape test stats to barplot-compliant data.frame
#'
#' Transform the association summary statistics from a single landscape
#' case:control burden test into a format expected by [barplot.by.phenotype()]
#'
#' @param stats data.frame with test statistics as output by [pedsv.glm()]
#' @param ci.mode mode for confidence interval calculation \[choices: `binomial`, `normal`; default: `normal`\]
#' @param order.by.cancer order the results by study-wide convention \[default: TRUE\]
#'
#' @seealso [barplot.by.phenotype()], [pedsv.glm()]
#'
#' @return data.frame
#'
#' @export stats2barplotdf
#' @export
stats2barplotdf <- function(stats, ci.mode="normal", order.by.cancer=TRUE){
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
                        "p" = pvals)
  if(order.by.cancer){
    rownames(plot.df) <- stats$disease
    plot.df <- plot.df[intersect(names(cancer.colors[1:4]), rownames(plot.df)), ]
  }
  return(plot.df)
}
