#!/usr/bin/env R

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>


# Plotting helper functions


#' Prepare Plot Area
#'
#' Prepare a standardized & formatted plot area
#'
#' @param xlims Range of values for X axis
#' @param ylims Range of values for Y axis
#' @param parmar Margin values passed to par()
#' @param xaxs Value of `xaxs` passed to plot()
#' @param yaxs Value of `yaxs` passed to plot()
#'
#' @examples
#' prep.plot.area(xlims=c(0, 5), ylims=(-10, 10), parmar=rep(3, 4));
#'
#' @export prep.plot.area
#' @export
prep.plot.area <- function(xlims, ylims, parmar, xaxs="i", yaxs="i"){
  par(mar=parmar, bty="n")
  plot(NA, xlim=xlims, ylim=ylims, type="n",
       xaxs=xaxs, xlab="", xaxt="n",
       yaxs=yaxs, ylab="", yaxt="n")
}


#' Clean axis
#'
#' Print a clean axis using visually pleasing defaults
#'
#' @param side Value passed to `axis()`. See `?axis` for details.
#' @param at Positions where axis ticks should be plotted \[default: `axTicks(side)`\]
#' @param labels Labels for axis ticks \[default: values of `at`\]
#' @param labels.at Positions for axis labels \[default: values of `at`\]
#' @param label.units Specify custom units for the axis label. Default of `NULL`
#' will display numeric values. Options currently include "percent" for percentages.
#' Can be overridden by supplying `labels` directly.
#' @param parse.labels Should `labels` be parsed as R expressions? \[default: FALSE\]
#' @param max.ticks Maximum number of axis ticks. Will be overridden by `at` \[default: 6\]
#' @param title Axis title
#' @param tck Value passed to `axis()`. See `?axis` for details. \[default: -0.025\]
#' @param cex.axis Value passed to `axis()`. See `?axis` for details. \[default: 5/6\]
#' @param label.line `line` parameter for axis labels \[default: -0.65\]
#' @param title.line `line` parameter for axis title \[default: 0.5\]
#' @param infinite Indicator for the axis to be extended infinitely (without ticks) \[default: FALSE\]
#' @param infinite.positive Indicator for the axis to be extended infinitely
#' in the positive direction (without ticks) \[default: FALSE\]
#' @param infinite.negative Indicator for the axis to be extended infinitely
#' in the negative direction (without ticks) \[default: FALSE\]
#'
#' @returns NULL
#'
#' @export clean.axis
#' @export
clean.axis <- function(side, at=NULL, labels=NULL, labels.at=NULL, label.units=NULL,
                       parse.labels=FALSE, max.ticks=6, title=NULL, tck=-0.025,
                       cex.axis=5/6, label.line=-0.65, title.line=0.5,
                       infinite=FALSE, infinite.positive=FALSE, infinite.negative=FALSE){
  if(infinite){axis(side, at=c(-10e10, 10e10), tck=0, labels=NA)}
  if(is.null(at)){
    at <- axTicks(side)
    if(length(at) > max.ticks){
      at <- at[c(TRUE, FALSE)]
    }
  }
  if(infinite.positive){axis(side, at=c(at[1], 10e10), tck=0, labels=NA)}
  if(infinite.negative){axis(side, at=c(-10e10, at[length(at)]), tck=0, labels=NA)}
  if(is.null(labels)){
    labels <- at
    if(!is.null(label.units)){
      if(label.units == "percent"){
        labels <- paste(100 * labels, "%", sep="")
      }
    }
  }
  if(is.null(labels.at)){labels.at <- at}
  if(side %in% c(1, 3)){
    las <- 1
    title.at <- mean(par("usr")[1:2])
  }else{
    las <- 2
    title.at <- mean(par("usr")[3:4])
  }
  axis(side, at=at, labels=NA, tck=tck)
  sapply(1:length(labels.at), function(i){
    if(parse.labels){
      label <- parse(text=labels[i])
    }else if(is.numeric(labels[i])){
      label <- prettyNum(labels[i], big.mark=",")
    }else{
      label <- labels[i]
    }
    axis(side, at=labels.at[i], labels=label, tick=F, cex.axis=cex.axis,
         las=las, line=label.line)
  })
  if(!is.null(title)){
    axis(side, at=title.at, tick=F, labels=title, line=title.line, xpd=T)
  }
}


#' Convert colors to greyscale
#'
#' Convert one or more HEX-encoded color(s) to its greyscale HEX equivalent(s)
#'
#' @param in.colors
#'
#' @return character vector of HEX color
#'
#' @examples
#' hex2grey(c("#FF0000", "#00008B", "#ffffe0"))
#'
#' @export hex2grey
#' @export
hex2grey <- function(in.colors){
  mean.rgbs <- round(apply(col2rgb(in.colors), 2, mean), 0)
  sapply(mean.rgbs, function(k){rgb(k, k, k, maxColorValue=255)})
}


#' Format P-value
#'
#' Format P-value for plotting
#'
#' @param p P-value
#' @param nsmall number of digits after the decimal to retain for scientific
#' notification \[default: 2\]
#' @param max.decimal convert all P-values requiring more digits after the decimal
#' to be converted to scientific notation \[default: 3\]
#' @param equality equality symbol to print after `P` \[default: '='\]
#' @param min.neg.log10.p minimum order of magnitude to process before considering
#' P-value to be arbitrarily/meaninglessly small \[default: 100\]
#'
#' @details Function borrowed from rCNV2 library (see Collins et al., Cell, 2022)
#'
#' @return formatted P-value as character
#'
#' @export format.pval
#' @export
format.pval <- function(p, nsmall=2, max.decimal=3, equality="=", min.neg.log10.p=100){
  if(-log10(p)>min.neg.log10.p){
    bquote(italic(P) %~~% 0)
  }else if(ceiling(-log10(p)) > max.decimal){
    parts <- unlist(strsplit(format(p, scientific=T), split="e"))
    base <- gsub(" ", "", format(round(as.numeric(parts[1]), nsmall), digits=1+nsmall), fixed=T)
    exp <- gsub(" ", "", as.numeric(parts[2]), fixed=T)
    if(base %in% c("1", "10")){
      if(base == "10"){
        exp <- as.character(as.numeric(exp) - 1)
      }
      bquote(italic(P) ~ .(equality) ~ 10 ^ .(exp))
    }else{
      bquote(italic(P) ~ .(equality) ~ .(base) ~ "x" ~ 10 ^ .(exp))
    }
  }else{
    bquote(italic(P) ~ .(equality) ~ .(formatC(round(p, max.decimal), digits=max.decimal)))
  }
}


#' Color points by density
#'
#' Generate colors for XY scatterplot based on point density
#'
#' @param x independent variable vector
#' @param y dependent variable vector
#' @param palette 256-color palette to be applied based on density \[default: viridis()\]
#' @param bandwidth `bandwidth` parameter passed to [densCols()]
#'
#' @details Inspired by heatscatter.R from Colby Chiang:
#'  https://github.com/cc2qe/voir/blob/master/bin/heatscatter.R
#'
#' @return dataframe of values to be plotted with density and colors
#'
#' @seealso [viridis()], [densCols()]
#'
#' @export
color.points.by.density <- function(x, y, palette=NULL, bandwidth=1){
  # Based on heatscatter.R from Colby Chiang
  # (https://github.com/cc2qe/voir/blob/master/bin/heatscatter.R)
  plot.df <- data.frame("x"=x, "y"=y)
  plot.df <- plot.df[which(!is.infinite(plot.df$x) & !is.infinite(plot.df$y)
                           & !is.na(plot.df$x) & !is.na(plot.df$y)), ]
  dens <- densCols(plot.df$x, plot.df$y, bandwidth=bandwidth,
                   colramp=colorRampPalette(c("black", "white")))
  plot.df$dens <- col2rgb(dens)[1, ] + 1L
  if(is.null(palette)){
    require(viridis, quietly=TRUE)
    palette <- viridis(256)
  }
  plot.df$col <- palette[plot.df$dens]
  plot.df[order(plot.df$dens), ]
}


#' Add phenotype barplots
#'
#' Add stacked barplots by phenotype to an existing plot
#'
#' @param plot.df Plotting dataframe. See [barplot.by.phenotype()].
#' @param bar.mids Midpoints for all bars (controls first)
#' @param bar.hex Bar width
#' @param control.sep Separation constant between cases & controls
#' @param horiz Should the bars be plotted horizontally? \[default: TRUE\]
#' @param color.by.sig Should the case bars be shaded by significance? \[default: TRUE\]
#' @param legend.on.bars Should a legend be written on the bars? \[default: FALSE\]
#'
#' @seealso [barplot.by.phenotype()]
#'
#' @export add.pheno.bars
#' @export
add.pheno.bars <- function(plot.df, bar.mids, bar.hex, control.sep, horiz=TRUE,
                           color.by.sig=TRUE, legend.on.bars=FALSE){
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
    legend.text.y <- c(staggered.bar.mids[longest.control]+(bar.hex/5)+control.sep,
                       staggered.bar.mids[longest.case]+(bar.hex/10)-control.sep)
  }else{
    bar.x.left <- bar.width.min
    bar.x.right <- bar.width.max
    bar.y.bottom <- bar.val.min
    bar.y.top <- bar.val.max
    ci.x0 <- staggered.bar.mids
    ci.x1 <- staggered.bar.mids
    ci.y0 <- ci.min
    ci.y1 <- ci.max
    legend.text.x <- c(staggered.bar.mids[longest.control]+(bar.hex/5)+control.sep,
                       staggered.bar.mids[longest.case]+(bar.hex/10)-control.sep)
    legend.text.y <- rep(0-(0.05*diff(par("usr")[3:4])), 2)
  }

  # Set coloring for case bars based on cancer type and value of color.by.sig
  sig.idx <- which(plot.df[, 7] < 0.05)
  bar.pals <- lapply(rownames(plot.df), function(pheno){
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
}




