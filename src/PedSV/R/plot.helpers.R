#!/usr/bin/env R

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, and the Van Allen Laboratory
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
                       max.ticks=6, title=NULL, tck=-0.025, cex.axis=5/6,
                       label.line=-0.65, title.line=0.5,
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
    if(is.numeric(labels[i])){
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
    base <- gsub(" ", "", formatC(round(as.numeric(parts[1]), nsmall), digits=1+nsmall), fixed=T)
    exp <- gsub(" ", "", as.numeric(parts[2]), fixed=T)
    bquote(italic(P) ~ .(equality) ~ .(base) ~ "x" ~ 10 ^ .(exp))
  }else{
    bquote(italic(P) ~ .(equality) ~ .(formatC(round(p, max.decimal), digits=max.decimal)))
  }
}
