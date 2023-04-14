#!/usr/bin/env R

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>


# Basic plotting functions to generate entire plots
# See plot.helpers.R for smaller subroutines used to generate certain plot elements


#' Plot principal components
#'
#' Generate a scatterplot of principal components with options for coloring
#'
#' @param pcs Matrix of principal components as loaded by [PedSV::load.pc.matrix]
#' @param pc.X PC to be plotted on the X-axis (see `Details`)
#' @param pc.Y PC to be plotted on the Y-axis (see `Details`)
#' @param colors Vector of colors. Must be supplied in the same order as rows in `pcs`
#' @param title Main title for plot \[default: NULL\]
#' @param x.title Title for X-axis
#' @param x.title.line Line for X-axis title (`title.line` parameter for [PedSV::clean.axis])
#' @param y.title Title for Y-axis
#' @param y.title.line Line for axis titles (`title.line` parameter for [PedSV::clean.axis])
#' @param legend.vals Named vector mapping category names to colors \[default: NULL\]
#' @param legend.labels Optional vector to overwrite names of `legend.vals`
#' @param cex Character expansion factor for individual points \[default: 0.3\]
#' @param parmar Numeric vector of values to pass to `par(mar)`
#'
#' @details `pc.X` and `pc.Y` can be specified as numeric column indexes (e.g., `1`),
#' or can be specified by column name (e.g., `PC1`).
#'
#' If `legend.vals` is not provided (or is `NULL`), no legend will be added.
#'
#' @export pc.scatterplot
#' @export
pc.scatterplot <- function(pcs, pc.X, pc.Y, colors, title=NULL, x.title=NULL,
                           x.title.line=0.5, y.title=NULL, y.title.line=0.5,
                           legend.vals=NULL, legend.labels=NULL,
                           cex=0.3, parmar=c(2.5, 2.5, 1, 1)){
  if(!is.numeric(pc.X)){
    pc.X <- which(colnames(pcs) == pc.X)
  }
  if(!is.numeric(pc.Y)){
    pc.Y <- which(colnames(pcs) == pc.Y)
  }
  x <- pcs[, pc.X]
  y <- pcs[, pc.Y]

  prep.plot.area(xlims=range(x), ylims=range(y), parmar=parmar, xaxs="r", yaxs="r")
  mtext(3, text=title, font=2)

  if(is.null(x.title)){
    x.title <- paste("Principal Component", pc.X)
  }
  clean.axis(1, title=x.title, infinite=T, title.line=x.title.line)

  if(is.null(y.title)){
    y.title <- paste("Principal Component", pc.Y)
  }
  clean.axis(2, title=y.title, infinite=T, title.line=y.title.line)

  points(x, y, pch=19, cex=0.3, col=colors, xpd=T)

  # Add legend, if optioned
  if(!is.null(legend.vals)){
    if(!is.null(legend.labels)){
      names(legend.vals) <- legend.labels
    }
    quad.counts <- table(x > mean(par("usr")[1:2]), y < mean(par("usr")[3:4]))
    least.dense.quad <- head(which(quad.counts == min(quad.counts)), 1)
    legend.pos <- if(least.dense.quad == 1){
      "topleft"
    }else if(least.dense.quad == 2){
      "topright"
    }else if(least.dense.quad == 3){
      "bottomleft"
    }else if(least.dense.quad == 4){
      "bottomright"
    }
    legend(legend.pos, legend=names(legend.vals), pch=21, pt.bg=legend.vals,
           pt.cex=1.5, bty="n", border=NA)
  }
}


#' Cowplot
#'
#' Generate a cowplot using base R syntax
#'
#' @param data list of [density] objects to plot
#' @param names optional list of names for Y axis \[default: take names from data\]
#' @param cow.overlap relative fraction of overlap bewtween adjacent cows \[default: 0.35\]
#' @param xlims custom X axis limits
#' @param ylims custom Y axis limits
#' @param fill vector of polygon fill colors \[default: "grey70"\]
#' @param border vector of cow border colors \[default: "grey35"\]
#' @param border.lwd line width for cow borders \[default: 2\]
#' @param parmar vector of values passed to par(mar)
#'
#' @export cowplot
#' @export
cowplot <- function(data, names=NULL, cow.overlap=0.35, xlims=NULL, x.axis=TRUE,
                    fill="grey70", border="grey35", border.lwd=2,
                    parmar=c(2.5, 3, 0.25, 0.25)){
  # Get names before manipulating data
  if(is.null(names)){
    names <- names(data)
    if(is.null(names)){
      names <- 1:length(data)
    }
  }

  # Scale Y values of data to [0, cow.overlap]
  for(i in 1:length(data)){
    y <- data[[i]]$y
    data[[i]]$y <- (1 + cow.overlap) * (y / max(y))
  }

  # Get plot dimensions
  if(is.null(xlims)){
    xlims <- c(min(sapply(data, function(d){min(d$x)})),
               max(sapply(data, function(d){max(d$x)})))
  }
  ylims <- c(0, length(data) + cow.overlap)

  # Prep plot area
  prep.plot.area(xlims, ylims, parmar, xaxs="i", yaxs="r")
  if(x.axis){
    clean.axis(1, title="Values", infinite=TRUE)
  }
  axis(2, at=(1:length(data)) - 0.5, tick=F, las=2, line=-0.8, labels=names)

  # Add cows
  sapply(length(data):1, function(i){
    x <- c(data[[i]]$x, rev(data[[i]]$x))
    y <- c(data[[i]]$y, rep(0, times=length(data[[i]]$y))) + (i-1)
    polygon(x, y, border=border, col=fill, lwd=2)
  })
}

