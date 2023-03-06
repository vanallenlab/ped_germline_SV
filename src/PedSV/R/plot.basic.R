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
#' @param legend.vals Named vector mapping category names to colors \[default: NULL\]
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
pc.scatterplot <- function(pcs, pc.X, pc.Y, colors, title=NULL, legend.vals=NULL,
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
  clean.axis(1, title=paste("Principal Component", pc.X), infinite=T)
  clean.axis(2, title=paste("Principal Component", pc.Y), infinite=T)
  points(x, y, pch=19, cex=0.3, col=colors, xpd=T)
  mtext(3, text=title, font=2)

  # Add legend, if optioned
  if(!is.null(legend.vals)){
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
