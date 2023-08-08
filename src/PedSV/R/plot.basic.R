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

  points(x, y, pch=19, cex=cex, col=colors, xpd=T)

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


#' Ridgeplot
#'
#' Generate a ridgeplot using base R syntax
#'
#' @param data list of [density] objects to plot
#' @param names optional list of names for Y axis \[default: take names from data\]
#' @param hill.overlap relative fraction of overlap bewtween adjacent hills \[default: 0.35\]
#' @param xlims custom X axis limits
#' @param ylims custom Y axis limits
#' @param fill vector of polygon fill colors \[default: "grey70"\]
#' @param border vector of hill border colors \[default: "grey35"\]
#' @param border.lwd line width for hill borders \[default: 2\]
#' @param parmar vector of values passed to par(mar)
#'
#' @export ridgeplot
#' @export
ridgeplot <- function(data, names=NULL, hill.overlap=0.35, xlims=NULL, x.axis=TRUE,
                    fill="grey70", border="grey35", border.lwd=2,
                    parmar=c(2.5, 3, 0.25, 0.25)){
  # Get names before manipulating data
  if(is.null(names)){
    names <- names(data)
    if(is.null(names)){
      names <- 1:length(data)
    }
  }

  # Scale Y values of data to [0, hill.overlap]
  for(i in 1:length(data)){
    y <- data[[i]]$y
    data[[i]]$y <- (1 + hill.overlap) * (y / max(y))
  }

  # Get plot dimensions
  if(is.null(xlims)){
    xlims <- c(min(sapply(data, function(d){min(d$x)})),
               max(sapply(data, function(d){max(d$x)})))
  }
  ylims <- c(0, length(data) + hill.overlap)

  # Prep plot area
  prep.plot.area(xlims, ylims, parmar, xaxs="i", yaxs="r")
  if(x.axis){
    clean.axis(1, title="Values", infinite=TRUE)
  }
  axis(2, at=(1:length(data)) - 0.5, tick=F, las=2, line=-0.8, labels=names)

  # Add hills
  sapply(length(data):1, function(i){
    abline(h=i-1, col="gray85")
    x <- c(data[[i]]$x, rev(data[[i]]$x))
    y <- c(data[[i]]$y, rep(0, times=length(data[[i]]$y)))+i-1
    polygon(x, y, border=fill, col=fill, lwd=border.lwd)
    points(data[[i]]$x, data[[i]]$y+i-1, type="l", lwd=border.lwd, col=border)
  })
}


#' Barplot of values per cancer type
#'
#' Generate a barplot of one value for each cancer type
#'
#' @param plot.df Data.frame of values to plot. See `Details`.
#' @param bar.hex Relative width for each bar \[default: 0.5\]
#' @param case.control.sep Relatve separation between overlapping case and control bars \[default: 0.33\]
#' @param color.by.sig Should bars be shaded differently by significance level? \[default: TRUE\]
#' @param add.top.axis Should the top axis be added? \[default: TRUE\]
#' @param top.axis.units Specify units for top axis. Options are NULL (for
#' numeric) or "percent" \[default: NULL\]
#' @param title Custom title for top axis
#' @param parmar Value of `mar` passed to `par()`
#'
#' @details `plot.df` is expected to adhere to the following specifications:
#' * One row per cancer type
#' * Row names indicate cancer type (or "control")
#' * Columns ordered as follows: (1) case value to plot, (2) case lower CI bound,
#' (3) case upper CI bound, (4) control value to plot, (5) control lower CI bound,
#' (6) control upper CI bound, (7) P-value
#'
#' @export barplot.by.phenotype
#' @export
barplot.by.phenotype <- function(plot.df, bar.hex=0.5, case.control.sep=1/3,
                                 color.by.sig=TRUE, add.top.axis=TRUE, top.axis.units=NULL,
                                 title="Value", orient.cases="top",
                                 parmar=c(0.2, 4.1, 2.1, 4)){
  # Get plot dimensions
  xlims <- c(0, min(c(2*max(plot.df[, c(1, 4)]), max(plot.df[, 1:6]))))
  ylims <- c(nrow(plot.df), 0)
  y.mids <- (1:nrow(plot.df))-0.5

  # Vertically order cases and controls
  if(orient.cases == "top"){
    control.sep <- case.control.sep/2
    case.sep <- -case.control.sep/2
  }else{
    control.sep <- -case.control.sep/2
    case.sep <- case.control.sep/2
  }

  # Prep plot area
  prep.plot.area(xlims, ylims, parmar=parmar, xaxs="i", yaxs="r")
  axis(2, at=c(-10e10, 10e10), tck=0, labels=NA)
  abline(v=plot.df["control", 1], col=control.colors["main"])
  if(add.top.axis){
    clean.axis(3, label.units=top.axis.units,
               title=title, infinite.positive=TRUE,
               label.line=-0.8, title.line=0.1, max.ticks=5)
  }
  axis(2, at=y.mids, tick=F, line=-0.85, las=2, cex.axis=5.5/6,
       labels=cancer.names.short[rownames(plot.df)])

  # Add control bars
  rect(ybottom=y.mids-(bar.hex/2)+control.sep,
       ytop=y.mids+(bar.hex/2)+control.sep,
       xleft=0, xright=plot.df[, 4],
       col=control.colors[["main"]], border=NA, bty="n")
  segments(y0=y.mids+control.sep,
           y1=y.mids+control.sep,
           x0=plot.df[, 5], x1=plot.df[, 6], lwd=2, lend="butt",
           col=cancer.palettes[["control"]]["dark1"])
  rect(ybottom=y.mids-(bar.hex/2)+control.sep,
       ytop=y.mids+(bar.hex/2)+control.sep,
       xleft=0, xright=plot.df[, 4],
       col=NA, xpd=T)

  # Add bars, CI ticks, and outer borders
  sig.idx <- which(plot.df[, 7] < 0.05)
  bar.cols <- sapply(rownames(plot.df), function(pheno){cancer.palettes[[pheno]]["light1"]})
  ci.cols <- cancer.colors[rownames(plot.df)]
  if(color.by.sig){
    if(length(sig.idx) > 0){
      bar.cols[sig.idx] <- cancer.colors[rownames(plot.df)[sig.idx]]
      ci.cols[sig.idx] <- sapply(rownames(plot.df)[sig.idx], function(pheno){cancer.palettes[[pheno]]["dark1"]})
    }
  }
  rect(ybottom=y.mids-(bar.hex/2)+case.sep,
       ytop=y.mids+(bar.hex/2)+case.sep,
       xleft=0, xright=plot.df[, 1],
       col=bar.cols,
       border=NA, bty="n")
  segments(y0=y.mids+case.sep,
           y1=y.mids+case.sep,
           x0=plot.df[, 2], x1=plot.df[, 3], lwd=2, lend="butt",
           col=ci.cols)
  rect(ybottom=y.mids-(bar.hex/2)+case.sep,
       ytop=y.mids+(bar.hex/2)+case.sep,
       xleft=0, xright=plot.df[, 1],
       col=NA, xpd=T)

  # Add P values
  sapply(1:nrow(plot.df), function(i){
    if(color.by.sig){
      pval <- plot.df[i, 7]
      p.col <- if(pval < 0.05){"black"}else{control.colors[["main"]]}
    }else{
      p.col <- "black"
    }
    axis(4, at=y.mids[i], tick=F, line=-0.9, las=2, cex.axis=5/6,
         labels=PedSV::format.pval(plot.df[i, 7], nsmall=0),
         col.axis=p.col)
  })
}

