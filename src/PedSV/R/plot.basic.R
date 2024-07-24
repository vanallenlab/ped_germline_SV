#!/usr/bin/env R

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
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
#' @param x.label.line Line for X-axis labels (`label.line` parameter for [PedSV::clean.axis])
#' @param x.title Title for X-axis
#' @param x.title.line Line for X-axis title (`title.line` parameter for [PedSV::clean.axis])
#' @param xlims Limits for X-axis
#' @param y.label.line Line for X-axis labels (`label.line` parameter for [PedSV::clean.axis])
#' @param y.title Title for Y-axis
#' @param y.title.line Line for axis titles (`title.line` parameter for [PedSV::clean.axis])
#' @param ylims Limits for Y-axis
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
pc.scatterplot <- function(pcs, pc.X, pc.Y, colors, title=NULL,
                           x.label.line=NULL, x.title=NULL, x.title.line=0.5, xlims=NULL,
                           y.label.line=NULL, y.title=NULL, y.title.line=0.5, ylims=NULL,
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
  if(is.null(xlims)){
    xlims <- range(x)
  }
  if(is.null(ylims)){
    ylims <- range(y)
  }

  prep.plot.area(xlims, ylims, parmar=parmar, xaxs="r", yaxs="r")
  mtext(3, text=title)

  if(is.null(x.title)){
    x.title <- paste("Principal component", pc.X)
  }
  clean.axis(1, title=x.title, infinite=T, label.line=x.label.line, title.line=x.title.line)

  if(is.null(y.title)){
    y.title <- paste("Principal component", pc.Y)
  }
  clean.axis(2, title=y.title, infinite=T, label.line=y.label.line, title.line=y.title.line)

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
           pt.cex=1.5, bty="n", border=NA, cex=5/6)
  }
}


#' Barplot of values per cancer type
#'
#' Generate a barplot of one value for each cancer type
#'
#' @param plot.df Data.frame of values to plot. See `Details`.
#' @param bar.hex Relative width for each bar \[default: 0.5\]
#' @param case.control.sep Relative separation between overlapping case
#' and control bars \[default: 0.375\]
#' @param color.by.sig Should bars be shaded differently by significance
#' level? \[default: TRUE\]
#' @param add.top.axis Should the top axis be added? \[default: TRUE\]
#' @param top.axis.units Specify units for top axis. Options are NULL (for
#' numeric) or "percent" \[default: NULL\]
#' @param title Custom title for top axis
#' @param orient.cases Should casees be plotted on the `top` or `bottom` of controls? \[default: "top"\]
#' @param custom.pheno.labels Custom phenotype labels for the groups on the
#' Y axis, if desired.
#' @param legend.on.bars Should "case" and "control" labels be printed on the
#' largest bars? \[default: FALSE\]
#' @param case.label Label for "case" bars if `legend.on.bars` is `TRUE`
#' \[default: "Case"\]
#' @param control.label Label for "control" bars if `legend.on.bars` is
#' `TRUE` \[default: "Control"\]
#' @param add.pvals Should P values be annotated on the opposite margin?
#' \[default: TRUE\]
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
barplot.by.phenotype <- function(plot.df, bar.hex=0.5, case.control.sep=0.375,
                                 color.by.sig=TRUE, add.top.axis=TRUE, top.axis.units=NULL,
                                 title="Value", orient.cases="top",
                                 custom.pheno.labels=NULL, legend.on.bars=FALSE,
                                 case.label="Case", control.label="Control",
                                 add.pvals=TRUE, parmar=c(0.2, 4.1, 2.1, 4)){
  # Get plot dimensions
  xlims <- c(0, min(c(2*max(plot.df[, c(1, 4)], na.rm=T),
                      max(plot.df[, 1:6], na.rm=T)),
                    na.rm=T))
  ylims <- c(nrow(plot.df), 0)
  y.mids <- (1:nrow(plot.df))-0.5

  # Vertically order cases and controls
  if(orient.cases == "top"){
    control.sep <- case.control.sep/2
  }else{
    control.sep <- -case.control.sep/2
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
  y.labels <- if(!is.null(custom.pheno.labels)){
    custom.pheno.labels
  }else{
    cancer.names.short[rownames(plot.df)]
  }
  axis(2, at=y.mids, tick=F, line=-0.85, las=2, cex.axis=5.5/6, labels=y.labels)

  # Add bars
  add.pheno.bars(plot.df, bar.mids=y.mids, bar.hex=bar.hex, control.sep=control.sep,
                 horiz=TRUE, color.by.sig=color.by.sig, legend.on.bars=legend.on.bars,
                 add.pvals=add.pvals, case.label=case.label, control.label=control.label)
}


#' Double-wide barpolot of values per cancer type
#'
#' Generate two side-by-side sets of barplots of two values for each cancer type
#'
#' @param plot.df.l Data.frame of values to plot on the left. See [barplot.by.phenotype()].
#' @param plot.df.r Data.frame of values to plot on the right. See [barplot.by.phenotype()].
#' @param bar.wex Relative width for each bar \[default: 0.5\]
#' @param case.control.sep Relative separation between overlapping case
#' and control bars \[default: 0.375\]
#' @param group.spacer Distance between left and right barplot groups \[default: 0.5\]
#' @param color.by.sig Should bars be shaded differently by significance
#' level? \[default: TRUE\]
#' @param add.left.axis Should the left axis be added? \[default: TRUE\]
#' @param left.axis.units Specify units for left axis. Options are NULL (for
#' numeric) or "percent" \[default: NULL\]
#' @param title Custom title for left axis
#' @param y.title.line `title.line` parameter passed to [PedSV::clean.axis()]
#' for left axis \[default: 0.35\]
#' @param orient.cases Should cases be plotted to the `left` or `right` of controls? \[default: "left"\]
#' @param label.l Optional label for left set of bars.
#' @param label.r Optional label for right set of bars.
#' @param group.label.cex Character expansion used for `label.l` and `label.r` \[default: 1\]
#' @param parmar Value of `mar` passed to `par()`
#'
#' @details This function is effectively a double-wide set of [barplot.by.phenotype()].
#' As such, please refer to [barplot.by.phenotype()] for more information.
#'
#' @seealso [barplot.by.phenotype()]
#'
#' @export doublewide.barplot.by.phenotype
#' @export
doublewide.barplot.by.phenotype <-
  function(plot.df.l, plot.df.r, bar.wex=0.5, case.control.sep=0.375,
           group.spacer=0.5, color.by.sig=TRUE, add.left.axis=TRUE,
           left.axis.units=NULL, title="Value", y.title.line=0.35,
           orient.cases="left", label.l=NULL, label.r=NULL, group.label.cex=1,
           parmar=c(1.15, 2.25, 0.75, 0.25)){
    # Get plot dimensions
    both.plot.df <- rbind(plot.df.l, plot.df.r)
    ylims <- c(0, 1.15*min(c(2*max(both.plot.df[, c(1, 4)], na.rm=T),
                             max(both.plot.df[, 1:6], na.rm=T)),
                           na.rm=T))
    xlims <- c(-0.5*group.spacer, nrow(both.plot.df) + (1.5*group.spacer))
    x.mids.l <- (1:nrow(plot.df.l))-0.5
    x.mids.r <- (1:nrow(plot.df.r))-0.5 + group.spacer + nrow(plot.df.l)
    group.x.mean.l <- mean(c(0, nrow(plot.df.l)))
    group.x.mean.r <- mean(c(0, nrow(plot.df.l))) + nrow(plot.df.l) + group.spacer

    # Horizontally order cases and controls
    if(orient.cases == "left"){
      control.sep <- case.control.sep/2
    }else{
      control.sep <- -case.control.sep/2
    }

    # Prep plot area
    prep.plot.area(xlims, ylims, parmar=parmar, xaxs="i", yaxs="r")
    axis(1, at=c(0, nrow(plot.df.l)), tck=0, labels=NA, col="gray70")
    axis(1, at=c(0, nrow(plot.df.r)) + nrow(plot.df.l) + group.spacer,
         tck=0, labels=NA, col="gray70")
    if(add.left.axis){
      clean.axis(2, label.units=left.axis.units,
                 title=title, infinite.positive=TRUE, tck=-0.015,
                 label.line=-0.8, title.line=y.title.line, max.ticks=5)
    }
    sapply(1:2, function(i){
      axis(1, at=c(group.x.mean.l, group.x.mean.r)[i], tick=F, line=-0.9,
           cex.axis=group.label.cex, labels=c(label.l, label.r)[i])
    })

    # Add bars
    add.pheno.bars(plot.df.l, bar.mids=x.mids.l,
                   bar.hex=bar.wex, control.sep=control.sep,
                   horiz=FALSE, color.by.sig=color.by.sig)
    add.pheno.bars(plot.df.r, bar.mids=x.mids.r,
                   bar.hex=bar.wex, control.sep=control.sep,
                   horiz=FALSE, color.by.sig=color.by.sig)

    # Add P value annotations
    p.spacer <- 0.05 * diff(par("usr")[3:4])
    sapply(1:nrow(plot.df.l), function(x){
      p <- plot.df.l[x, 7]
      xpos <- x.mids.l[x] - control.sep
      ypos <- max(plot.df.l[x, 1:6], na.rm=T) + p.spacer
      p.label <- if(p<0.0005){"***"}else if(p<0.005){"**"}else if(p<0.05){"*"}else{NULL}
      text(x=xpos, y=ypos, labels=p.label, font=2, xpd=T)
    })
    sapply(1:nrow(plot.df.r), function(x){
      p <- plot.df.r[x, 7]
      xpos <- x.mids.r[x] - control.sep
      ypos <- max(plot.df.r[x, 1:6], na.rm=T) + p.spacer
      p.label <- if(p<0.0005){"***"}else if(p<0.005){"**"}else if(p<0.05){"*"}else{NULL}
      text(x=xpos, y=ypos, labels=p.label, font=2, xpd=T)
    })
  }


#' Plot values vs. SV length
#'
#' Generic function to plot one or more (x, y) datasets versus SV length
#'
#' @param x.svlen List of SV length values to be plotted on X axis; one element per stratum
#' @param y.value List of values to be plotted on Y axis; one element per stratum
#' @param colors Vector of colors for the list elements in `x.svlen` and `y.value`
#' @param ci.lower (Optional) List of lower confidence interval bounds per `y.value`; one element per stratum
#' @param ci.upper (Optional) List of lower confidence interval bounds per `y.value`; one element per stratum
#' @param group.names (Optional) group names to assign to each list element  in `x.svlen` and `y.value`
#' @param step Should the line be rendered as a step function (a la Kaplan-Meyer)? \[default: TRUE\]
#' @param lwds (Optional) vector of `lwd` for each line \[default: 3\]
#' @param ci.alpha Transparency value `alpha` for confidence interval shading \[default: 0.15\]
#' @param legend Should a legend be plotted? \[default: FALSE\]
#' @param legend.names (Optional) mapping of `values` to labels for legend
#' @param legend.pos (Optional) position for legend. See [legend()] \[default: "topright"]
#' @param legend.label.spacing Minimum vertical spacing between legend labels \[default: 0.075\]
#' @param title (Optional) Title for plot
#' @param xlab (Optional) Title for X axis
#' @param xlab.line Line parameter for `xlab` \[default: 1\]
#' @param xlims (Optional) two-element vector of start and stop values for X-axis, in log10(SVLEN)
#' @param ylab (Optional) Title for Y axis
#' @param y.title.line Value for `title.line` passed to [PedSV::clean.axis()] \[default: 1\]
#' @param y.axis.units Specify units for Y axis. Options are NULL (for
#' numeric) or "percent" \[default: NULL\]
#' @param ylims (Optional) two-element vector of start and stop values for Y-axis
#' @param x.axis.labels (Optional) custom labels for X axis
#' @param x.axis.labels.at (Optional) custom positions for X axis labels
#' @param x.tck Value of `tck` passed to [clean.axis()]
#' @param parmar Margin values passed to par()
#'
#' @export svlen.line.plot
#' @export
svlen.line.plot <- function(x.svlen, y.value, colors, ci.lower=NULL, ci.upper=NULL,
                            group.names=NULL, step=TRUE, lwds=NULL, ci.alpha=NULL,
                            legend=FALSE, legend.names=NULL, legend.pos="topright",
                            legend.label.spacing=0.075, title=NULL, xlab=NULL, xlab.line=0,
                            xlims=log10(c(10000, 1000000)), ylab="Value", y.title.line=1,
                            y.axis.units=NULL, ylims=NULL, x.axis.labels=NULL,
                            x.axis.labels.at=NULL, x.tck=-0.015, parmar=c(2, 3, 1, 1)){
  # Ensure PedSV scale constants are loaded within function scope
  PedSV::load.constants("scales", envir=environment())

  # Get plotting values
  if(is.null(group.names)){
    group.names <- names(x.svlen)
  }
  n.groups <- length(x.svlen)
  if(is.null(lwds)){
    lwds <- rep(1, n.groups)
  }
  if(is.null(ci.alpha)){
    ci.alpha <- rep(0.15, n.groups)
  }
  if(is.null(legend.names)){
    legend.names <- names(x.svlen)
  }
  if(is.null(xlab)){
    xlab <- "Size of Largest SV"
  }
  if(is.null(ylab)){
    ylab <- "Value"
  }
  if(is.null(xlims)){
    xlims <- c(0, max(unlist(x.svlen), na.rm=T))
  }
  if(is.null(ylims)){
    ylims <- range(unlist(y.value), na.rm=T)
  }
  if(is.null(x.axis.labels)){
    x.axis.labels <- logscale.demi.bp.labels
  }
  if(is.null(x.axis.labels.at)){
    x.axis.labels.at <- log10(logscale.demi.bp)
  }

  # Prep plot area
  prep.plot.area(xlims, ylims, parmar)

  # Add confidence intervals
  if(!is.null(ci.lower) & !is.null(ci.upper)){
    # Loop over this twice: first to lay white backgrounds, then add colors
    for(layer in c("white", "colors")){
      sapply(1:n.groups, function(i){
        n.times <- length(x.svlen[[i]])
        if(n.times > 1){
          if(step){
            x.bottom <- c(0, PedSV::stretch.vector(x.svlen[[i]], 2)[-2*n.times])
            x.top <- rev(x.bottom)
            y.bottom <- c(1, 1, PedSV::stretch.vector(ci.lower[[i]], 2)[-c(2*n.times-c(0, 1))])
            y.top <- rev(c(1, 1, PedSV::stretch.vector(ci.upper[[i]], 2)[-c(2*n.times-c(0, 1))]))
          }else{
            x.bottom <- x.svlen[[i]]
            x.top <- rev(x.svlen[[i]])
            y.bottom <- ci.lower[[i]]
            y.top <- rev(ci.upper[[i]])
          }
          polygon(x=c(x.bottom, x.top), y=c(y.bottom, y.top), border=NA, bty="n",
                  col=if(layer == "white"){"white"}else{adjustcolor(colors[[i]], alpha=ci.alpha[i])})
        }
      })
    }
  }

  # Add lines
  sapply(1:n.groups, function(i){
    n.times <- length(x.svlen[[i]])
    if(n.times > 0){
      if(step){
        x <- c(0, PedSV::stretch.vector(x.svlen[[i]], 2))
        y <- c(1, 1, PedSV::stretch.vector(y.value[[i]], 2))[1:length(x)]
      }else{
        x <- x.svlen[[i]]
        y <- y.value[[i]]
      }
      points(x, y, type="l", col=colors[[i]], lwd=lwds[i])
    }
  })

  # Add axes
  clean.axis(1, at=log10(logscale.minor), labels=NA, infinite=TRUE,
             title=xlab, label.line=-0.75, title.line=xlab.line, tck=x.tck)
  clean.axis(1, at=x.axis.labels.at, labels=x.axis.labels,
             infinite=FALSE, title=NA, label.line=-0.75, title.line=0, tck=-0.0225)
  clean.axis(2, title=ylab, infinite=TRUE, tck=-0.0175,
             label.units=y.axis.units, title.line=y.title.line)
  mtext(title, side=3, line=0)

  # Add legend
  if(legend){
    legend(legend.pos, legend=rev(legend.names), lwd=rev(lwds),
           col=rev(colors), bty="n", cex=5/6)
  }
}
