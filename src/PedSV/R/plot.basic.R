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
#' @param x.title Title for X-axis
#' @param x.title.line Line for X-axis title (`title.line` parameter for [PedSV::clean.axis])
#' @param xlims Limits for X-axis
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
                           x.title=NULL, x.title.line=0.5, xlims=NULL,
                           y.title=NULL, y.title.line=0.5, ylims=NULL,
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
           pt.cex=1.5, bty="n", border=NA, cex=5/6)
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
#' @param case.control.sep Relative separation between overlapping case
#' and control bars \[default: 0.375\]
#' @param color.by.sig Should bars be shaded differently by significance
#' level? \[default: TRUE\]
#' @param add.top.axis Should the top axis be added? \[default: TRUE\]
#' @param top.axis.units Specify units for top axis. Options are NULL (for
#' numeric) or "percent" \[default: NULL\]
#' @param title Custom title for top axis
#' @param custom.pheno.labels Custom phenotype labels for the groups on the
#' Y axis, if desired.
#' @param legend.on.bars Should "case" and "control" labels be printed on the
#' largest bars? \[default: FALSE\]
#' @param case.label Label for "case" bars if `legend.on.bars` is `TRUE`
#' \[default: "Case"\]
#' @param control.label Label for "control" bars if `legend.on.bars` is
#' `TRUE` \[default: "Control"\]
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
                                 parmar=c(0.2, 4.1, 2.1, 4)){
  # Get plot dimensions
  xlims <- c(0, min(c(2*max(plot.df[, c(1, 4)], na.rm=T),
                      max(plot.df[, 1:6], na.rm=T)),
                    na.rm=T))
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
  y.labels <- if(!is.null(custom.pheno.labels)){
    custom.pheno.labels
  }else{
    cancer.names.short[rownames(plot.df)]
  }
  axis(2, at=y.mids, tick=F, line=-0.85, las=2, cex.axis=5.5/6, labels=y.labels)

  # Add control bars
  rect(ybottom=y.mids-(bar.hex/2)+control.sep,
       ytop=y.mids+(bar.hex/2)+control.sep,
       xleft=0, xright=plot.df[, 4],
       col=control.colors[["main"]], border=NA, bty="n")
  segments(y0=y.mids+control.sep,
           y1=y.mids+control.sep,
           x0=plot.df[, 5], x1=plot.df[, 6], lwd=2, lend="butt",
           col=cancer.palettes[["control"]]["dark1"])
  if(legend.on.bars){
    longest.control <- head(which(plot.df[, 4] == max(plot.df[, 4], na.rm=T)), 1)
    text(x=0-(0.05*diff(par("usr")[1:2])),
         y=y.mids[longest.control]+(bar.hex/5)+control.sep,
         pos=4, label=control.label, cex=4/6, col=control.colors[["dark2"]])
  }
  rect(ybottom=y.mids-(bar.hex/2)+control.sep,
       ytop=y.mids+(bar.hex/2)+control.sep,
       xleft=0, xright=plot.df[, 4],
       col=NA, xpd=T)

  # Add bars, CI ticks, and outer borders
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
  rect(ybottom=y.mids-(bar.hex/2)+case.sep,
       ytop=y.mids+(bar.hex/2)+case.sep,
       xleft=0, xright=plot.df[, 1],
       col=bar.cols,
       border=NA, bty="n")
  segments(y0=y.mids+case.sep,
           y1=y.mids+case.sep,
           x0=plot.df[, 2], x1=plot.df[, 3], lwd=2, lend="butt",
           col=ci.cols)
  if(legend.on.bars){
    longest.case <- head(which(plot.df[, 1] == max(plot.df[, 1], na.rm=T)), 1)
    text(x=0-(0.05*diff(par("usr")[1:2])),
         y=y.mids[longest.case]+(bar.hex/10)+case.sep,
         pos=4, label=case.label, cex=4/6, col=bar.pals[[longest.case]][["dark2"]])
  }
  rect(ybottom=y.mids-(bar.hex/2)+case.sep,
       ytop=y.mids+(bar.hex/2)+case.sep,
       xleft=0, xright=plot.df[, 1],
       col=NA, xpd=T)

  # Add P values
  sapply(1:nrow(plot.df), function(i){
    if(color.by.sig){
      pval <- plot.df[i, 7]
      p.col <- if(!is.na(pval) & pval < 0.05){"black"}else{control.colors[["main"]]}
    }else{
      p.col <- "black"
    }
    p.label <- if(is.na(pval)){"NA"}else{PedSV::format.pval(plot.df[i, 7], nsmall=0)}
    axis(4, at=y.mids[i], tick=F, line=-0.9, las=2, cex.axis=5/6,
         labels=p.label, col.axis=p.col)
  })
}


#' Kaplan-Meyer plots of SV length
#'
#' Customized Kaplan-Meyer "survival" plots of one or more strata by longest SV
#' per sample
#'
#' @param surv.models List of one or more [`survival::summary.survfit`] objects
#' @param colors Vector of colors for the list elements in `surv.models`
#' @param group.names (Optional) group names to assign to each list element in `surv.models`
#' @param lwds (Optional) vector of `lwd` for each of surv.models \[default: 3\]
#' @param ci.alpha Transparency value `alpha` for confidence interval shading \[default: 0.15\]
#' @param legend Should a legend be plotted? \[default: FALSE\]
#' @param legend.names (Optional) mapping of `values` to labels for legend
#' @param legend.pos (Optional) position for legend. See [legend()] \[default: "topright"]
#' @param legend.label.spacing Minimum vertical spacing between legend labels \[default: 0.075\]
#' @param title (Optional) Title for plot
#' @param xlab (Optional) Title for X axis
#' @param xlims (Optional) two-element vector of start and stop values for X-axis, in log10(SVLEN)
#' @param ylab (Optional) Title for Y axis
#' @param ylims (Optional) two-element vector of start and stop values for Y-axis
#' @param parmar Margin values passed to par()
#'
#' @seealso [`survival::Surv`], [`survival::survfit`], [`survival::summary.survfit`]
#'
#' @export svlen.km.plot
#' @export
svlen.km.plot <- function(surv.models, colors, group.names=NULL, lwds=NULL, ci.alpha=0.15,
                          legend=FALSE, legend.names=NULL, legend.pos="topright",
                          legend.label.spacing=0.075, title=NULL, xlab=NULL,
                          xlims=log10(c(10000, 1000000)), ylab=NULL, ylims=NULL,
                          parmar=c(2, 3, 1, 1)){
  # Ensure survival library and PedSV scale constants are loaded within function scope
  require(survival, quietly=TRUE)
  PedSV::load.constants("scales", envir=environment())

  # Get plotting values
  if(is.null(group.names)){
    group.names <- names(surv.models)
  }
  n.groups <- length(surv.models)
  if(is.null(lwds)){
    lwds <- rep(1, n.groups)
  }
  if(is.null(legend.names)){
    legend.names <- names(surv.models)
  }
  if(is.null(xlab)){
    xlab <- "Size of Largest SV"
  }
  if(is.null(ylab)){
    ylab <- bquote("Samples with" >= 1 ~ "SV")
  }
  if(is.null(xlims)){
    xlims <- c(0, max(sapply(surv.models, function(ss){max(ss$time, na.rm=T)})))
  }
  if(is.null(ylims)){
    ylims <- c(0, max(sapply(surv.models, function(ss){max(ss$surv[which(ss$time>xlims[1])], na.rm=T)}), na.rm=T) + 0.025)
  }

  # Prep plot area
  prep.plot.area(xlims, ylims, parmar)

  # Add confidence intervals
  # Loop over this twice: first to lay white backgrounds, then add colors
  for(layer in c("white", "colors")){
    sapply(1:n.groups, function(i){
      n.times <- length(surv.models[[i]]$time)
      if(n.times > 1){
        x.bottom <- c(0, PedSV::stretch.vector(surv.models[[i]]$time, 2)[-2*n.times])
        x.top <- rev(x.bottom)
        y.bottom <- c(1, 1, PedSV::stretch.vector(surv.models[[i]]$lower, 2)[-c(2*n.times-c(0, 1))])
        y.top <- rev(c(1, 1, PedSV::stretch.vector(surv.models[[i]]$upper, 2)[-c(2*n.times-c(0, 1))]))
        polygon(x=c(x.bottom, x.top), y=c(y.bottom, y.top), border=NA, bty="n",
                col=if(layer == "white"){"white"}else{adjustcolor(colors[[i]], alpha=ci.alpha)})
      }
    })
  }

  # Add K-M curves
  sapply(1:n.groups, function(i){
    n.times <- length(surv.models[[i]]$time)
    # If summary.survfit returns no data, this is either because
    # there are no patients in this group or nobody died.
    # If the latter, we can plot as a flat line at Y=1 until rmean.endtime (I think?)
    if(surv.models[[i]]$n > 0){
      if(n.times == 0){
        x <- c(0, surv.models[[i]]$rmean.endtime)
        y <- c(1, 1)
      }else{
        x <- c(0, PedSV::stretch.vector(surv.models[[i]]$time, 2))
        y <- c(1, 1, PedSV::stretch.vector(surv.models[[i]]$surv, 2))[1:length(x)]
      }
      points(x, y, type="l", col=colors[[i]], lwd=lwds[i])
    }
  })

  # Add axes
  clean.axis(1, at=log10(logscale.minor), labels=NA, infinite=TRUE,
             title=xlab, label.line=-0.75, title.line=0, tck=-0.01)
  clean.axis(1, at=log10(logscale.demi.bp), labels=logscale.demi.bp.labels,
             infinite=FALSE, title=NA, label.line=-0.75, title.line=0, tck=-0.0225)
  clean.axis(2, title=ylab, infinite=TRUE, tck=-0.0175,
             label.units="percent", title.line=1)
  mtext(title, side=3, line=0, font=2)

  # Add legend
  if(legend){
    legend(legend.pos, legend=rev(legend.names), lwd=rev(lwds),
           col=rev(colors), bty="n", cex=5/6)
  }
}
