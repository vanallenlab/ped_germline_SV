#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Identify outlier samples based on SV counts


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(PedSV, quietly=TRUE)
require(beeswarm, quietly=TRUE)
require(argparse, quietly=TRUE)
PedSV::load.constants("all")


##################
# Data functions #
##################
# Load population labels
load.pop.labels <- function(pops.in){
  df <- read.table(pops.in, sep="\t", header=F)
  pops <- df[, 2]
  names(pops) <- df[, 1]
  return(pops)
}

# Load counts, group by SV type, and add population labels
load.counts <- function(counts.in, pop.labels, smallest.max=10){
  # Read data
  counts <- read.table(counts.in, sep="\t", header=F)
  colnames(counts) <- c("sample", "svtype", "count")

  # Add ancestry labels
  counts$pop <- pop.labels[counts$sample]
  missing.pop <- which(is.na(counts$pop))
  if(length(missing.pop) > 0){
    counts$pop[missing.pop] <- "OTH"
  }

  # Drop SV types with fewer than smallest.max SVs in all samples
  all.svtypes <- sort(unique(counts$svtype))
  max.by.svt <- sapply(all.svtypes, function(svt){
    max(counts$count[which(counts$svtype == svt)], na.rm=T)
  })
  low.n.idx <- which(max.by.svt < smallest.max)
  if(length(low.n.idx) > 0){
    counts <- counts[-which(counts$svtype %in% names(low.n.idx)), ]
  }

  # Split by SV type and pop
  svtypes <- sort(unique(counts$svtype))
  counts <- lapply(svtypes, function(svtype){
    subcounts <- counts[which(counts$svtype == svtype), ]
    pop.counts <- lapply(names(pop.colors), function(pop){
      subcounts[which(subcounts$pop == pop), ]
    })
    names(pop.counts) <- names(pop.colors)
    return(pop.counts)
  })
  names(counts) <- svtypes
  return(counts)
}

# Calculate cutoffs on an ancestry-specific basis
calculate.cutoffs <- function(counts, MAD=5){
  lapply(counts, function(clist){
    res.by.pop <- lapply(clist, function(df){
      vals <- as.numeric(df$count)
      d <- MAD * mad(vals, na.rm=T)
      lower <- quantile(vals, prob=0.25) - d
      upper <- quantile(vals, prob=0.75) + d
      list("df" = df, "cutoffs" = c(lower, upper), "MAD" = d)
    })

    # Update OTH population to use stricter cutoffs
    non.OTH.pops <- setdiff(names(pop.colors), "OTH")
    pop.lower <- min(sapply(res.by.pop[non.OTH.pops],
                        function(plist){plist$cutoffs[1]}), na.rm=T)
    pop.upper <- max(sapply(res.by.pop[non.OTH.pops],
                            function(plist){plist$cutoffs[2]}), na.rm=T)
    oth.cutoffs <- res.by.pop[["OTH"]][["cutoffs"]]
    new.oth.cutoffs <- c(max(pop.lower, oth.cutoffs[1]),
                         min(pop.upper, oth.cutoffs[2]))
    res.by.pop[["OTH"]][["cutoffs"]] <- new.oth.cutoffs
    return(res.by.pop)
  })
}

# Define outliers per condition using precomputed cutoffs
define.outliers <- function(counts.anno){
  lapply(counts.anno, function(clist){
    res.by.pop <- lapply(clist, function(plist){
      df <- plist$df
      cutoffs <- plist$cutoffs
      vals <- as.numeric(df$count)
      if(plist$MAD > 0){
        out.idx <- which(vals < cutoffs[1] | vals > cutoffs[2])
      }else{
        out.idx <- c()
      }
      df$outlier <- FALSE
      if(length(out.idx) > 0){
        outliers <- df$sample[out.idx]
        df$outlier[out.idx] <- TRUE
      }else{
        outliers <- c()
      }
      list("df" = df, "outliers" = outliers, "cutoffs" = cutoffs)
    })
  })
}


# Build a data frame of outlier samples and the reasons they were excluded
gather.outliers <- function(counts.anno){
  outlier.df <- data.frame("sample"=character(), "failures"=character())
  for(k in names(counts.anno)){
    for(pop in names(counts.anno[[k]])){
      for(sid in counts.anno[[k]][[pop]]$outliers){
        if(sid %in% outlier.df$sample){
          hit.idx <- which(outlier.df$sample == sid)
          outlier.df$failures[hit.idx] <- paste(outlier.df$failures[hit.idx], k, sep=";")
        }else{
          outlier.df <- rbind(outlier.df, data.frame("sample"=sid, "failures"=k))
        }
      }
    }
  }
  colnames(outlier.df)[1] <- "#sample"
  outlier.df[order(sapply(strsplit(outlier.df$failures, split=";"), length), decreasing=T), ]
}


######################
# Plotting functions #
######################
# Plot counts by ancestry
plot.counts <- function(count.list, title){
  # Get plot dimensions
  n.pops <- length(count.list)
  ylims <- range(do.call("rbind", lapply(count.list, function(l){l$df}))$count, na.rm=T)

  # Prepare plot
  prep.plot.area(c(0, n.pops), ylims, parmar=c(1.25, 4, 1, 0.25), yaxs="r")
  clean.axis(1, at=(1:n.pops) - 0.5, tck=0, labels=names(count.list), infinite=T)
  clean.axis(2, title="SV Count", infinite=T, title.line=1.5)
  mtext(title, 3, font=2)

  # Add one swarm for each population
  sapply(1:length(count.list), function(i){
    df <- count.list[[i]]$df
    vals <- df$count
    pwpch <- c("TRUE"=4, "FALSE"=20)[as.character(df$outlier)]
    # pwcex <- c("TRUE"=0.2, "FALSE"=0.3)[as.character(df$outlier)]
    beeswarm(vals, at=i-0.5, corral="wrap", corralWidth=0.9,
             col=pop.colors[names(count.list[i])],
             pwpch=pwpch, cex=0.15, add=T)
    segments(x0=rep(i-0.85, 2), x1=rep(i-0.15, 2),
             y0=count.list[[i]]$cutoffs, y1=count.list[[i]]$cutoffs,
             col="gray70")
  })
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Determine outlier samples")
parser$add_argument("counts", metavar="counts.tsv", type="character",
                    help="SV count .tsv")
parser$add_argument("pop_labels", metavar="populations.tsv", type="character",
                    help="Two-column .tsv of sample ID & population")
parser$add_argument("--MAD", metavar="numeric", type="numeric", default=5,
                    help="Number of median absolute deviations for defining outliers")
parser$add_argument("--out-prefix", metavar="path", type="character", default="./",
                    help="path/prefix for all output files [default: pwd]")
args <- parser$parse_args()

# # DEV:
# args <- list("counts" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/data/outliers/v2.5.1/PedSV.v2.5.1.all_outlier_counts.tsv.gz",
#              "pop_labels" = "~/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/data/outliers/v2.5.1/PedSV.v2.best_guess_ancestries.for_outlier_exclusion.tsv",
#              "MAD" = 5,
#              "out_prefix" = "~/scratch/PedSV.v2.5.1.")

# Read population labels
pop.labels <- load.pop.labels(args$pop_labels)

# Load and format counts
counts <- load.counts(args$counts, pop.labels)

# Define outliers
counts.anno <- calculate.cutoffs(counts, args$MAD)
counts.anno <- define.outliers(counts.anno)

# Compile list of all outlier samples across all conditions
outlier.df <- gather.outliers(counts.anno)
write.table(outlier.df, paste(args$out_prefix, "all_outliers.tsv", sep=""),
            sep="\t", col.names=T, row.names=F, quote=F)

# Generate plots of each condition
for(k in names(counts.anno)){
  pdf(paste(args$out_prefix, paste(k, "counts_by_pop.pdf", sep="."), sep=""),
      height=3, width=4)
  plot.counts(counts.anno[[k]], k)
  dev.off()
}
