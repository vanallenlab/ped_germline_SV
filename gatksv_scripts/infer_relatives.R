#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Infer relatives from PC-Relate metrics


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(PedSV, quietly=TRUE)
require(argparse, quietly=TRUE)
require(e1071, quietly=TRUE)
PedSV::load.constants("all")
metric.labels <- c("kin" = "Kinship Coefficient",
                   "ibd0" = "Proportion IBD0",
                   "ibd1" = "Proportion IBD1",
                   "ibd2" = "Proportion IBD2")


##################
# Data functions #
##################
# Load training labels
load.training.labels <- function(path){
  labs <- read.table(path, header=T, sep="\t", comment.char="")
  colnames(labs)[3] <- "rel"
  rownames(labs) <- apply(labs[, 1:2], 1,
                          function(pair){paste(sort(pair), collapse="|")})
  labs[, -c(1:2), drop=FALSE]
}

# Invert label pairing
mirror.pair <- function(pair){
  parts <- unlist(strsplit(pair, split="|", fixed=T))
  paste(parts[2], parts[1], collapse="|")
}

# Identify & exclude outliers from training set
prune.outliers <- function(train, n.iqr=5){
  outliers <- unique(unlist(sapply(unique(train$rel), function(label){
    r.idxs <- which(train$rel == label)
    unlist(sapply(setdiff(colnames(train), "rel"), function(metric){
      vals <- as.numeric(train[r.idxs, metric])
      med <- median(vals, na.rm=T)
      iqr <- abs(diff(quantile(vals, probs=c(0.25, 0.75))))
      bounds <- med + (c(-n.iqr, n.iqr) * iqr)
      outliers <- which(vals < bounds[1] | vals > bounds[2])
      rownames(train[r.idxs, ])[outliers]
    }))
  })))
  train[setdiff(rownames(train), outliers), ]
}

# Train SVM classifier
train.classifier <- function(train, seed=2023){
  set.seed(seed)
  train$rel <- factor(train$rel)
  svm(rel ~ ., data=train,
      type="C-classification",
      kernel="linear",
      cross=10,
      probability=T,
      na.action=na.omit)
}

# Assign relationship labels
assign.labels <- function(kdf, classifier){
  predictions <- predict(classifier, kdf)
  predicted.labels <- as.character(predictions)
  names(predicted.labels) <- rownames(kdf)
  return(predicted.labels)
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Infer related samples")
parser$add_argument("--metrics", metavar=".tsv", type="character",
                    help="PC-Relate output .tsv", required=TRUE)
parser$add_argument("--training-labels", metavar=".tsv", type="character",
                    help=paste("three-column .tsv linking pairs of sample IDs to",
                               "asserted ground truth relationship labels. Used",
                               "for training."),
                    required=TRUE)
parser$add_argument("--out-prefix", metavar="path", type="character", required=TRUE,
                    help="path/prefix for all output files")
parser$add_argument("--plot", action="store_true", help="generate diagnostic plots")
args <- parser$parse_args()

# # DEV
# setwd("/Users/collins/Desktop/Collins/VanAllen/pediatric/riaz_pediatric_SV_collab/data/ancestry_and_relatedness")
# args <- list("metrics" = "PedSV.merged.kinship.tsv.gz",
#           "training_labels" = "PedSV.merged.known_relationships.tsv.gz",
#           "out_prefix" = "~/scratch/kinship_test",
#           "plot" = TRUE)

# Load PC-Relate metrics
kdf <- load.kinship.metrics(args$metrics)

# Load training data
train.labels <- load.training.labels(args$training_labels)
train <- merge(kdf, train.labels, by="row.names", sort=F, all=F)
rownames(train) <- train$Row.names
train$Row.names <- NULL
train.recip <- sapply(rownames(train), mirror.pair)

# Prune outliers from training data
train <- prune.outliers(train)

# Train classifier
classifier <- train.classifier(train)

# Apply classifier
pred.labels <- assign.labels(kdf, classifier)

# Evaluate accuracy vs. training labels
train.acc <- length(which(train$rel == pred.labels[rownames(train)])) / nrow(train)
cat(paste("\nClassification accuracy vs. training labels: ",
          formatC(round(100 * train.acc, 2), digits=4), "%\n\n", sep=""))

# Write labels to file
write.table(pred.labels,
            paste(args$out_prefix, "relatedness_labels.tsv", sep="."),
            row.names=T, col.names=F, sep="\t", quote=F)

# If optioned, generate diagnostic plots
png.dim <- 3.5
plot.legend.vals <- relative.colors
names(plot.legend.vals) <- relative.names
if(args$plot){
  apply(combn(colnames(kdf), 2), 2, function(col.pair){
    # All samples
    png(paste(args$out_prefix, ".all_pairs.", col.pair[1], "_vs_",
              col.pair[2], ".png", sep=""), res=300, height=png.dim*300,
        width=png.dim*300)
    pc.scatterplot(kdf,
                   which(colnames(kdf) == col.pair[1]),
                   which(colnames(kdf) == col.pair[2]),
                   colors=relative.colors[pred.labels[rownames(kdf)]],
                   title="All Pairs w/Predicted Labels",
                   x.title=metric.labels[col.pair[1]],
                   y.title=metric.labels[col.pair[2]],
                   legend.vals=plot.legend.vals)
    dev.off()

    # Training samples colored by training labels
    png(paste(args$out_prefix, ".training_pairs.", col.pair[1], "_vs_",
              col.pair[2], ".png", sep=""), res=300, height=png.dim*300,
        width=png.dim*300)
    pc.scatterplot(train,
                   which(colnames(train) == col.pair[1]),
                   which(colnames(train) == col.pair[2]),
                   colors=relative.colors[train$rel],
                   title="Training Pairs w/Training Labels",
                   x.title=metric.labels[col.pair[1]],
                   y.title=metric.labels[col.pair[2]],
                   legend.vals=plot.legend.vals)
    dev.off()

    # Non-training samples with inferred labels
    non.train <- kdf[setdiff(rownames(kdf), c(rownames(train), mirror.pair)), ]
    png(paste(args$out_prefix, ".new_pairs.", col.pair[1], "_vs_",
              col.pair[2], ".png", sep=""), res=300, height=png.dim*300,
        width=png.dim*300)
    pc.scatterplot(non.train,
                   which(colnames(non.train) == col.pair[1]),
                   which(colnames(non.train) == col.pair[2]),
                   colors=relative.colors[pred.labels[rownames(non.train)]],
                   title="New Pairs w/Predicted Labels",
                   x.title=metric.labels[col.pair[1]],
                   y.title=metric.labels[col.pair[2]],
                   legend.vals=plot.legend.vals)
    dev.off()
  })
}

