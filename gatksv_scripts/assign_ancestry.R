#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Assign sample ancestry based on genetic principal components


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(PedSV, quietly=TRUE)
require(argparse, quietly=TRUE)
require(e1071, quietly=TRUE)
PedSV::load.constants("all")


##################
# Data functions #
##################
# Load training labels
load.training.labels <- function(path){
  labs <- read.table(path, header=T, sep="\t", comment.char="")
  rownames(labs) <- labs[, 1]
  labs[, 1] <- NULL
  return(labs)
}

# Train SVM classifier
train.classifier <- function(train, seed=2023){
  set.seed(seed)
  train$pop <- factor(train$pop)
  svm(pop ~ ., data=train,
      type="C-classification",
      kernel="radial",
      cross=10,
      probability=T,
      na.action=na.omit)
}

# Assign ancestry labels
assign.ancestries <- function(pcs, classifier, min.prob=0.8){
  predictions <- predict(classifier, pcs, probability=T)
  predicted.labels <- as.character(predictions)
  names(predicted.labels) <- rownames(pcs)
  prediction.pvals <- attr(predictions, "probabilities")
  max.pred.pvals <- apply(prediction.pvals, 1, max)
  predicted.labels[which(max.pred.pvals < min.prob)] <- "OTH"
  return(predicted.labels)
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Assign sample ancestry")
parser$add_argument("--PCs", metavar=".tsv", type="character",
                    help="principal components .tsv", required=TRUE)
parser$add_argument("--training-labels", metavar=".tsv", type="character",
                    help=paste("two-column .tsv linking sample IDs to asserted",
                               "ground truth ancestry labels. Used for training."),
                    required=TRUE)
parser$add_argument("--out-prefix", metavar="path", type="character", required=TRUE,
                    help="path/prefix for all output files")
parser$add_argument("--use-N-PCs", default=10, metavar="integer", type="integer",
                    help="number of PCs to use [default: 10]")
parser$add_argument("--min-probability", metavar="float", type="double", default=0.8,
                    help=paste("minimum prediction probability required to assign",
                               "a sample to an ancestry [default: 0.8]"))
parser$add_argument("--plot", action="store_true", help="generate diagnostic plots")
args <- parser$parse_args()

# Load PCs
pcs <- load.pc.matrix(args$PCs, args$use_N_PCs)

# Load training data
train.labels <- load.training.labels(args$training_labels)
train <- merge(pcs, train.labels, by="row.names", sort=F, all=F)
rownames(train) <- train$Row.names
train$Row.names <- NULL

# Train classifier
classifier <- train.classifier(train)

# Apply classifier
pred.labels <- assign.ancestries(pcs, classifier, args$min_probability)

# Evaluate accuracy vs. training labels
acc <- length(which(train$pop == pred.labels[rownames(train)])) / nrow(train)
cat(paste("\nClassification accuracy vs. training labels: ",
            formatC(round(100 * acc, 2), digits=4), "%\n\n", sep=""))

# Write labels to file
write.table(pred.labels,
            paste(args$out_prefix, "ancestry_labels.tsv", sep="."),
            row.names=T, col.names=F, sep="\t", quote=F)

# If optioned, generate diagnostic plots
png.dim <- 3.5
if(args$plot){
  sapply(list(c(1, 2), c(3, 4)), function(pc.idxs){
    # All samples
    png(paste(args$out_prefix, ".all_samples.pc", pc.idxs[1], "_vs_pc",
              pc.idxs[2], ".png", sep=""), res=300, height=png.dim*300,
        width=png.dim*300)
    pc.scatterplot(pcs, pc.idxs[1], pc.idxs[2],
                   colors=pop.colors[pred.labels[rownames(pcs)]],
                   title="All Samples w/Predicted Labels",
                   legend.vals=pop.colors)
    dev.off()

    # Training samples colored by training labels
    png(paste(args$out_prefix, ".training_samples.pc", pc.idxs[1], "_vs_pc",
              pc.idxs[2], ".png", sep=""), res=300, height=png.dim*300,
        width=png.dim*300)
    pc.scatterplot(train, pc.idxs[1], pc.idxs[2], colors=pop.colors[train$pop],
                   title="Training Samples w/Training Labels",
                   legend.vals=pop.colors)
    dev.off()

    # Non-training samples with inferred labels
    non.train <- pcs[setdiff(rownames(pcs), rownames(train)), ]
    png(paste(args$out_prefix, ".new_samples.pc", pc.idxs[1], "_vs_pc",
              pc.idxs[2], ".png", sep=""), res=300, height=png.dim*300,
        width=png.dim*300)
    pc.scatterplot(non.train, pc.idxs[1], pc.idxs[2],
                   colors=pop.colors[pred.labels[rownames(non.train)]],
                   title="New Samples w/Predicted Labels",
                   legend.vals=pop.colors)
    dev.off()
  })
}

