#!/usr/bin/env Rscript

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Diploid-normalize an allele dosage matrix
# Reads to/from stdin/stdout

options(scipen=1000, stringsAsFactors=F)

ad <- read.table(file("stdin"), header=F, sep="\t")

ad[, -c(1:4)] <- ad[, -c(1:4)] - 2

write.table(ad, "", col.names=F, row.names=F, sep="\t", quote=F)
