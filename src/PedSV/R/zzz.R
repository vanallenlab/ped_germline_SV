#!/usr/bin/env R

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>


# Ensure non-standard dependencies are available when package is loaded


.onLoad <- function(libname, pkgname){

  # Set useful global constants
  options(scipen=1000, stringsAsFactors=F, family="sans")

  # Check to make sure RLCtools is available
  if(!require(RLCtools)){
    stop(paste("Dependency `RLCtools` is required for `PedSV` but is not found.\n",
              "For more info, see https://github.com/RCollins13/RLCtools\n", sep=""))
  }
}
