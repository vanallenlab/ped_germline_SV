#!/usr/bin/env R

###############################
#    Study of Germline SVs    #
#    in Pediatric Cancers     #
###############################

# Copyright (c) 2023-Present Ryan L. Collins, Riaz Gillani, Jett Crowdis, and the Van Allen Laboratory
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Project-wide constants


#' Load Study Constants
#'
#' Load a subset of constants for germline SV analyses in pediatric cancers
#'
#' @param susbet Vector of constant groups to load See `Details` for options.
#' @param envir Environment passed to [base::assign] \[default: .GlobalEnv\]
#'
#' @details Recognized values for `subset` include:
#' * `colors` : all color palettes used
#' * `scales` : all scales and scale labels
#' * `names` : names of various variables
#' * `all` : load all constants
#'
#' @examples
#' # Load list of color palettes
#' get.constants("colors");
#'
#' # Load scales and names colors
#' get.constants(c("scales", "names"))
#'
#' @seealso [base::assign]
#'
#' @export load.constants
#' @export
load.constants <- function(subset, envir=.GlobalEnv){
  # Define colors
  pancan.colors <- c("dark3" = "#290A0D",
                     "dark2" = "#52141A",
                     "dark1" = "#7A1D26",
                     "main" = "#A32733",
                     "light1" = "#B5525C",
                     "light2" = "#C87D85",
                     "light3" = "#EDD4D6")
  EWS.colors <- c("dark3" = "#141F3B",
                  "dark2" = "#273D76",
                  "dark1" = "#3B5CB0",
                  "main" = "#4E7AEB",
                  "light1" = "#7195EF",
                  "light2" = "#95AFF3",
                  "light3" = "#DCE4FB")
  NBL.colors <- c("dark3" = "#27132B",
                  "dark2" = "#4F2757",
                  "dark1" = "#763A82",
                  "main" = "#9D4DAD",
                  "light1" = "#B171BD",
                  "light2" = "#C494CE",
                  "light3" = "#EBDBEF")
  OS.colors <- c("dark3" = "#072607",
                  "dark2" = "#0E4B0E",
                  "dark1" = "#147114",
                  "main" = "#1B961B",
                  "light1" = "#49AB49",
                  "light2" = "#76C076",
                  "light3" = "#D1EAD1")
  control.colors <- c("dark3" = "#363636",
                      "dark2" = "#6B6B6B",
                      "dark1" = "#A1A1A1",
                      "main" = "#D6D6D6",
                      "light1" = "#DEDEDE",
                      "light2" = "#E6E6E6",
                      "light3" = "#F7F7F7")
  male.colors <- c("dark2" = "#2A5869",
                   "dark1" = "#3F839D",
                   "main" = "#54AFD1",
                   "light1" = "#76BFDA",
                   "light2" = "#98CFE3")
  female.colors <- c("dark2" = "#73475E",
                     "dark1" = "#AD6B8C",
                     "main" = "#E68EBB",
                     "light1" = "#EBA5C9",
                     "light2" = "#F0BBD6")
  AFR.colors <- c("dark2" = "#796B28",
                  "dark1" = "#B5A13C",
                  "main" = "#F1D650",
                  "light1" = "#F4DE73",
                  "light2" = "#F7E696")
  AMR.colors <- c("dark2" = "#543726",
                  "dark1" = "#7E5239",
                  "main" = "#A86D4C",
                  "light1" = "#B98A70",
                  "light2" = "#CBA794")
  EAS.colors <- c("dark2" = "#5B651F",
                  "dark1" = "#89982E",
                  "main" = "#B6CA3D",
                  "light1" = "#C5D564",
                  "light2" = "#D3DF8B")
  EUR.colors <- c("dark2" = "#4C626B",
                  "dark1" = "#7292A1",
                  "main" = "#98C3D6",
                  "light1" = "#ADCFDE",
                  "light2" = "#C1DBE6")
  SAS.colors <- c("dark2" = "#6D3C3A",
                  "dark1" = "#A35A56",
                  "main" = "#D97873",
                  "light1" = "#E1938F",
                  "light2" = "#E8AEAB")
  DEL.colors <- c("dark2" = "#6A1D13",
                  "dark1" = "#9F2B1C",
                  "main" = "#D43925",
                  "light1" = "#DD6151",
                  "light2" = "#EEB0A8")
  DUP.colors <- c("dark2" = "#123B59",
                  "dark1" = "#1A5985",
                  "main" = "#2376B2",
                  "light1" = "#4F91C1",
                  "light2" = "#A7C8E0")
  CNV.colors <- c("dark2" = "#3A2D59",
                  "dark1" = "#574385",
                  "main" = "#7459B2",
                  "light1" = "#907AC1",
                  "light2" = "#C7BDE0")
  INS.colors <- c("dark2" = "#6A3A70",
                  "dark1" = "#9F57A8",
                  "main" = "#D474E0",
                  "light1" = "#DD90E6",
                  "light2" = "#EEC7F3")
  INV.colors <- c("dark2" = "#7D4A0F",
                  "dark1" = "#BB6E16",
                  "main" = "#FA931E",
                  "light1" = "#FBA94B",
                  "light2" = "#FDD4A5")
  CPX.colors <- c("dark2" = "#397246",
                  "dark1" = "#55AA69",
                  "main" = "#71E38C",
                  "light1" = "#8DE9A3",
                  "light2" = "#C6F4D1")
  CTX.colors <- c("dark2" = "#1D3923",
                  "dark1" = "#2B5635",
                  "main" = "#397246",
                  "light1" = "#618E6B",
                  "light2" = "#B0C7B5")
  OTH.colors <- BND.colors <- CTX.colors
  sv.colors <- c("DEL" = DEL.colors[["main"]],
                 "DUP" = DUP.colors[["main"]],
                 "CNV" = CNV.colors[["main"]],
                 "INS" = INS.colors[["main"]],
                 "INV" = INV.colors[["main"]],
                 "CPX" = CPX.colors[["main"]],
                 "BND" = BND.colors[["main"]],
                 "CTX" = CTX.colors[["main"]],
                 "OTH" = OTH.colors[["main"]])
  colors <- list(
    "cancer.colors" = c("pancan" = pancan.colors[["main"]],
                        "EWS" = EWS.colors[["main"]],
                        "NBL" = NBL.colors[["main"]],
                        "OS" = OS.colors[["main"]],
                        "control" = control.colors[["main"]]),
    "cancer.palettes" = list("pancan" = pancan.colors,
                             "EWS" = EWS.colors,
                             "NBL" = NBL.colors,
                             "OS" = OS.colors,
                             "control" = control.colors),
    "control.colors" = control.colors,
    "EWS.colors" = EWS.colors,
    "NBL.colors" = NBL.colors,
    "OS.colors" = OS.colors,
    "sex.colors" = c("MALE" = male.colors[["main"]],
                     "FEMALE" = female.colors[["main"]],
                     "OTHER" = "#666245"),
    "MALE.colors" = male.colors,
    "FEMALE.colors" = female.colors,
    "pop.colors" = c("AFR" = AFR.colors[["main"]],
                     "AMR" = AMR.colors[["main"]],
                     "EAS" = EAS.colors[["main"]],
                     "EUR" = EUR.colors[["main"]],
                     "SAS" = SAS.colors[["main"]],
                     "OTH" = "gray"),
    "pop.palettes" = list("AFR" = AFR.colors,
                          "AMR" = AMR.colors,
                          "EAS" = EAS.colors,
                          "EUR" = EUR.colors,
                          "SAS" = SAS.colors,
                          "OTH" = "gray"),
    "AFR.colors" = AFR.colors,
    "AMR.colors" = AMR.colors,
    "EAS.colors" = EAS.colors,
    "EUR.colors" = EUR.colors,
    "SAS.colors" = SAS.colors,
    "DEL.colors" = DEL.colors,
    "DUP.colors" = DUP.colors,
    "CNV.colors" = CNV.colors,
    "INS.colors" = INS.colors,
    "INV.colors" = INV.colors,
    "CPX.colors" = CPX.colors,
    "BND.colors" = BND.colors,
    "CTX.colors" = CTX.colors,
    "OTH.colors" = OTH.colors,
    "sv.colors" = sv.colors,
    "sv.palettes" = list("DEL" = DEL.colors,
                         "DUP" = DUP.colors,
                         "CNV" = CNV.colors,
                         "INS" = INS.colors,
                         "INV" = INV.colors,
                         "CPX" = CPX.colors,
                         "BND" = BND.colors,
                         "CTX" = CTX.colors,
                         "OTH" = OTH.colors),
    "stage.colors" = c("0" = "white",
                       "1" = "#F8FAA7",
                       "2" = "#FFCC66",
                       "3" = "#FE8002",
                       "4" = "#ED3823"),
    "relative.colors" = c("duplicates" = "#9D1309",
                          "parent-child" = "#FF6103",
                          "siblings" = "#FFB14D",
                          "unrelated" = "#AAAAAA"))

  # Define scales
  logscale.major <- 10^(-10:10)
  contig.lengths <- c("chr1" = 248956422,
                      "chr2" = 242193529,
                      "chr3" = 198295559,
                      "chr4" = 190214555,
                      "chr5" = 181538259,
                      "chr6" = 170805979,
                      "chr7" = 159345973,
                      "chr8" = 145138636,
                      "chr9" = 138394717,
                      "chr10" = 133797422,
                      "chr11" = 135086622,
                      "chr12" = 133275309,
                      "chr13" = 114364328,
                      "chr14" = 107043718,
                      "chr15" = 101991189,
                      "chr16" = 90338345,
                      "chr17" = 83257441,
                      "chr18" = 80373285,
                      "chr19" = 58617616,
                      "chr20" = 64444167,
                      "chr21" = 46709983,
                      "chr22" = 50818468,
                      "chrX" = 156040895,
                      "chrY" = 57227415)
  scales <- list(
    "logscale.major" = logscale.major,
    "logscale.major.bp" = 10^(0:9),
    "logscale.major.bp.labels" = c(sapply(c("bp", "kb", "Mb"),
                                          function(suf){paste(c(1, 10, 100), suf, sep="")}),
                                   "1 Gb"),
    "logscale.demi" = as.numeric(sapply(logscale.major, function(e){c(1, 5)*e})),
    "logscale.demi.bp" = as.numeric(sapply(10^(0:9), function(e){c(1, 5)*e})),
    "logscale.demi.bp.labels" = c(paste(c(1, 5, 10, 50, 100, 500), "bp", sep=""),
                                  paste(c(1, 5, 10, 50, 100, 500), "kb", sep=""),
                                  paste(c(1, 5, 10, 50, 100, 500), "Mb", sep=""),
                                  paste(c(1, 5), "Gb", sep="")),
    "logscale.minor" = as.numeric(sapply(logscale.major, function(e){(1:9)*e})),
    "yearscale.major" = 0:100 * 365,
    "yearscale.demi" = seq(0, 100, 0.5) * 365,
    "yearscale.minor" = seq(0, 100, 1/12) * 365,
    "contig.lengths" = contig.lengths
  )

  # Define names
  all.names <- list(
    "cancer.names.short" = c("pancan" = "All cases",
                             "EWS" = "Ewing",
                             "NBL" = "Neuroblast.",
                             "OS" = "Osteosarc.",
                             "control" = "Control"),
    "cancer.names.long" = c("pancan" = "All cases",
                            "EWS" = "Ewing sarcoma",
                            "NBL" = "Neuroblastoma",
                            "OS" = "Osteosarcoma",
                            "control" = "Controls"),
    "cancer.names.vshort" = c("pancan" = "All cases",
                             "EWS" = "Ewing",
                             "NBL" = "Neuro.",
                             "OS" = "Osteo.",
                             "control" = "Controls"),
    "metadata.cancer.label.map" = c("case" = "pancan",
                                    "ewing" = "EWS",
                                    "neuroblastoma" = "NBL",
                                    "osteosarcoma" = "OS",
                                    "control" = "control"),
    "pop.abbreviations" = c("AFR" = "Afr.",
                            "AMR" = "Amer.",
                            "EAS" = "E. Asn.",
                            "EUR" = "Eur.",
                            "SAS" = "S. Asn.",
                            "OTH" = "Other"),
    "pop.names.short" = c("AFR" = "African",
                          "AMR" = "American",
                          "EAS" = "E. Asian",
                          "EUR" = "European",
                          "SAS" = "S. Asian",
                          "OTH" = "Other"),
    "pop.names.long" = c("AFR" = "African/African-American",
                         "AMR" = "Latino/admixed American",
                         "EAS" = "East Asian",
                         "EUR" = "European",
                         "SAS" = "South Asian",
                         "OTH" = "Other/unknown"),
    "sex.names" = c("MALE" = "Male",
                    "FEMALE" = "Female",
                    "OTHER" = "Other"),
    "stage.names" = c("0" = "",
                      "1" = "I",
                      "2" = "II",
                      "3" = "III",
                      "4" = "IV"),
    "relative.names" = c("duplicates" = "Identical",
                         "parent-child" = "Parent-child",
                         "siblings" = "Siblings",
                         "unrelated" = "Unrelated"),
    "sv.abbreviations" = c("DEL" = "Del.",
                           "DUP" = "Dup.",
                           "CNV" = "mCNV",
                           "INS" = "Ins.",
                           "INV" = "Inv.",
                           "CPX" = "Complex",
                           "CTX" = "Tloc.",
                           "OTH" = "Other"),
    "sv.names" = c("DEL" = "Deletion",
                   "DUP" = "Duplication",
                   "CNV" = "mCNV",
                   "INS" = "Insertion",
                   "INV" = "Inversion",
                   "CPX" = "Complex SV",
                   "CTX" = "Translocation",
                   "OTH" = "Other SV"),
    "freq.names" = c("rare" = "Rare",
                     "vrare" = "V. rare",
                     "singleton" = "Singleton"),
    "cohort.names" = c("GMKF" = "GMKF",
                       "ICGC" = "ICGC",
                       "StJude" = "St. Jude",
                       "Topmed_BIOME" = "BioMe",
                       "Topmed_MESA" = "MESA")
  )

  # Assign constants to global environment
  if(length(intersect(subset, c("colors", "all"))) > 0){
    for(variable in names(colors)){
      assign(variable, colors[[variable]], envir=envir)
    }
  }
  if(length(intersect(subset, c("scales", "all"))) > 0){
    for(variable in names(scales)){
      assign(variable, scales[[variable]], envir=envir)
    }
  }
  if(length(intersect(subset, c("names", "all"))) > 0){
    for(variable in names(all.names)){
      assign(variable, all.names[[variable]], envir=envir)
    }
  }
}
