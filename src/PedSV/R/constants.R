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
  pancan.colors <- c("dark3" = "#2A120C",
                     "dark2" = "#542319",
                     "dark1" = "#7E3525",
                     "main" = "#A84631",
                     "light1" = "#B96B5A",
                     "light2" = "#CB9083",
                     "light3" = "#EEDAD6")
  OS.colors <- c("dark3" = "#372C10",
                 "dark2" = "#6E5721",
                 "dark1" = "#A48332",
                 "main" = "#DBAE42",
                 "light1" = "#E2BE68",
                 "light2" = "#E9CE8E",
                 "light3" = "#F8EFD9")
  NBL.colors <- c("dark3" = "#0D260D",
                  "dark2" = "#1B4B1B",
                  "dark1" = "#287128",
                  "main" = "#359635",
                  "light1" = "#5DAB5D",
                  "light2" = "#86C086",
                  "light3" = "#D7EAD7")
  EWS.colors <- c("dark3" = "#16233B",
                  "dark2" = "#2D4676",
                  "dark1" = "#4368B0",
                  "main" = "#598BEB",
                  "light1" = "#7AA2EF",
                  "light2" = "#9BB9F3",
                  "light3" = "#DEE8FB")
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
  female.colors <- c("dark2" = "#692A4A",
                     "dark1" = "#9D3F6F",
                     "main" = "#D15494",
                     "light1" = "#DA76A9",
                     "light2" = "#E398BF")
  AFR.colors <- c("dark2" = "#796624",
                  "dark1" = "#B59836",
                  "main" = "#F1CB48",
                  "light1" = "#F4D56D",
                  "light2" = "#F7E091")
  AMR.colors <- c("dark2" = "#652223",
                  "dark1" = "#973234",
                  "main" = "#C94345",
                  "light1" = "#D4696A",
                  "light2" = "#DF8E8F")
  EAS.colors <- c("dark2" = "#5B6519",
                  "dark1" = "#889826",
                  "main" = "#B5CA32",
                  "light1" = "#C4D55B",
                  "light2" = "#D3DF84")
  EUR.colors <- c("dark2" = "#4E6774",
                  "dark1" = "#749AAE",
                  "main" = "#9BCDE8",
                  "light1" = "#AFD7ED",
                  "light2" = "#C3E1F1")
  SAS.colors <- c("dark2" = "#4D2D4E",
                  "dark1" = "#734474",
                  "main" = "#995A9B",
                  "light1" = "#AD7BAF",
                  "light2" = "#C29CC3")
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
                        "OS" = OS.colors[["main"]],
                        "NBL" = NBL.colors[["main"]],
                        "EWS" = EWS.colors[["main"]],
                        "control" = control.colors[["main"]]),
    "cancer.palettes" = list("pancan" = pancan.colors,
                             "OS" = OS.colors,
                             "NBL" = NBL.colors,
                             "EWS" = EWS.colors,
                             "control" = control.colors),
    "control.colors" = control.colors,
    "OS.colors" = OS.colors,
    "NBL.colors" = NBL.colors,
    "EWS.colors" = EWS.colors,
    "sex.colors" = c("MALE" = male.colors[["main"]],
                     "FEMALE" = female.colors[["main"]]),
    "MALE.colors" = male.colors,
    "FEMALE.colors" = female.colors,
    "pop.colors" = c("AFR" = AFR.colors[["main"]],
                     "AMR" = AMR.colors[["main"]],
                     "EAS" = EAS.colors[["main"]],
                     "EUR" = EUR.colors[["main"]],
                     "SAS" = SAS.colors[["main"]],
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
    "yearscale.minor" = seq(0, 100, 1/12) * 365
  )

  # Define names
  all.names <- list(
    "cancer.names.short" = c("pancan" = "Pan-Can.",
                             "OS" = "Osteosarc.",
                             "NBL" = "Neuroblast.",
                             "EWS" = "Ewing",
                             "control" = "Control"),
    "cancer.names.long" = c("pancan" = "Pan-Cancer",
                            "OS" = "Osteosarcoma",
                            "NBL" = "Neuroblastoma",
                            "EWS" = "Ewing Sarcoma",
                            "control" = "Control"),
    "metadata.cancer.label.map" = c("osteosarcoma" = "OS",
                                    "neuroblastoma" = "NBL",
                                    "ewing" = "EWS",
                                    "case" = "pancan",
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
                         "AMR" = "Latino/Admixed American",
                         "EAS" = "East Asian",
                         "EUR" = "European",
                         "SAS" = "South Asian",
                         "OTH" = "Other/Unknown"),
    "stage.names" = c("0" = "",
                      "1" = "I",
                      "2" = "II",
                      "3" = "III",
                      "4" = "IV"),
    "relative.names" = c("duplicates" = "Identical",
                         "parent-child" = "Parent-Child",
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
                     "vrare" = "V. Rare",
                     "singleton" = "Singleton")
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
