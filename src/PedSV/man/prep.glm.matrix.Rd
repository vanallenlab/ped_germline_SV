% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assoc.R
\name{prep.glm.matrix}
\alias{prep.glm.matrix}
\title{Construct association dataframe}
\usage{
prep.glm.matrix(meta, X, Y, use.N.pcs = 3, extra.terms = NULL)
}
\arguments{
\item{meta}{Sample metadata loaded with \link{load.sample.metadata}}

\item{X}{Vector of values for primary independent variable. See \code{Details}.}

\item{Y}{Vector of values for dependent variable. See \code{Details}.}

\item{use.N.pcs}{Specify how many principal components should be adjusted in
model [default: 3]}

\item{extra.terms}{Specify if any extra terms should be added to the model.
Named options include:  "batch", "coverage", "insert.size", and "wgd".
Custom terms can be passed using their exact column names in \code{meta}.}
}
\value{
data.frame
}
\description{
Prepare a dataframe with dependent and independent variables for association
testing with additional covariates as optioned
}
\details{
There are several options for providing \code{X} and \code{Y} values:
\itemize{
\item As an unnamed vector. In this case, the values are assumed to be in the
same order as the samples in \code{meta}.
\item As a named vector. In this case, the vector names are assumed to be sample
IDs, and any sample ID failing to match in \code{meta} will be dropped.
\item As a data frame. In this case, the row names are assumed to be sample IDs,
and any sample ID failing to match in \code{meta} will be dropped.
}
}
\seealso{
\link{pedsv.glm}
}
