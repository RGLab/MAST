\name{filter}
\alias{filter}
\title{Filter a SingleCellAssay or Fluidigm Assay}
\usage{
  filter(sc, groups = NULL, filt_control = NULL,
    apply_filter = TRUE)
}
\arguments{
  \item{sc}{The \code{SingleCellAssay} object}

  \item{groups}{The \code{character} naming the grouping
  variable (optional)}

  \item{filt_control}{The \code{list} with configuration
  parameters for the filter.}

  \item{apply_filter}{\code{logical} should the filter be
  applied?}
}
\value{
  A filtered result
}
\description{
  Function that filters a single cell assay
}
\details{
  The function filters wells that don't pass filtering
  criteria described in filt_control. filt_control is a
  list with named elements nOutliers (minimum nmber of
  outlier cells for a cell to be filtered. sigmaContinuous
  (the z-score outlier threshold for the continuous part of
  the signal), and sigmaProportion (the z-score outlier
  threshold for the discrete part of the signal).
}
\author{
  Andrew McDavid
}

