\docType{methods}
\name{LRT}
\alias{LRT}
\alias{LRT,SingleCellAssay,character,character,character-method}
\alias{LRT,SingleCellAssay,character,character,missing-method}
\alias{LRT,SingleCellAssay,character,character,NULL-method}
\title{Likelihood Ratio Tests for SingleCellAssays}
\arguments{
  \item{sca}{A \code{SingleCellAssay} class object}

  \item{comparison}{A \code{character} specifying the
  factor for comparison}

  \item{referent}{A \code{character} specifying the
  reference level of \code{comparison}.}

  \item{groups}{A optional \code{character} specifying a
  variable on which to stratify the test.  For each level
  of \code{groups}, there will be a separate likelihood
  ratio test.}

  \item{returnall}{A \code{logical} specifying if
  additional columns should be returned with information
  about the different componnetns of the test.}
}
\description{
  Tests for a change in ET binomial proportion or mean of
  positive ET Likelihood Ratio Test for SingleCellAssay
  objects
}
\details{
  Combined Likelihood ratio test (binomial and normal) for
  SingleCellAssay and derived objects
}

