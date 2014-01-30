\docType{methods}
\name{split,DataLayer-method}
\alias{split,DataLayer-method}
\title{Split a DataLayer correctly}
\usage{
\S4method{split}{DataLayer}(x, f, drop = FALSE, ...)
}
\arguments{
  \item{x}{vector or data frame containing values to be
  divided into groups.}

  \item{f}{a \sQuote{factor} in the sense that
  \code{\link{as.factor}(f)} defines the grouping, or a
  list of such factors in which case their interaction is
  used for the grouping.}

  \item{drop}{logical indicating if levels that do not
  occur should be dropped (if \code{f} is a \code{factor}
  or a list).}

  \item{...}{further potential arguments passed to
  methods.}
}
\description{
split.default doesn't seem to do the right thing for an
array
}

