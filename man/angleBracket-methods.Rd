\name{[}
\alias{[}
\alias{[,DataLayer,ANY-method}
\alias{[,SingleCellAssay,ANY-method}
\title{subset methods}
\description{
  subset methods
}
\details{
  \code{signature(x="SingleCellAssay", i="ANY")}:
  \code{x[i]}, where \code{i} is a logical, integer, or
  character vector, recycled as necessary to match
  \code{nrow(x)}. Optional \code{x[[i,j]]} where j is a
  logical, integer or character vector selecting the
  features based on ``primerid'' which is unique, while
  ``geneid'' or gene name is not necessarily unique.
}
\keyword{transform}

