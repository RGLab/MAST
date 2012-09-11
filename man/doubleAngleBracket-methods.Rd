\docType{methods}
\name{[[}
\alias{[[}
\alias{[[,SCASet,ANY-method}
\alias{[[,SingleCellAssay,ANY-method}
\title{Index into an SCASet object}
\arguments{
  \item{x}{object to be subscripted}

  \item{i}{index}

  \item{...}{Ignored}
}
\value{
  subscripted SingleCellAssay or derived class
}
\description{
  Return a \code{SingleCellAssay} at index/indices i
}
\details{
  \code{signature(x="SCASet", i="ANY")}: \code{x[[i]]},
  where \code{i} is length-1 integer or character matching
  sampleNames.  Returns the SingleCellAssay at position or
  with sampleName \code{i}.

  \code{signature(x="SingleCellAssay", i="ANY")}:
  \code{x[[i]]}, where \code{i} is a logical, integer, or
  character vector, recycled as necessary to match
  \code{nrow(x)}. Optional \code{x[[i,j]]} where j is a
  logical, integer or character vector selecting the
  features based on ``primerid'' which is unique, while
  ``geneid'' or gene name is not necessarily unique.
}
\keyword{transform}

