\docType{methods}
\name{Extract.SCASet}
\alias{Extract.SCASet}
\alias{[[,SCASet,ANY-method}
\title{Index into an SCASet object}
\usage{
\S4method{[[}{SCASet}(x, i, j, ...)
}
\arguments{
  \item{x}{object to be subscripted}

  \item{i}{index}

  \item{j}{Ignored}

  \item{...}{Ignored}
}
\value{
subscripted SingleCellAssay or derived class
}
\description{
Return a \code{SingleCellAssay} at index i
}
\details{
\code{signature(x="SCASet", i="ANY")}: \code{x[[i]]}, where
\code{i} is length-1 integer or character matching
sampleNames.  Returns the SingleCellAssay at position or
with sampleName \code{i}.
}
\keyword{transform}

