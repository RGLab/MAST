\docType{methods}
\name{subset}
\alias{subset}
\alias{subset,SingleCellAssay-method}
\title{Evaluate an expression (returning a \code{logical} vector in the context of
\code{x}}
\arguments{
  \item{x}{A matrix or \code{data.frame} typically}

  \item{thesubset}{An expression, which evaluated in either
  the current environment, or the parent, yields a logical
  vector}
}
\value{
  The same class as \code{x} typically, with entries in
  which thesubset was \code{TRUE}
}
\description{
  Really just syntactic sugar.
}
\details{
  \code{signature(x='SingleCellAssay', thesubset='ANY')}:
  Return a new SingleCellAssay consisting of cells in which
  thesubset is TRUE
}
\keyword{transformation}

