\docType{methods}
\name{split}
\alias{split}
\alias{split,SingleCellAssay,ANY-method}
\title{Split into SCASet}
\usage{
\S4method{split}{SingleCellAssay}(x, f, drop = FALSE, ...)
}
\arguments{
  \item{x}{SingleCellAssay}

  \item{f}{length-1 character or factor of length nrow(x)}

  \item{drop}{logical indicating if levels that do not
  occur should be dropped (if \code{f} is a \code{factor}
  or a list).}

  \item{...}{further potential arguments passed to
  methods.}
}
\value{
SCASet
}
\description{
Splits a \code{SingleCellAssay} into a \code{SCASet} by a
factor (or something coercible into a factor) or a
character giving a column of the melted SingleCellAssay
}

