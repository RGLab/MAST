\docType{methods}
\name{split}
\alias{split}
\alias{split,DataLayer,ANY-method}
\alias{split,SingleCellAssay,ANY-method}
\title{Split into SCASet}
\arguments{
  \item{x}{SingleCellAssay}

  \item{f}{length-1 character or factor of length nrow(x)}
}
\value{
  SCASet
}
\description{
  Splits a \code{SingleCellAssay} into a \code{SCASet} by a
  factor (or something coercible into a factor) or a
  character giving a column of the melted SingleCellAssay
}

