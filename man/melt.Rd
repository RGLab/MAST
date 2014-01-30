\docType{methods}
\name{melt.SingleCellAssay}
\alias{melt}
\alias{melt,SingleCellAssay-method}
\alias{melt.SingleCellAssay}
\title{melt}
\usage{
melt.SingleCellAssay(data, ...)

\S4method{melt}{SingleCellAssay}(data, ...)
}
\arguments{
  \item{data}{A rectangular array, with attributes attached
  to its rows and columns}

  \item{...}{ignored}
}
\value{
A \code{data.frame} typically, with the cartesian product
of the row and column attributes and the values from the
rectangular array
}
\description{
Melt a rectangular array
}
\details{
Return a 'molten' (flat) representation of a rectangular
array

\code{signature(data="SingleCellAssay")}: return a
\code{data.frame}, which contains a melted version of
\code{data}.
}
\keyword{transformation}

