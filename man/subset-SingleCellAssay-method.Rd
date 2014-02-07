\docType{methods}
\name{subset,SingleCellAssay-method}
\alias{subset,SingleCellAssay-method}
\title{Subset a SingleCellAssay by data in cellData}
\usage{
\S4method{subset}{SingleCellAssay}(x, thesubset, ...)
}
\arguments{
  \item{x}{SingleCellAssay}

  \item{thesubset}{expression, which when evaluated in
  cellData environment which returns a logical}

  \item{...}{further arguments to be passed to or from
  other methods.}
}
\description{
Subset a SingleCellAssay by data in cellData
}
\examples{
data(vbetaFA)
subset(vbetaFA, ncells==1)
}

