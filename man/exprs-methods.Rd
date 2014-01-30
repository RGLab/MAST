\docType{methods}
\name{exprs}
\alias{exprs}
\alias{exprs,DataLayer-method}
\alias{exprs<-,DataLayer,ANY-method}
\title{exprs}
\usage{
\S4method{exprs}{DataLayer}(object)
}
\arguments{
  \item{object}{DataLayer}
}
\value{
numeric matrix

a \code{matrix} of measurement values with wells on the
rows and features on the columns of the default layer
}
\description{
Get or set a matrix of measurement values in a
\code{SingleCellAssay}
}
\details{
Return or set a matrix of the measurement: cells by
primerids
}

