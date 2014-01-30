\docType{methods}
\name{copy}
\alias{copy}
\alias{copy,SingleCellAssay-method}
\title{Return a deep copy of an object}
\usage{
copy(object)

\S4method{copy}{SingleCellAssay}(object)
}
\arguments{
  \item{sc}{object}
}
\value{
copy of sc
}
\description{
Because SingleCellAssay keeps the datastore in an
environment, it's not sufficient to create a new object to
make a copy on an object
}

