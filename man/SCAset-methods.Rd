\name{SCASet}
\alias{SCASet}
\title{SCASet constructor}
\usage{
  SCASet(dataframe, splitby, idvars = NULL,
    primerid = NULL, measurement = NULL, ...)
}
\arguments{
  \item{dataframe}{flat data.frame ala SingleCellAssay}

  \item{splitby}{either a character vector naming columns
  or a factor or a list of factors used to split dataframe
  into SingleCellAssays}

  \item{idvars}{character vector naming columns that
  uniquely identify a cell}

  \item{primerid}{character vector of length 1 that names
  the column that containing what feature was measured}

  \item{measurement}{character vector of length 1 that
  names the column containing the measurement}

  \item{...}{passed up to SingleCellAssay}
}
\value{
  SCASet
}
\description{
  Constructs a SCASet, which is a list of SingleCellAssays
}
\details{
  FIXME
}

