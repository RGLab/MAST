\name{FluidigmAssay}
\alias{FluidigmAssay}
\title{Fluidigm Assay Constructor}
\usage{
  FluidigmAssay(dataframe = NULL, idvars, primerid,
    measurement, ncells, geneid = NULL, id = numeric(0),
    cellvars = NULL, featurevars = NULL, phenovars = NULL,
    ...)
}
\arguments{
  \item{dataframe}{A data.table containing the raw data}

  \item{idvars}{See \code{\link{SingleCellAssay}}}

  \item{primerid}{See \code{\link{SingleCellAssay}}}

  \item{measurement}{See \code{\link{SingleCellAssay}}}

  \item{ncells}{A \code{character} specifying the column
  which gives the number of cells per well}

  \item{geneid}{See \code{\link{SingleCellAssay}}}

  \item{id}{An identifier for the resulting object. Should
  be a meaningful name.}

  \item{cellvars}{See \code{\link{SingleCellAssay}}}

  \item{featurevars}{See \code{\link{SingleCellAssay}}}

  \item{phenovars}{See \code{\link{SingleCellAssay}}}

  \item{...}{Additional parameters passed to
  \code{SingleCellAssay} constructor}
}
\value{
  A FluidigmAssay object
}
\description{
  Constructor for a FluidigmAssay
}
\details{
  Constructs a FluidigmAssay object. Differs little from
  the SingleCellAssay constructor. Only the \code{ncells}
  parameter is additionally required. Mapping argument has
  been removed for simplicity.
}
\author{
  Andrew McDavid and Greg Finak
}

