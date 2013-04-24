\name{NanoStringAssay}
\alias{NanoStringAssay}
\title{NanoString Assay Constructor}
\usage{
  NanoStringAssay(rccfiles = NULL, keyfile = NULL, idvars,
    primerid, measurement = "Count", ncells = NULL,
    id = NULL, cellvars = NULL, featurevars = NULL,
    phenovars = NULL, post.process.function = NULL, ...)
}
\arguments{
  \item{rccfiles}{A list of rcc files with full paths}

  \item{keyfile}{A full path to a keyfile}

  \item{idvars}{See \code{\link{SingleCellAssay}}}

  \item{primerid}{See \code{\link{SingleCellAssay}}}

  \item{measurement}{The measurement column for raw data.}

  \item{ncells}{A \code{character} specifying the column
  which gives the number of cells per well}

  \item{id}{An identifier for the resulting object. Should
  be a meaningful name}

  \item{cellvars}{See \code{\link{SingleCellAssay}}}

  \item{featurevars}{See \code{\link{SingleCellAssay}}}

  \item{phenovars}{See \code{\link{SingleCellAssay}}}

  \item{post.process.function}{function applied to
  \code{data.frame} of all rcc files, before the
  NanostringAssay object is constructed.}

  \item{...}{Additional parameters passed to
  \code{SingleCellAssay} constructor}
}
\value{
  A NanoStringAssay object
}
\description{
  Constructor for a NanoStringAssay
}
\details{
  Constructs a NanoStringAssay object. Differs little from
  the FluidigmAssay constructor. Accepts a function for
  post-processing. mapping argument has been removed as
  well for simplicity.
}
\author{
  Andrew McDavid and Greg Finak
}

