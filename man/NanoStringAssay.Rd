\name{NanoStringAssay}
\alias{NanoStringAssay}
\title{NanoString Assay Constructor}
\usage{
  NanoStringAssay(rccfiles = NULL, keyfile = NULL, idvars,
    primerid, measurement, ncells = NULL, geneid = NULL,
    id = NULL, cellvars = NULL, featurevars = NULL,
    phenovars = NULL, post.process.function = NULL,
    ltrans = FALSE, ...)
}
\arguments{
  \item{rccfiles}{A list of rcc files with full paths}

  \item{keyfile}{A full path to a keyfile}

  \item{idvars}{See \code{\link{SingleCellAssay}}}

  \item{primerid}{See \code{\link{SingleCellAssay}}}

  \item{measurement}{The measurement column for raw data.
  Will be placed in a variable named \code{raw} in the
  environment storing the data. The transformed and
  thresholded data will be placed in data.}

  \item{ncells}{A \code{character} specifying the column
  which gives the number of cells per well}

  \item{geneid}{See \code{\link{SingleCellAssay}}}

  \item{id}{An identifier for the resulting object. Should
  be a meaningful name}

  \item{cellvars}{See \code{\link{SingleCellAssay}}}

  \item{featurevars}{See \code{\link{SingleCellAssay}}}

  \item{phenovars}{See \code{\link{SingleCellAssay}}}

  \item{post.process.function}{function applied to
  \code{data.frame} of all rcc files, before the
  NanostringAssay object is constructed.}

  \item{ltrans}{Should the counts be log2 + 1 transformed?}

  \item{...}{Additional parameters passed to
  \code{SingleCellAssay} constructor}
}
\value{
  A FluidigmAssay object
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

