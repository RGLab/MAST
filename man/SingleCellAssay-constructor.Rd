\docType{methods}
\name{SingleCellAssay-constructor}
\alias{SingleCellAssay}
\alias{SingleCellAssay-constructor}
\title{SingleCellAssay: A constructor for an object of type SingleCellAssay.}
\usage{
  SingleCellAssay(dataframe = NULL, idvars = NULL,
    primerid = NULL, measurement = NULL, geneid = NULL,
    id = numeric(0), mapping = NULL, cellvars = NULL,
    featurevars = NULL, phenovars = NULL, ...)
}
\arguments{
  \item{dataframe}{A 'flattened' data.table containing
  columns giving cell and feature identifiers and a
  measurement column}

  \item{idvars}{character vector naming columns that
  uniquely identify a cell}

  \item{primerid}{character vector of length 1 that names
  the column that identifies what feature (i.e. gene) was
  measured}

  \item{measurement}{character vector of length 1 that
  names the column containing the measurement}

  \item{geneid}{character vector of length 1 that names a
  'gene' column.  This could be placed into ...?}

  \item{id}{optional numeric (not sure what this is
  supposed to do)}

  \item{mapping}{named list.  Names are identifiers used by
  find special columns in the dataframe.  This shouldn't be
  modified directly by the user}

  \item{cellvars}{Character vector naming columns
  containing additional cellular metadata}

  \item{featurevars}{Character vector naming columns
  containing additional feature metadata}

  \item{phenovars}{Character vector naming columns
  containing additional phenotype metadata}

  \item{...}{Additional keywords to be added to mapping}
}
\value{
  SingleCellAssay object
}
\description{
  This is the constructor for the class. This class intends
  to ease the analysis of single cell assays, in which
  multiple, exchangible, cells from an experimental unit
  (patient, or organism) are assayed along several (or
  many) dimensions, such as genes. A few examples of this
  might be Fluidigm gene expression chips, or single cell
  sequencing experiments.  The chief functionality is to
  make it easy to keep cellular-level metadata linked to
  the measurements through \code{cellData} and
  \code{phenoData}.  There are also subsetting and
  splitting measures to coerce between a SingleCellAssay,
  and a \link{SCASet}.
}

