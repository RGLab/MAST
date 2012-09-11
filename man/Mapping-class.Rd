\docType{class}
\name{Mapping-class}
\alias{Mapping}
\alias{Mapping-class}
\title{Column Mapping for SingleCellAssays}
\description{
  Mapping class for SingleCellAssay package
}
\details{
  A class that represents a mapping of columns in a raw
  data file to cell-leve, feature-level, and
  phenotype-level metadata, as well as unique identifiers
  for individual cells. mapNames for the SingleCellAssay
  class are in the object
  \code{SingleCellAssay:::SingleCellAssayMapNames} mapNames
  for the FluidigmAssay class are in the object
  \code{SingleCellAssay:::FluidigmMapNames} }
  \section{Slots}{ \describe{ \item{mapping}{A named list
  providing the mapping from required fields in a
  SingleCellAssay or FluidigmAssay to column names in the
  data file} }
}
\seealso{
  \code{\link{SingleCellAssay}},\code{\link{FluidigmAssay}},
  \code{\link{getMapping}},\code{\link{addMapping}},\code{\link{getMapNames}},\code{\link{isEmpty}},\code{\link{mappingIntersection}}
}

