\docType{methods}
\name{getMapping}
\alias{getMapping}
\alias{getMapping,Mapping,character-method}
\alias{getMapping,Mapping,missing-method}
\alias{getMapping,SCA,character-method}
\alias{getMapping,SingleCellAssay,missing-method}
\title{getMapping}
\arguments{
  \item{object}{A \code{Mapping} object or a
  \code{SingleCellAssay} or object inheriting from
  \code{SingleCellAssay}}

  \item{mapnames}{A \code{character} vector of map names.
  Can be empty}
}
\value{
  A \code{Mapping} object if \code{object} is a
  \code{SingleCellAssay}. A character vector of mapped
  columns if \code{mapnames} is provided.
}
\description{
  Accessor for mapping
}
\details{
  This returns the mapping, which is a named list.  The
  names of the list are keywords, while the contents of the
  list are column names in the melted
  \code{SingleCellAssay}.  The \code{mapping} identifies
  columns in the melted which trigger special object
  behavior. Names recognized by the class are contained in
  \link[=getMapNames]{mapNames}, but additional names may
  be defined.
}
\keyword{accessor}

