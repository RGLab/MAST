\docType{methods}
\name{mappingIntersection-methods}
\alias{mappingIntersection}
\alias{mappingIntersection,Mapping,character,character-method}
\alias{mappingIntersection,Mapping,missing,missing-method}
\alias{mappingIntersection-methods}
\title{mappingIntersection}
\arguments{
  \item{object}{A \code{Mapping} object}

  \item{name1}{An optional \code{character} map name}

  \item{name2}{Another optional \code{character} map name
  to be compared against \code{name1}}
}
\value{
  A matrix of intersecitons between different mappings if
  \code{name1} and \code{name2} are omitted, or a numeric
  value otherwise.
}
\description{
  mappingIntersection checks the intersection of two maps
  in a \code{Mapping} object.
}

