\docType{methods}
\name{combine,DataLayer,DataLayer-method}
\alias{combine,DataLayer,DataLayer-method}
\alias{combine,SingleCellAssay,SingleCellAssay-method}
\title{Combine two SingleCellAssay or derived classes}
\usage{
\S4method{combine}{DataLayer,DataLayer}(x, y, ...)

\S4method{combine}{SingleCellAssay,SingleCellAssay}(x, y, ...)
}
\arguments{
  \item{x}{One of the values.}

  \item{y}{A second value.}

  \item{...}{Any other objects of the same class as
  \code{x} and \code{y}.}
}
\description{
Combines two Single Cell-like objects provided they have
the same number of Features and Layers. The union of
columns from featureData will be taken The union (padded if
necessary with NA) will be taken from cellData.

combines two single cell assays provided they share a
common mapping TODO combining based on a common mapping may
be too restrictive. This may change depending on needs.
}

