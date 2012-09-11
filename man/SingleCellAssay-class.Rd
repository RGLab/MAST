\docType{class}
\name{SingleCellAssay-class}
\alias{FluidigmAssay-class}
\alias{SingleCellAssay-class}
\title{SingleCellAssay represents an arbitrary single cell assay
It is meant to be flexible and can (and will) be subclassed to represent specific assay
types like Fluidigm and others. It should be constructed using the \code{SingleCellAssay} constructor, or ideally the \code{SCASet} constructor.
mapNames for the SingleCellAssay class are in the object \code{SingleCellAssay:::SingleCellAssayMapNames}
mapNames for the FluidigmAssay class are in the object \code{SingleCellAssay:::FluidigmMapNames}
}
\section{Slots}{
\describe{
  \item{featureData}{an \code{AnnotatedDataFrame} that describes feature-level metadata (i.e. genes)}
  \item{phenoData}{an \code{AnnotatedDataFrame} that describes the phenotype-level metadata (i.e. subject or experimental unit)}
  \item{cellData}{an \code{AnnotatedDataFrame} that describes the cell-level metadata (i.e. per individual cell)}
  \item{mapNames}{a \code{character} vector that describes some mandatory fields used by the class to map data from the raw file to the object. These are defined in the package, class definintion and subclasses.}
  \item{mapping}{a named \code{character} vector that maps \code{mapNames} to column names in the raw data file or data frame. This provides some flexibility for changing file formats and future assay types.}
  \item{description}{a \code{data.frame}}
  \item{wellKey}{A unique key that identifies the well, which USUALLY contains a single cell being measured.}
  \item{id}{a vector of type \code{character} that identifies the set of columns acting as a primary key to uniquely identify a single-cell or single-well across all wells / cells / assays / subjects / conditions in the data set.}
  \item{env}{an environment that will hold the data}
}}
\description{
  SingleCellAssay represents an arbitrary single cell assay
  It is meant to be flexible and can (and will) be
  subclassed to represent specific assay types like
  Fluidigm and others. It should be constructed using the
  \code{SingleCellAssay} constructor, or ideally the
  \code{SCASet} constructor. mapNames for the
  SingleCellAssay class are in the object
  \code{SingleCellAssay:::SingleCellAssayMapNames} mapNames
  for the FluidigmAssay class are in the object
  \code{SingleCellAssay:::FluidigmMapNames} }
  \section{Slots}{ \describe{ \item{featureData}{an
  \code{AnnotatedDataFrame} that describes feature-level
  metadata (i.e. genes)} \item{phenoData}{an
  \code{AnnotatedDataFrame} that describes the
  phenotype-level metadata (i.e. subject or experimental
  unit)} \item{cellData}{an \code{AnnotatedDataFrame} that
  describes the cell-level metadata (i.e. per individual
  cell)} \item{mapNames}{a \code{character} vector that
  describes some mandatory fields used by the class to map
  data from the raw file to the object. These are defined
  in the package, class definintion and subclasses.}
  \item{mapping}{a named \code{character} vector that maps
  \code{mapNames} to column names in the raw data file or
  data frame. This provides some flexibility for changing
  file formats and future assay types.}
  \item{description}{a \code{data.frame}} \item{wellKey}{A
  unique key that identifies the well, which USUALLY
  contains a single cell being measured.} \item{id}{a
  vector of type \code{character} that identifies the set
  of columns acting as a primary key to uniquely identify a
  single-cell or single-well across all wells / cells /
  assays / subjects / conditions in the data set.}
  \item{env}{an environment that will hold the data} }
}

