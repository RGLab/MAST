\name{SCASet}
\alias{SCASet}
\title{SCASet constructor}
\usage{
SCASet(dataframe, splitby, idvars = NULL, primerid = NULL,
  measurement = NULL, contentClass = "SingleCellAssay", ...)
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

  \item{contentClass}{a character, the name of the class
  being constructed within this SCASet. Defaults to
  SingleCellAssay. Other methods may pass in other classes,
  i.e. FluidigmAssay.}

  \item{...}{passed up to SingleCellAssay or other
  dynamically called constructor.}
}
\value{
SCASet
}
\description{
Constructs a SCASet
}
\details{
An SCASet is a list of SingleCellAssays or objects
inheriting from SingleCellAssay. The type of constructor
called is determined by the value of contentClass, which
should be the class of the SCA inheriting object contained
in this SCASet. Both the class and the constructor should
exist and have the same name. The code dynamically looks to
see if the a function with the same name exists, and
ASSUMES it is the constructor for the class. ##' ##' TODO
SCASet constructor should perhaps take a SingleCellAssay
class or FluidigmClass rather than a dataframe. Then we can
learn the class type for construction.
}
\note{
The dynamic lookup of the constructor could be made more
robust.
}

