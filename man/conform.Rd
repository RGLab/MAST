\name{conform}
\alias{conform}
\title{conform}
\arguments{
  \item{dl}{DataLayer}

  \item{other}{Another object}

  \item{compareLayers}{Should the number of layers be
  tested?}
}
\value{
  bitmask containing number of dimensions that agree
}
\description{
  Do two objects conform in dimension and type?
}
\details{
  Returns false if: other does not inherit from Matrix
  Otherwise returns a bitmask, which is the sum of: 1:
  objects have same number of rows 2: objects have same
  number of columns 4: objects have same number of layers
  (when compareLayers is TRUE)
}
\author{
  andrew
}

