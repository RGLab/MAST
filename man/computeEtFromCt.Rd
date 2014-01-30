\name{computeEtFromCt}
\alias{computeEtFromCt}
\title{computeEtFromCt}
\usage{
computeEtFromCt(df, column = "Ct", Cmax = 40)
}
\arguments{
  \item{df}{a \code{data.frame}}

  \item{column}{The name of the \code{Ct} column. A
  \code{character}. 'Ct' by default.}

  \item{Cmax}{the maximum number of cycles performed. 40 by
  default.}
}
\value{
A copy of \code{df} with the 'Et' column appended
}
\description{
Compute the Et from the Ct
}
\details{
Computes the Et value from the Ct value in an existing data
frame and returns a new data frame with the Et column
appended
}
\author{
Greg Finak
}

