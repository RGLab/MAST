\name{plot.threshold}
\alias{plot.threshold}
\title{Plot NanoString log Counts vs cluster ID}
\usage{
\method{plot}{threshold}(thresholdedNanoString, primerids)
}
\arguments{
  \item{thresholdedNanoString}{output from
  thresholdNanoString, debug=TRUE}

  \item{primerids}{character vector of primerids to plot}
}
\value{
ggplot object
}
\description{
Histogram of log Counts by primerid, colored by cluster
}

