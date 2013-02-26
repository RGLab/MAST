\name{summary}
\alias{summary}
\title{Summarize expression parameters}
\arguments{
  \item{fd}{object inheriting from SingleCellAssay}

  \item{groups}{character vector of grouping variables}
}
\value{
  list of mu, pi and num.
}
\description{
  The mean of positive cells, mu, proportion of gene
  expression pi, and number of expressing cells per
  \code{groups} per gene is returned as a list
}

