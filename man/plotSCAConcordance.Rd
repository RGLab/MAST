\name{plotSCAConcordance}
\alias{plotSCAConcordance}
\title{Concordance plots of filtered single vs n-cell assays}
\usage{
  plotSCAConcordance(SCellAssay, NCellAssay,
    filterCriteria = list(nOutlier = 2, sigmaContinuous = 9, sigmaProportion = 9),
    groups = NULL)
}
\arguments{
  \item{SCellAssay}{is a SingleCellAssay for the 1-cell per
  well assay}

  \item{NCellAssay}{is a SingleCellAssay for the n-cell per
  well assay}

  \item{filterCriteria}{is a list of filtering criteria to
  apply to the SCellAssay and NCellAssay}

  \item{groups}{is a character vector naming the group
  within which to perform filtering. NULL by default.}
}
\description{
  This will generate the types of plots shown in the
  bioinformatics manuscript
}

