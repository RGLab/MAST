\name{threshold}
\alias{threshold}
\title{Estimate thresholds for positive expression}
\arguments{
  \item{nsa}{NanostringAssay object}

  \item{groups}{groups to apply thresholding}

  \item{thresholds}{data.frame of thresholds of groups x
  genes user may specify}

  \item{posteriorprob}{Min posterior probability of cluster
  membership for an observation to be truncated}

  \item{clip}{Should values with uncertain posterior probs
  be clipped}
}
\value{
  modifies nsa in place
}
\description{
  Estimates per-gene x unit thresholds for positive
  expression and truncates values below this threshold
  Uncertain values (in terms of posterior probability of
  membership) can be set to NA or rounded left or right
  Thresholds are estimated using a Gaussian mixture model
  with prior supplied by population estimates.
}

