\name{threshold}
\alias{threshold}
\title{Estimate thresholds for positive expression}
\arguments{
  \item{nsa}{NanostringAssay object}

  \item{include.primers}{primeris to use for population
  estimates.  If missing, then all primers are used.}

  \item{exclude.primers}{primers to exclude from population
  estimates}

  \item{posteriorprob}{currently ignored}

  \item{clip}{currently ignored}

  \item{debug}{should extra fields helpful for examining
  the thresholding be returned as a list?}

  \item{location.strength}{scaling of prior on location}

  \item{pseudo.counts}{total strength of prior vs data}

  \item{hard.threshold}{location estimates less than this
  value will be treated as belonging to the "noise
  cluster"}
}
\value{
  modified nsa or list with elements nsa and debugging info
  regarding the clustering
}
\description{
  Estimates per-gene x unit thresholds for positive
  expression and truncates values below this threshold
  Uncertain values (in terms of posterior probability of
  membership) can be set to NA or rounded left or right
  Thresholds are estimated using a Gaussian mixture model
  with prior supplied by population estimates.
}

