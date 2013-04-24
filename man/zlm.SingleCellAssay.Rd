\name{zlm.SingleCellAssay}
\alias{zlm.SingleCellAssay}
\title{zlm.SingleCellAssay}
\usage{
  zlm.SingleCellAssay(formula, sca, scope,
    keep.zlm = FALSE, ...)
}
\arguments{
  \item{formula}{a formula with the measurement variable on
  the LHS and predictors present in cData on the RHS}

  \item{sca}{SingleCellAssay object}

  \item{scope}{(optional) a formula giving the size of the
  smaller model to be fit.  If omitted, each predictor will
  be dropped in turn.}

  \item{...}{passed to lm and glm.}
}
\value{
  a \code{list} of \code{data.frame}, one per tested
  predictor.  See details.
}
\description{
  zero-inflated regression for SingleCellAssay
}
\details{
  Fits a hurdle model in \code{formula} (linear for et>0),
  logistic for et==0 vs et>0. Conducts likelihood ratio
  tests for each predictor in \code{formula} that does not
  appear in \code{scope}.

  A \code{list} of \code{data.frame}s, is returned, with
  one \code{data.frame} per tested predictor. Rows of each
  \code{data.frame} are genes, the columns give the value
  of the LR test and its P-value, and the sum of the
  T-statistics for each level of the factor (when the
  predictor is categorical).
}

