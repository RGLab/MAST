\name{zlm.SingleCellAssay}
\alias{zlm.SingleCellAssay}
\title{zlm.SingleCellAssay}
\usage{
  zlm.SingleCellAssay(formula, sca, lm.fun = glm,
    hypothesis.matrix, type = "Wald", hypo.fun = NULL,
    keep.zlm = FALSE, .parallel = FALSE, .drop = TRUE,
    .inform = FALSE, ...)
}
\arguments{
  \item{formula}{a formula with the measurement variable on
  the LHS and predictors present in cData on the RHS}

  \item{sca}{SingleCellAssay object}

  \item{lm.fun}{a function accepting lm-style arguments and
  a family argument}

  \item{hypothesis.matrix}{names of coefficients to test in
  lht form}

  \item{hypo.fun}{a function taking a model as input and
  returning output suitable for hypothesis.matrix}

  \item{keep.zlm}{should the model objects be kept}

  \item{.parallel}{run fits using parallel processing.
  must have doParallel}

  \item{...}{passed to lm.fun}
}
\value{
  either an array of tests (one per primer) or a list
}
\description{
  zero-inflated regression for SingleCellAssay
}
\details{
  Fits a hurdle model in \code{formula} (linear for et>0),
  logistic for et==0 vs et>0. Conducts likelihood ratio
  tests using hypothesis.matrix.

  A \code{list} of \code{data.frame}s, is returned, with
  one \code{data.frame} per tested predictor. Rows of each
  \code{data.frame} are genes, the columns give the value
  of the LR test and its P-value, and the sum of the
  T-statistics for each level of the factor (when the
  predictor is categorical).
}

