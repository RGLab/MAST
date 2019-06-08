##' MAST: Model-based Analysis of Single- cell Transcriptomics
##'
##' Methods for analysing single cell assay data using hurdle models.
##'
##' This packages provides data structures and functions for statistical analysis of single-cell assay data such as Fluidigm single cell gene expression assays.
##' @references Finak, et al.  MAST: a flexible statistical framework for assessing transcriptional changes and characterizing heterogeneity in single-cell RNA sequencing data.  Genome Biology (2015).
"_PACKAGE"

##' Vbeta Data Set
##' @docType data
##' @name vbeta
##' @rdname vbeta-dataset
##' @format a data frame with 11 columns.
##' Column \code{Ct} contains the cycle threshold, with NA denoting that the threshold was never crossed.  So it is inversely proportional to the log2 mRNA, and should be negated (and NAs set to zero) if it is used as a expression measurement for a \code{FluidigmAssay}.
NULL


##' Vbeta Data Set, FluidigmAssay
##' @docType data
##' @name vbetaFA
##' @rdname vbetaFA-dataset
##' @format a \code{FluidigmAssay} of the vbeta data set.
##' @seealso \code{\link{vbeta}}, \code{\link{FromFlatDF}}
NULL


##' MAITs data set, RNASeq
##' @docType data
##' @name maits
##' @rdname maits-dataset
##' @format a \code{list} containing an expression matrix (\code{expressionmat}), cell \code{cdat} and feature \code{fdat}.
##' @seealso \code{\link{FromMatrix}}
NULL


##' Predicted signatures
##' @docType data
##' @name predicted_sig
##' @rdname predicted_sig-dataset
##' @format A data frame of predicted gene expresion signatures for stimulated and unstimulated cells.
NULL

#' Defunct functions in package `MAST`
#' 
#' These functions are defunct or have been renamed.
#' @section Functions (and replacements, if available):
#' \describe{
#'   \item{filter}{mast_filter}
#' }
#' 
#' @aliases filter
#' @name MAST-defunct
NULL