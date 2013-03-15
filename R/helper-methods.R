##' Compute the Et from the Ct
##'
##' Computes the Et value from the Ct value in an existing data frame and returns a new data frame with the Et column appended
##' @title computeEtFromCt
##' @param df a \code{data.frame}
##' @param column The name of the \code{Ct} column. A \code{character}. 'Ct' by default.
##' @param Cmax the maximum number of cycles performed. 40 by default.
##' @return A copy of \code{df} with the 'Et' column appended
##' @author Greg Finak
##' @export computeEtFromCt
##' @aliases computeEtFromCt
computeEtFromCt<-function(df,column='Ct',Cmax=40){
  within.data.frame(df, {Et <- Cmax-get(column); Et <- ifelse(is.na(Et), 0, Et)})
}

#'Convert a pair of character names to a list that can be placed into a data.table
#'
#'@param v a character vector
#'@return a call that can be evaluated
#'@export
todt<-function(v){
  as.call(c(list,sapply(v,as.name)))
}

#'Recursive combine
#'
#'Combines single cell assays recursively
#' @param dfs is a list of single-cell assays from an SCASet
#' @param ... ignored
#'@export
combine_recurse<-function (dfs=NULL, ...)
{
  if (length(dfs) == 2) {
    combine(dfs[[1]], y=dfs[[2]])
  }
  else {
    combine(dfs[[1]], Recall(dfs[-1]))
  }
}
