##' Compute the Et from the Ct
##'
##' Computes the Et value from the Ct value in an existing data frame and returns a new data frame with the Et column appended
##' @title computeEtFromCt
##' @param df a \code{data.frame}
##' @param column The name of the \code{Ct} column. A \code{character}. 'Ct' by default.
##' @return A copy of \code{df} with the 'Et' column appended
##' @author Greg Finak
##' @export computeEtFromCt
##' @aliases computeEtFromCt
computeEtFromCt<-function(df,column='Ct'){
  within.data.frame(df, {Et <- 40-get(column); Et <- ifelse(is.na(Et), 0, Et)})
}
