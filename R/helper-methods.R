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

reraise <- function(err, convertToWarning=FALSE, silent=FALSE){
    if(silent) return(err)
    if(convertToWarning){
        warning(simpleWarning(message=err$message, call=err$call))
    } else{
        stop(simpleError(message=err$message, call=err$call))
    }
    return(err)
}

##' Selectively muffle warnings based on output
##'
##' @param expr an expression
##' @param regexp a regexp to be matched (with str_detect)
##' @return the result of expr
##' @export
hushWarning <- function(expr, regexp){
    withCallingHandlers(expr, warning=function(w){
        if(str_detect(conditionMessage(w), regexp)) invokeRestart("muffleWarning")
    })
}

##' Remove the left hand side (response) from a formula
##'
##' The order of terms will be rearrange to suit R's liking for hierarchy but otherwise the function should be idempotent for
##' @param Formula formula
##' @param warn Issue a warning if a response variable is found?
##' @return formula
##' @author Andrew
removeResponse <- function(Formula, warn=TRUE){
    charForm <- paste0(deparse(Formula, width.cutoff=500), collapse='')
    fsplit <- str_split_fixed(charForm, fixed('~'), 2)
    if(nchar(fsplit[1,1])>0 && warn) message("Ignoring LHS of formula (", fsplit[1,1], ') and using exprs(sca)')
    as.formula(paste0('~', fsplit[1,2])) 
}

## removeResponse <- function(Formula, warn=TRUE){
##     trm <- terms(Formula)
##     if(attr(trm, 'response')==1 && warn) warning("Ignoring LHS variable from formula", Formula)
##     reformulate(labels(trm), intercept=attr(trm, 'intercept'))
## }
