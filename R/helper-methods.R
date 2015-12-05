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

##' Return a dependency-free version of a SingleCellAssay Object
##'
##' Breaks \code{sca} into a list of components `exprs` (expression matrix, cell X genes),
##' `cdata` (\code{data.frame} of cell annotations)
##' `fdata` (\code{data.frame} of feature annotations).
##' Optionally writes this list to an RDS file.
##' 
##' @param sca \code{SingleCellAssay} object
##' @param file optional \code{character} naming a file to which a .RDS will be written
##' @return a list (invisibly if \code{file} is set)
##' @export
jailBreak <- function(sca, file){
    components = list(exprs=exprs(sca), cdata=cData(sca), fdata=fData(sca))
    if(!missing(file)){
        saveRDS(components, file)
        invisible(components)
    } else{
        components
    }
}
