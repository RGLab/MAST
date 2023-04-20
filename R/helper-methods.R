##' Compute the Et from the Ct
##'
##' Computes the Et value from the Ct value in an existing data frame and returns a new data frame with the Et column appended
##' @param df a \code{data.frame}
##' @param column The name of the \code{Ct} column. A \code{character}. 'Ct' by default.
##' @param Cmax the maximum number of cycles performed. 40 by default.
##' @return A copy of \code{df} with the 'Et' column appended
##' @author Greg Finak
##' @export
##' @examples
##' data(vbeta)
##' vbeta <- computeEtFromCt(vbeta)
computeEtFromCt<-function(df,column='Ct',Cmax=40){
    within.data.frame(df, {Et <- Cmax-get(column); Et <- ifelse(is.na(Et), 0, Et)})
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
##' @examples
##' hushWarning(warning('Beware the rabbit'), 'rabbit')
##' hushWarning(warning('Beware the rabbit'), 'hedgehog')
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
    if(nchar(fsplit[1,1])>0 && warn) message("Ignoring LHS of formula (", fsplit[1,1], ') and using assay(sca)')
    as.formula(paste0('~', fsplit[1,2])) 
}

#' Combine lists, preferentially taking elements from x if there are duplicate names
#'
#' @param x list
#' @param y list
#' @examples
#' MAST:::meld_list_left(list(A=1, B=2), list(A = 0))
meld_list_left = function(x, y){
    unite = c(x, y)
    dups = nchar(names(unite)) & duplicated(names(unite))
    unite[!dups]
}

#' Instantiate a class, but warn rather than error for badly named slots
#'
#' @param classname `character` naming a class
#' @param ... slots in `classname`
#' @param extra named list giving other slots in `classname`
#'
#' @return `new(classname)`
#' @examples
#' MAST:::new_with_repaired_slots("SimpleList",  listData = list(x = LETTERS), 
#' extra = list(elementType = 'character', food = "tasty", beer = "cold"))
new_with_repaired_slots = function(classname, ..., extra){
  #... were actually named in caller, so assumed good.
  # extra are dot args, and might have been passed in erroneously by a caller.
  bad_slots = setdiff(names(extra), slotNames(classname))
  if(length(bad_slots) > 0){
    warning(sprintf("Dropping illegal slot(s) %s for class %s.  
                    This likely indicates a bug in an upstream package.", 
                    paste0(bad_slots, collapse = ', '), classname))
  }
  do.call(new, c(Class = classname, list(...), extra[setdiff(names(extra), bad_slots)]))
}