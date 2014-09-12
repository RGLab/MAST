##' @importFrom limma makeContrasts
##' @importFrom stringr str_replace
callName <- function(n=1){
    sc <- sys.calls()
    cl <- deparse(sc[[length(sc)-n]])
    str_replace(cl, '\\(.*\\)', '')
}

##' Describe a linear model hypothesis to be tested
##'
##' A \code{Hypothesis} can be any linear combination of coefficients, compared to zero.  Specify it as a character vector that can be parsed to yield the desired equalities ala \code{makeContrasts}.
##' A \code{CoefficientHypothesis} is a hypothesis for which terms are singly or jointly tested to be zero (generally the case in a t-test or F-test), by dropping coefficients from the model.
##' @param hypothesis a character vector specifying a hypothesis, following makeContrasts, or a character vector naming coefficients to be dropped.
##' @return a Hypothesis with a "transformed" component
##' @export Hypothesis
##' @export CoefficientHypothesis
##' @aliases Hypothesis CoefficientHypothesis
##' @seealso zlm.SingleCellAssay waldTest lrTest linearHypothesis
##  Eliminate boilerplate by dynamically inferring what our callname was 
Hypothesis <- CoefficientHypothesis <- function(hypothesis){
    whoami <- callName()
    new(whoami, .Data=hypothesis)
}

generateHypothesis <- function(h, terms){
    stopifnot(inherits(h, 'Hypothesis') | inherits(h, 'CoefficientHypothesis'))

    if(inherits(h, 'Hypothesis')){
        trans <- makeContrasts(contrasts=h, levels=terms)
        sd <- setdiff(rownames(h), terms)
    } else{                             #CoefficientHypothesis
        trans <- h
        sd <- setdiff(h, terms)
    }
    if(length(sd)>0) stop("Term(s) '", paste(sd, ','), "' not found.\nTerms available: ", paste(terms, ", "))

    h@transformed <- trans
    h
}

listType <- function(alist){
    types <- lapply(alist, function(x) class(x)[1])
    if(length(unique(types))>1) stop("Not 'atomic' list")
    types[1]
}


## want a contrast matrix in order to find reduced design matrix
## but design matrix size changes per gene
## so represent contrasts symbolically until we evaluate the gene

## Or..get maximum contrast matrix, then downscale as needed.


