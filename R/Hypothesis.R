##' @importFrom car makeHypothesis
setMethod('initialize', c(.Object='Hypothesis'), function(.Object, characterHypothesis, Terms, ...){
    .Object <- callNextMethod()
    #cat('Initialize', class(.Object), 'with .Data\n')
    #str(.Object@.Data)
    if(!missing(characterHypothesis)){
        mh <- makeHypothesis(Terms, characterHypothesis)
        if(is.null(dim(mh))) mh <- t(as.matrix(mh))
        .Object@.Data <- mh
    }
    validObject(.Object)
    return(.Object)
})

##' @importFrom stringr str_replace
callName <- function(n=1){
    sc <- sys.calls()
    cl <- deparse(sc[[length(sc)-n]])
    cl
    str_replace(cl, '\\(.*\\)', '')
}

##' Describe a linear model hypothesis to be tested
##'
##' A \code{Hypothesis} can be any linear combination of coefficients, compared to an arbitrary constant.
##' A \code{SimpleHypothesis} is a hypothesis for which terms are singly or jointly tested to be zero (generally the case in a t-test or F-test).
##' A \code{TermHypothesis} is a \code{SimpleHypothesis} for which entire factors are tested jointly.
##' @param hypothesis a character vector specifying a hypothesis, following the format of car::linearHypothesis
##' @return a function to be called with argument 'terms' giving the actual model terms
##' @export Hypothesis
##' @export SimpleHypothesis
##' @export TermHypothesis
##' @aliases SimpleHypothesis TermHypothesis
##' @seealso zlm.SingleCellAssay waldTest lrTest linearHypothesis
##' @importFrom plyr laply
Hypothesis <- function(hypothesis){
    whoami <- callName()
    if(!is(hypothesis, 'list')) hypothesis <- list(hypothesis)
    if(!all(laply(hypothesis, is, class2='character')  | laply(hypothesis, is,class2= 'matrix'))) stop("'hypothesis' must be list of symbolic character vectors or contrast matrices")
    f <- function(terms, Assign){
        lapply(hypothesis, function(x) new(whoami, characterHypothesis=x, Terms=terms, Assign=Assign))
    }
    class(f) <- c(class(f), 'hypothesisGenerator')
    f
}

TermHypothesis <- SimpleHypothesis <- Hypothesis

getTerms <- function(hypo){
    stopifnot(is(hypo, 'Hypothesis'))
    tested <- apply(hypo!=0, 2, sum)>0
    colnames(hypo)[tested]
}

contrastMatrix <- function(hypo, drop=FALSE){
    stopifnot(is(hypo, 'Hypothesis'))
    hypo[,seq_len(ncol(hypo)-1), drop=drop]
}

rhs <- function(hypo, drop=FALSE){
    stopifnot(is(hypo, 'Hypothesis'))
    hypo[,ncol(hypo), drop=drop]
}
