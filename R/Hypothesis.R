##' @importFrom limma makeContrasts
##' @importFrom stringr str_replace
callName <- function(n=1){
    sc <- sys.calls()
    cl <- deparse(sc[[length(sc)-n]])
    str_replace(cl, '\\(.*\\)', '')
}

##' Describe a linear model hypothesis to be tested
##'
##' A \code{Hypothesis} can be any (set) linear combination of coefficients, (jointly) compared to zero.  Specify it as a character vector that can be parsed to yield the desired equalities ala \code{makeContrasts}.
##' Note that non-syntactic names are allowed, but you'll need to escape them with backticks "`".  In this case, your contrast needs to be parseable by R as an expression representing the comparison you desire.
##' @param hypothesis a character vector specifying a hypothesis, following makeContrasts, or a character vector naming coefficients to be dropped.
##' @return a Hypothesis
##' @export
##' @seealso zlm.SingleCellAssay waldTest lrTest makeContrasts
##  Eliminate boilerplate by dynamically inferring what our callname was 
Hypothesis <- function(hypothesis){
    whoami <- callName()
    new('Hypothesis', .Data=hypothesis)
}

generateHypothesis <- function(h, terms){
    stopifnot(inherits(h, 'Hypothesis') | inherits(h, 'CoefficientHypothesis'))

    if(inherits(h, 'Hypothesis')){
       trans <- makeContrasts2(contrasts=h@.Data, levels=terms)
       rownames(trans) <- terms
        
        ## makeContrasts can't handle non-syntactic names :-/
        ## So we'll use this instead
        ## trans <- makeHypothesis(cnames=terms,h@.Data)
        ## if(is.null(dim(trans))) trans <- as.matrix(trans)
        ## ## kill *rhs* column
        ## trans <- t(trans[-nrow(trans),,drop=FALSE])
        sd <- setdiff(rownames(trans), terms)
    } else{                             #CoefficientHypothesis
        trans <- h@.Data
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

##' @importFrom stringr str_replace_all str_c
## escape symbols with backticks
escapeSymbols <- function(text, warn=TRUE){
    hasBT <- str_detect(text, fixed('`'))
    if(any(hasBT)){
        if(warn) warning("Some symbols already contain backticks ('`').  Deleting backticks and hoping for the best.")
        text <- str_replace_all(text, fixed('`'), '')
    }
    hasSymbols <- str_detect(text, '[():+*/]|-')
    text[hasSymbols] <- str_c('`', text[hasSymbols], '`')
    text

}

## Adapted from limma, but allowing for non-syntactic names
##' @importFrom stringr str_detect
makeContrasts2 <- function (contrasts = NULL, levels, warn=TRUE) 
{
    if (is.factor(levels)) 
        levels <- levels(levels)
    if (!is.character(levels)) 
        levels <- colnames(levels)
    symbols <- str_detect(levels, '[():+*/]|-')
    if (any(symbols) && warn) 
        warning("Some levels contain symbols.  Be careful to escape these names with backticks ('`') when specifying contrasts.")
    n <- length(levels)
    if (n < 1) 
        stop("No levels to construct contrasts from")
    indicator <- function(i, n) {
        out <- rep(0, n)
        out[i] <- 1
        out
    }
    levelsenv <- new.env()
    for (i in 1:n) assign(levels[i], indicator(i, n), pos = levelsenv)
    if (!is.null(contrasts)) {
        e <- as.character(contrasts)
        ne <- length(e)
        cm <- matrix(0, n, ne, dimnames = list(Levels = levels, 
            Contrasts = e))
        if (ne == 0) 
            return(cm)
        for (j in 1:ne) {
            tryCatch( ej <- parse(text = e[j]), error=function(E) stop('Could not parse contrast ', e[j]))
            cm[, j] <- eval(ej, envir = levelsenv)
        }
        return(cm)
    }
}
