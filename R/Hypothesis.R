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


## Courtesy of car, John Fox <jfox at mcmaster.ca>
## But not exported there
makeHypothesis <- function (cnames, hypothesis, rhs = NULL) 
{
    parseTerms <- function(terms) {
        component <- gsub("^[-\\ 0-9\\.]+", "", terms)
        component <- gsub(" ", "", component, fixed = TRUE)
        component
    }
    stripchars <- function(x) {
        x <- gsub("\\n", " ", x)
        x <- gsub("\\t", " ", x)
        x <- gsub(" ", "", x, fixed = TRUE)
        x <- gsub("*", "", x, fixed = TRUE)
        x <- gsub("-", "+-", x, fixed = TRUE)
        x <- strsplit(x, "+", fixed = TRUE)[[1]]
        x <- x[x != ""]
        x
    }
    char2num <- function(x) {
        x[x == ""] <- "1"
        x[x == "-"] <- "-1"
        as.numeric(x)
    }
    constants <- function(x, y) {
        with.coef <- unique(unlist(sapply(y, function(z) which(z == 
            parseTerms(x)))))
        if (length(with.coef) > 0) 
            x <- x[-with.coef]
        x <- if (is.null(x)) 
            0
        else sum(as.numeric(x))
        if (any(is.na(x))) 
            stop("The hypothesis \"", hypothesis, "\" is not well formed: contains bad coefficient/variable names.")
        x
    }
    coefvector <- function(x, y) {
        rv <- gsub(" ", "", x, fixed = TRUE) == parseTerms(y)
        if (!any(rv)) 
            return(0)
        if (sum(rv) > 1) 
            stop("The hypothesis \"", hypothesis, "\" is not well formed.")
        rv <- sum(char2num(unlist(strsplit(y[rv], x, fixed = TRUE))))
        if (is.na(rv)) 
            stop("The hypothesis \"", hypothesis, "\" is not well formed: contains non-numeric coefficients.")
        rv
    }
    if (!is.null(rhs)) 
        rhs <- rep(rhs, length.out = length(hypothesis))
    if (length(hypothesis) > 1) 
        return(rbind(Recall(cnames, hypothesis[1], rhs[1]), Recall(cnames, 
            hypothesis[-1], rhs[-1])))
    cnames_symb <- sapply(c("@", "#", "~"), function(x) length(grep(x, 
        cnames)) < 1)
    if (any(cnames_symb)) {
        cnames_symb <- head(c("@", "#", "~")[cnames_symb], 1)
        cnames_symb <- paste(cnames_symb, seq_along(cnames), 
            cnames_symb, sep = "")
        hypothesis_symb <- hypothesis
        for (i in order(nchar(cnames), decreasing = TRUE)) hypothesis_symb <- gsub(cnames[i], 
            cnames_symb[i], hypothesis_symb, fixed = TRUE)
    }
    else {
        stop("The hypothesis \"", hypothesis, "\" is not well formed: contains non-standard coefficient names.")
    }
    lhs <- strsplit(hypothesis_symb, "=", fixed = TRUE)[[1]]
    if (is.null(rhs)) {
        if (length(lhs) < 2) 
            rhs <- "0"
        else if (length(lhs) == 2) {
            rhs <- lhs[2]
            lhs <- lhs[1]
        }
        else stop("The hypothesis \"", hypothesis, "\" is not well formed: contains more than one = sign.")
    }
    else {
        if (length(lhs) < 2) 
            as.character(rhs)
        else stop("The hypothesis \"", hypothesis, "\" is not well formed: contains a = sign although rhs was specified.")
    }
    lhs <- stripchars(lhs)
    rhs <- stripchars(rhs)
    rval <- sapply(cnames_symb, coefvector, y = lhs) - sapply(cnames_symb, 
        coefvector, y = rhs)
    rval <- c(rval, constants(rhs, cnames_symb) - constants(lhs, 
        cnames_symb))
    names(rval) <- c(cnames, "*rhs*")
    rval
}
