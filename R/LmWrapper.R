## Generics

##' fit a zero-inflated regression
##'
##' Given a design and formula, fit the zero inflated regression, storing the fits in slots
##' \code{fitC} and \code{fitD}
##' @param object inheriting from \code{LMlike}
##' @param response a vector, same length as the design, or if missing then use the current response
##' @param ... currently ignored
##' @return LMlike or subclass
##' @export
setGeneric('fit', function(object, response, ...) standardGeneric('fit'))

##' Coefficients of zero-infated
##'
##' Given a fitted LMlike, return the coefficients from discrete or continuous
##' @param object LMlike
##' @param which  character vector, one of "C" (continuous) or "D" (discrete) specifying which component should be returned
##' @param ... passed to methods
##' @return numeric vector
##' @export
setGeneric('coef', function(object, ...) standardGeneric('coef'))
#setGeneric('coefD', function(object) standardGeneric('coefD'))

##' Run a likelihood-ratio test
##'
##' Compares the change in likelihood between the current \code{formula} and one dropping terms in \code{drop.terms}.
##' Only complete terms can be tested at this time
##' @param object LMlike or subclass
##' @param drop.terms character vector of \code{formula} terms
##' @return array giving test statistics
##' @export
##' @seealso fit
##' @seealso waldTest
setGeneric('lrTest', function(object, drop.terms) standardGeneric('lrTest'))

##' Run a Wald test
##'
##' Run a Wald tests on discrete and continuous components
##' @param object LMlike or subclass
##' @param hypothesis.matrix argument suitable to be passed to car::lht
##' @return array giving test statistics
##' @export
##' @seealso fit
##' @seealso lrTest
##' @seealso lht
##' @importFrom car linearHypothesis.default
setGeneric('waldTest', function(object, hypothesis.matrix) standardGeneric('waldTest'))
## setGeneric('vcovC', function(object) standardGeneric('vcovC'))
## setGeneric('vcovD', function(object) standardGeneric('vcovD'))

##' Variance-covariance matrix for zero inflated
##'
##' Given a fitted LMlike, return the variance-covariance from discrete or continuous
##' @param object LMlike
##' @param which character vector, one of "C" (continuous) or "D" (discrete) specifying which component should be returned
##' @return matrix
##' @export
setGeneric('vcov', function(object) standardGeneric('vcov'))

##' Degrees of freedom of Zero inflated model
##'
##' @param object LMlike or subclass
##' @return vector giving the model degrees of freedom for continuous and discrete
##' @export
setGeneric('dof', function(object) standardGeneric('dof'))

## Classes
##' Linear Model-like Class
##'
##' Wrapper around modeling function to make them behave enough alike that Wald tests and Likelihood ratio are easy to do.
##' To implement a new type of zero-inflated model, extend this class.
##'
##' @section Slots:
##' \describe{
##' \item{design}{a data.frame from which variables are taken for the right hand side of the regression}
##' \item{fitC}{The continuous fit}
##' \item{fitD}{The discrete fit}
##' \item{response}{The left hand side of the regression}
##' \item{fitted}{A \code{logical} with components "C" and "D", TRUE if the respective component has converge}
##' \item{formula}{A \code{formula} for the regression}
##' \item{fitArgsC}{}
##' \item{fitArgsD}{Both \code{list}s giving arguments that will be passed to the fitter (such as convergence criteria or case weights)}
##' }
##' @seealso fit
##' @seealso coef
##' @seealso lrTest
##' @seealso waldTest
##' @seealso vcov
##' @seealso dof
##' @seealso logLik
##' @name LMlike-class
##' @docType class
setClass('LMlike', slots=c(design='data.frame', fitC='ANY', fitD='ANY', response='ANY', fitted='logical', formula='formula', fitArgsD='list', fitArgsC='list'),     prototype=list(fitted =c(C=FALSE, D=FALSE)), validity=function(object){
    stopifnot( all(c("C", "D") %in% names(object@fitted)))
    if(length(object@response)>0 && any(is.na(object@response))) stop('NAs not permitted in response')
})

setClass('GLMlike', contains='LMlike', slots=c(modelMatrix='matrix'), validity=function(object){
    if(length(object@response)>0){
        stopifnot(length(object@response)==nrow(object@design))
        #stopifnot(length(object@response)==nrow(object@modelMatrix))
    }},
    prototype=list(modelMatrix=matrix(nrow=0, ncol=0),  fitted =c(C=FALSE, D=FALSE))
)
setClass('LMERlike', contains='LMlike')

## Methods for LMlike
setMethod('show',  signature=c(object='LMlike'), function(object){
    if(all(object@fitted)){
        cat('Fitted Continuous and discrete')
    } else if(object@fitted['C']){
        cat('Fitted Continuous')
    }else if(object@fitted['D']){
        cat('Fitted Discrete')
    }else {
        cat('Unfitted')
    }
    cat(':', as.character(object@formula), '\n')
    cat(class(object), ':', nrow(object@design), ' cases\n', sep='')
})

.fit <- function(object){
    frame <- sys.frame(-1)
    positive <- object@response>0
    object@fitted <- c(C=FALSE, D=FALSE)
    assign('pos', positive, pos=frame)
    assign('object', object, pos=frame)
    
    return(any(positive))
}


setMethod('fit', signature=c(object='LMlike', response='vector'), function(object, response, silent=TRUE, fitArgsC=list(), fitArgsD=list(), ...){
    object@response <- response
    object@fitArgsC <- fitArgsC
    object@fitArgsD <- fitArgsD
    validObject(object)
    fit(object, silent=silent, ...)
})

setMethod('coef', signature=c(object='LMlike'), function(object, which, singular=TRUE, ...){
    stopifnot(which %in% c('C', 'D'))
    co <- if(which=='C') coef(object@fitC) else coef(object@fitD)
    if(!singular) co <- co[!is.na(co)]
    co
})



setMethod('summary', signature=c(object='LMlike'), function(object){
    print('===========Discrete============')
    print(coef(object, 'D'))
    print('==========Continuous===========')
    print(coef(object, 'C'))
})

setMethod('update', signature=c(object='LMlike'), function(object, formula., ...){
    object@formula <- update.formula(object@formula, formula.)
    object@fitC <- object@fitD <- numeric(0)
    object@fitted <- c(C=FALSE, D=FALSE)
    object
})


makeChiSqTable <- function(lambda, df, test){
    stopifnot(all(names(lambda) == c('C', 'D')))
    stopifnot(all(names(df) == c('C', 'D')))
    lambdaC <- c(lambda, sum(lambda))
    dfC <- c(df, sum(df))
    tab <- cbind(lambda=lambdaC,
               df=dfC, 'Pr(>Chisq)'=pchisq(lambdaC, df=dfC, lower.tail=FALSE))
    row.names(tab) <- c('cont', 'disc', 'hurdle')
    structure(tab, test=test)
}


setMethod('waldTest', signature=c(object='LMlike', hypothesis.matrix='ANY'), function(object, hypothesis.matrix){
    if(object@fitted['C']){
            C <- car::linearHypothesis.default(object@fitC, hypothesis.matrix=hypothesis.matrix, test='Chisq', vcov.=vcov(object, which='C'), coef.=coef(object, which='C', singular=FALSE), singular.ok=TRUE)[2,c('Df', 'Chisq'), drop=TRUE]
    }else{
        C <- list(Chisq=0, Df=0)
    }
    if(object@fitted['D']){
        D <- car::linearHypothesis.default(object@fitD, hypothesis.matrix=hypothesis.matrix, test='Chisq', vcov.=vcov(object, which='D'), coef.=coef(object, which='D', singular=FALSE), singular.ok=TRUE)[2,c('Df', 'Chisq'), drop=TRUE]
    }else{
        D <- list(Chisq=0, Df=0)
    }
    makeChiSqTable(c(C[['Chisq']], D[['Chisq']]), c(C[['Df']],D[['Df']]),hypothesis.matrix)
})

setMethod('lrTest', signature=c(object='LMlike', drop.terms='character'), function(object, drop.terms){
    l0 <- logLik(object)
    F <- update.formula(object@formula, formula(sprintf(' ~. - %s', drop.terms)))
    U <- update(object, F)
    fitnew <- fit(U)
    l1 <- logLik(fitnew)
    bothfitted <- object@fitted & fitnew@fitted
    dl <- ifelse(bothfitted, -2*(l1-l0), c(0, 0))
    df <- ifelse(bothfitted, dof(object) - dof(fitnew), c(0, 0))
    makeChiSqTable(dl, df, drop.terms)
})
