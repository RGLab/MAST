## Generics
setGeneric('fit', function(object, response, ...) standardGeneric('fit'))
setGeneric('coef', function(object, ...) standardGeneric('coef'))
#setGeneric('coefD', function(object) standardGeneric('coefD'))
setGeneric('lrTest', function(object, drop.terms) standardGeneric('lrTest'))
setGeneric('waldTest', function(object, hypothesis.matrix) standardGeneric('waldTest'))
## setGeneric('vcovC', function(object) standardGeneric('vcovC'))
## setGeneric('vcovD', function(object) standardGeneric('vcovD'))
setGeneric('vcov', function(object) standardGeneric('vcov'))
setGeneric('dof', function(object) standardGeneric('dof'))

## Classes
setClass('LMlike', slots=c(design='data.frame', fitC='ANY', fitD='ANY', response='ANY', fitted='logical', formula='formula'),     prototype=list(fitted =c(C=FALSE, D=FALSE)), validity=function(object){
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


setMethod('fit', signature=c(object='LMlike', response='vector'), function(object, response, ...){
    object@response <- response
    validObject(object)
    fit(object, ...)
})

setMethod('coef', signature=c(object='LMlike'), function(object, which, singular=TRUE, ...){
    stopifnot(which %in% c('C', 'D'))
    co <- if(which=='C') coef(object@fitC) else coef(object@fitD)
    if(!singular) co <- co[!is.na(co)]
    co
})



setMethod('summary', signature=c(object='LMlike'), function(object){
    print('Discrete\n============')
    print(coefD(object))
    print('Continuous\n============')
    print(coefC(object))
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
               df=dfC, 'P(Chisq>lambda)'=pchisq(lambdaC, df=dfC, lower.tail=FALSE))
    row.names(tab) <- c('cont', 'disc', 'hurdle')
    structure(tab, test=test)
}


##' Wald test for hurdle regression
##'
##' .. content for \details{} ..
##' @param object class \code{GLMlike}
##' @param hypothesis.matrix character vector or matrix following car::linearHypothesis
##' @seealso linearHypothesis
##' @return matrix with columns
##' @importFrom car linearHypothesis.default
##' @importFrom plyr laply
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
    dl <- -2*(l1-l0)
    df <- dof(object) - dof(fitnew)
    makeChiSqTable(dl, df, drop.terms)
})
