## Generics
setGeneric('fit', function(object, response, ...) standardGeneric('fit'))
setGeneric('coefC', function(object) standardGeneric('coefC'))
setGeneric('coefD', function(object) standardGeneric('coefD'))
setGeneric('lrTest', function(object, drop.terms) standardGeneric('lrTest'))
setGeneric('waldTest', function(object, hypothesis.matrix) standardGeneric('waldTest'))
setGeneric('vcovC', function(object) standardGeneric('vcovC'))
setGeneric('vcovD', function(object) standardGeneric('vcovD'))
setGeneric('dof', function(object) standardGeneric('dof'))

## Classes
setClass('LMlike', slots=c(design='data.frame', fitC='ANY', fitD='ANY', response='ANY', fitted='logical', formula='formula'),     prototype=list(fitted =c(C=FALSE, D=FALSE)), validity=function(object) stopifnot( all(c("C", "D") %in% names(object@fitted))))

setClass('GLMlike', contains='LMlike', slots=c(modelMatrix='matrix'), validity=function(object){
    if(length(object@response)>0){
        stopifnot(length(object@response)==nrow(object@design))
        stopifnot(length(object@response)==nrow(object@modelMatrix))
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
    cat(' ', as.character(obj@formula), '\n')
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
    fit(object, ...)
})

setMethod('coefC', signature=c(object='LMlike'), function(object){
    coef(object@fitC)
})


setMethod('coefD', signature=c(object='LMlike'), function(object){
    coef(object@fitD)
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
    data.frame(test=test, lambda=lambda, df=df, 'P(Chisq>lambda)'=pchisq(lambda, df=df, lower.tail=FALSE), check.names=FALSE)

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
            C <- car::linearHypothesis.default(object@fitC, hypothesis.matrix=hypothesis.matrix, test='Chisq', vcov.=vcovC(object), coef.=coefC(object))[2,c('Df', 'Chisq'), drop=TRUE]
    }else{
        C <- list(Chisq=0, Df=0)
    }
    if(object@fitted['D']){
        D <- car::linearHypothesis.default(object@fitD, hypothesis.matrix=hypothesis.matrix, test='Chisq', vcov.=vcovD(object), coef.=coefD(object))[2,c('Df', 'Chisq'), drop=TRUE]
    }else{
        D <- list(Chisq=0, Df=0)
    }
    makeChiSqTable(C[['Chisq']] + D[['Chisq']], C[['Df']]+D[['Df']],hypothesis.matrix)
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
