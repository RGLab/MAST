##' @include AllClasses.R
##' @include AllGenerics.R
setMethod('update', signature=c(object='GLMlike'), function(object, formula., ...){
    object <- callNextMethod(object, formula., ...)
    object@modelMatrix <- model.matrix(object@formula, object@design, ...)
    object
})

setMethod('vcov', signature=c(object='GLMlike'), function(object, which, ...){
    stopifnot(which %in% c('C', 'D'))
    if(which=='C') stats:::summary.glm(object@fitC, dispersion=object@fitC$dispersion)$cov.scaled else stats:::summary.glm(object@fitD)$cov.scaled
})

.dispersion <- function(object){
    if(object@fitted['C']){
        object@fitC$dispersionMLE <- object@fitC$deviance/(object@fitC$df.null+1)
        object@fitC$dispersion <- object@fitC$deviance/object@fitC$df.residual
    } else{
        object@fitC$dispersion <- NA
    }

    object
}


setMethod('fit', signature=c(object='GLMlike', response='missing'), function(object, response, silent=TRUE, ...){
    prefit <- .fit(object)
    if(!prefit){
        if(!silent) warning('No positive observations')
        return(object)
    }

    fitArgsC <- object@fitArgsC
    fitArgsD <- object@fitArgsD
    object@fitC <- do.call(glm.fit, c(list(x=object@modelMatrix[pos,], y=object@response[pos]), fitArgsC))
    object@fitD <- do.call(glm.fit, c(list(x=object@modelMatrix, y=pos*1, family=binomial()), fitArgsD))
    object@fitted <- c(C=object@fitC$converged & object@fitC$df.residual>0, D=object@fitD$converged & object@fitD$df.residual>0)

    object <- .dispersion(object)
    
    if(!silent & !all(object@fitted)) warning('At least one component failed to converge')    
    object
})

setMethod('initialize', 'GLMlike', function(.Object, ...){
    .Object <- callNextMethod()
    .Object@modelMatrix <- model.matrix(.Object@formula, .Object@design)
    .Object
})

setMethod('logLik', signature=c(object='GLMlike'), function(object){
    L <- c(C=0, D=0)
    if(object@fitted['C']){
        s2 <- object@fitC$dispersionMLE
        dev <- object@fitC$deviance
        N <- (object@fitC$df.null+1)
        L['C'] <- -.5*N*(log(s2*2*pi) +1)
    }

    if(object@fitted['D']){
         dev <- object@fitD$deviance
         L['D'] <- -dev/2
    }
    return(L)
})

setMethod('dof', signature=c(object='GLMlike'), function(object){
    c(C=length(coef(object, 'C', singular=FALSE)), D=length(coef(object, 'D', singular=FALSE)))
})
