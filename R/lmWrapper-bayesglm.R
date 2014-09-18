##' @include AllClasses.R
##' @include AllGenerics.R

setMethod('fit', signature=c(object='BayesGLMlike', response='missing'), function(object, response, silent=TRUE, ...){
    prefit <- .fit(object)
    if(!prefit){
        if(!silent) warning('No positive observations')
        return(object)
    }

     fitArgsC <- object@fitArgsC
    fitArgsD <- object@fitArgsD
    object@fitC <- do.call(glm.fit, c(list(x=object@modelMatrix[pos,,drop=FALSE], y=object@response[pos]), fitArgsC))
    object@fitD <- do.call(bayesglm.fit, c(list(x=object@modelMatrix, y=pos*1, family=binomial()), fitArgsD))

    ## bayesglm doesn't correctly set the residual DOF
    object@fitC$df.residual <- sum(pos) - object@fitC$rank
    object@fitD$df.residual <- length(pos) - object@fitD$rank
    
    object@fitted <- c(C=object@fitC$converged & object@fitC$df.residual>0, D=object@fitD$converged & object@fitD$df.residual>0)
    object <- .dispersion(object)
    
    if(!silent & !all(object@fitted)) warning('At least one component failed to converge')
    object
})

