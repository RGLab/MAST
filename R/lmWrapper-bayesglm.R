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
    if(length(object@coefPrior)>0){
        fitArgsD$prior.mean <- object@coefPrior['loc', 'D',]
        fitArgsD$prior.scale <- object@coefPrior['scale', 'D',]
        fitArgsD$prior.df <- object@coefPrior['df', 'D', ]
        if(object@useContinuousBayes){
                 fitArgsC$prior.mean <- object@coefPrior['loc', 'C',]
                 fitArgsC$prior.scale <- object@coefPrior['scale', 'C',]
                 fitArgsC$prior.df <- object@coefPrior['df', 'C', ]
             }
    }
    
    contFit <- if(object@useContinuousBayes) .bayesglm.fit else glm.fit
    
    object@fitC <- do.call(contFit, c(list(x=object@modelMatrix[pos,,drop=FALSE], y=object@response[pos],  weights=object@weights[pos]), fitArgsC))
    object@fitD <- hushWarning(
        do.call(.bayesglm.fit, c(list(x=object@modelMatrix, y=object@weights, family=binomial()), fitArgsD)),
        fixed("non-integer #successes in a binomial glm"))
 

    object <- .glmDOF(object, pos)
    object <- .dispersion(object)
    
    if(!silent & !all(object@fitted)) warning('At least one component failed to converge')
    object
})

setReplaceMethod('model.matrix', 'BayesGLMlike', function(object, value){
    object <- callNextMethod()
    oldcols <- dimnames(object@coefPrior)[[3]]
    newcols <- colnames(model.matrix(object))
    keepcols <- intersect(oldcols, newcols)
    if(length(object@coefPrior)>0){
        newprior <- defaultPrior(newcols)
        newprior[,,keepcols] <- object@coefPrior[,,keepcols]
        object@coefPrior <- newprior
    }
    object
})
