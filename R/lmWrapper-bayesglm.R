##' @include AllClasses.R
##' @include AllGenerics.R
##' @import arm
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
    }
    
    
    object@fitC <- do.call(glm.fit, c(list(x=object@modelMatrix[pos,,drop=FALSE], y=object@response[pos]), fitArgsC))
    object@fitD <- do.call(bayesglm.fit, c(list(x=object@modelMatrix, y=pos*1, family=binomial()), fitArgsD))

    ## bayesglm doesn't correctly set the residual DOF
    object@fitC$df.residual <- sum(pos) - object@fitC$rank
    object@fitD$df.residual <- length(pos) - object@fitD$rank
    
    object@fitted <- c(C=object@fitC$converged &
                           object@fitC$df.residual>0, #kill unconverged or empty
                       D=object@fitD$converged &      #kill unconverged
                           (object@fitD$df.residual>0) & #note that we technically get a fit here, but it's probably not worth using
                               (min(sum(!pos), sum(pos))-object@fitD$rank)>0)
    object <- .dispersion(object)
    
    if(!silent & !all(object@fitted)) warning('At least one component failed to converge')
    object
})

setReplaceMethod('model.matrix', 'BayesGLMlike', function(object, value){
    oldcols <- dimnames(object@coefPrior)[[3]]
    newcols <- colnames(value)
    keepcols <- intersect(oldcols, newcols)
    if(length(object@coefPrior)>0){
        newprior <- .defaultPrior(newcols)
        newprior[,,keepcols] <- object@coefPrior[,,keepcols]
        object@coefPrior <- newprior
    }
    callNextMethod()
})
