##' @include lmWrapper.R
##' @include AllClasses.R
setMethod('update', signature=c(object='GLMlike'), function(object, formula., ...){
    object <- callNextMethod(object, formula., ...)
    object@modelMatrix <- model.matrix(object@formula, object@design, ...)
    object
})

setMethod('vcov', signature=c(object='GLMlike'), function(object, which, ...){
    stopifnot(which %in% c('C', 'D'))
    if(which=='C') stats:::summary.glm(object@fitC)$cov.scaled else stats:::summary.glm(object@fitD)$cov.scaled
})

## setMethod('vcovC', signature=c(object='GLMlike'), function(object){
##     stats:::summary.glm(object@fitC)$cov.scaled
## })

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
    object@fitted <- c(C=object@fitC$converged, D=object@fitD$converged)
    if(!silent & !all(object@fitted)) warning('At least one component failed to converge')
    object
})

setMethod('initialize', 'GLMlike', function(.Object, ...){
    .Object <- callNextMethod()
    .Object@modelMatrix <- model.matrix(.Object@formula, .Object@design)
    .Object
})

setMethod('logLik', signature=c(object='GLMlike'), function(object){

    setNames(ifelse(object@fitted, dof(object) - c(object@fitC$aic/2-1,    #AIC is has extra DOF penalty for gaussian??  See getS3method('logLik', 'glm')
                                                object@fitD$aic/2), c(0,0)), c('C', 'D'))
})

setMethod('dof', signature=c(object='GLMlike'), function(object){
    c(C=length(coef(object, 'C', singular=FALSE)), D=length(coef(object, 'D', singular=FALSE)))
})
