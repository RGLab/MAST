setMethod('update', signature=c(object='GLMlike'), function(object, formula., ...){
    object <- callNextMethod(object, formula., ...)
    object@modelMatrix <- model.matrix(object@formula, object@design)
    object
})

setMethod('vcov', signature=c(object='GLMlike'), function(object, which, ...){
    stopifnot(which %in% c('C', 'D'))
    if(which=='C') stats:::summary.glm(object@fitC)$cov.scaled else stats:::summary.glm(object@fitD)$cov.scaled
})

## setMethod('vcovC', signature=c(object='GLMlike'), function(object){
##     stats:::summary.glm(object@fitC)$cov.scaled
## })

setMethod('fit', signature=c(object='GLMlike', response='missing'), function(object, response, ...){
    prefit <- .fit(object)
    if(!prefit) return(object)
    
    object@fitC <- glm.fit(object@modelMatrix[pos,], object@response[pos])
    object@fitD <- glm.fit(object@modelMatrix, pos*1, family=binomial())
    object@fitted <- c(C=object@fitC$converged, D=object@fitD$converged)
    object
})

setMethod('initialize', 'GLMlike', function(.Object, ...){
    .Object <- callNextMethod()
    .Object@modelMatrix <- model.matrix(.Object@formula, .Object@design)
    .Object
})

setMethod('logLik', signature=c(object='GLMlike'), function(object){
    setNames(-.5*(ifelse(object@fitted, c(object@fitD$deviance, object@fitC$deviance), c(0,0))), c('C', 'D'))
})

setMethod('dof', signature=c(object='GLMlike'), function(object){
    c(C=length(coef(object, 'C', singular=FALSE)), D=length(coef(object, 'D', singular=FALSE)))
})
