setMethod('update', signature=c(object='GLMlike'), function(object, formula., ...){
    object <- callNextMethod(object, formula., ...)
    object@modelMatrix <- model.matrix(object@formula, object@design)
    object
})

setMethod('vcovD', signature=c(object='GLMlike'), function(object){
    CP <- crossprod(object@modelMatrix, diag(object@fitD$weights) %*% object@modelMatrix)
    solve(CP)
})

setMethod('vcovC', signature=c(object='GLMlike'), function(object){
    s2 <- sum(object@fitC$residuals^2)/object@fitC$df.residual
    mm <- object@modelMatrix[object@response>0,]
    solve(crossprod(mm)) * s2
})

setMethod('fit', signature=c(object='GLMlike', response='missing'), function(object, response, ...){
    prefit <- .fit(object)
    if(!prefit) return(object)
    
    object@fitC <- glm.fit(object@modelMatrix[pos,], object@response[pos])
    object@fitD <- glm.fit(object@modelMatrix, pos*1, family=binomial())
    object@fitted <- c(C=object@fitC$converged, D=object@fitD$converged)
    object
})

setMethod('initialize', 'GLMlike', function(.Object, design, formula, ...){
    .Object <- callNextMethod()
    if(!missing(formula)){
        .Object@modelMatrix <- model.matrix(formula, design)
        .Object@formula <- formula
    }
    .Object@design <- design
    .Object
})

setMethod('logLik', signature=c(object='GLMlike'), function(object){
    -.5* sum(ifelse(object@fitted, c(object@fitC$deviance, object@fitD$deviance), c(0,0)))
})

setMethod('dof', signature=c(object='GLMlike'), function(object){
    sum(object@fitted)*ncol(object@design)
})
