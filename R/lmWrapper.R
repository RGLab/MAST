setOldClass('glm')

setClass('GLMlike', representation=representation(modelMatrix='matrix', design='data.frame', formula='formula', fitC='ANY', fitD='ANY', response='ANY'), validity=function(object){
    if(length(object@response)>0){
        stopifnot(length(object@response)==nrow(object@design))
        stopifnot(length(object@response)==nrow(object@modelMatrix))
    }
})

setMethod('initialize', 'GLMlike', function(.Object, modelMatrix, design, formula, ...){
    .Object <- callNextMethod()
    if(missing(modelMatrix)){
        .Object@modelMatrix <- model.matrix(formula, design)
    }
    .Object@design <- design
    .Object@formula <- formula
    .Object
})

setGeneric('fit', function(object, response, ...) standardGeneric('fit'))
setGeneric('coefC', function(object) standardGeneric('coefC'))
setGeneric('coefD', function(object) standardGeneric('coefD'))
setGeneric('lrTest', function(object, drop.terms) standardGeneric('lrTest'))
setGeneric('waldTest', function(object, hypothesis.matrix) standardGeneric('waldTest'))


setMethod('show',  signature=c(object='GLMlike'), function(object){
    if(length(coefD(object)>0)) cat('Fitted' ) else cat('Unfitted')
    cat(class(object), ':', nrow(object@design), 'cases\n')
})

setMethod('fit', signature=c(object='GLMlike', response='missing'), function(object, response, ...){
    ## Assume design exists:
    ## Call lmFit with response and object
    pos <- object@response>0
    ## test if pos>0
    object@fitC <- glm.fit(object@modelMatrix[pos,], object@response[pos])
    object@fitD <- glm.fit(object@modelMatrix, pos*1, family=binomial())
    object@response
    object
    ## object@fitC <- glm(response[pos]>0~.,data=object@design[pos,])
    ## object@fitD <- glm(pos~.,data=object@design)
})

setMethod('fit', signature=c(object='GLMlike', response='vector'), function(object, response, ...){
    object@response <- response
    fit(object)
})

setMethod('coefC', signature=c(object='GLMlike'), function(object){
    coef(object@fitC)
})


setMethod('coefD', signature=c(object='GLMlike'), function(object){
    coef(object@fitD)
})

setMethod('summary', signature=c(object='GLMlike'), function(object){
    print(summary(object@fitD))
    print(summary(object@fitC))
})

setMethod('update', signature=c(object='GLMlike'), function(object, formula., ...){
object@modelMatrix <- model.matrix(formula., object@design)
object@formula <- formula.
object@fitC <- object@fitD <- numeric(0)
object
})

setMethod('lrTest', signature=c(object='GLMlike', drop.terms='character'), function(object, drop.terms){
    l0 <- logLik(object)
    F <- update.formula(object@formula, formula(sprintf(' ~. - %s', drop.terms)))
    U <- update(object, F)
    fitnew <- fit(U)
    l1 <- logLik(fitnew)
    dl <- -2*(l1-l0)
    df <- ncol(F@modelMatrix)-ncol(object@modelMatrix)
    ## Call table generating routine
})

setMethod('logLik', signature=c(object='GLMlike'), function(object){
    -.5*(object@fitC$deviance+object@fitD$deviance)
})

setMethod('waldTest', signature=c(object='GLMlike', hypothesis.matrix='character'), function(object, hypothesis.matrix){


})
