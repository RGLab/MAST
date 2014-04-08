setOldClass('glm')

setClass('GLMlike', representation=representation(modelMatrix='data.frame', design='matrix', formula='formula', fitC='glm', fitD='glm')) ##, validity=function(object){
##     nrow(object@response)==nrow(object@design)
## })


setMethod('fit', signature=c(object='GLMlike', response='vector'), function(object, ...){
    ## Assume design exists:
    ## Call lmFit with response and object
    ## set fitC and fitD


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

setMethod('lrTest', signature=c(object='GLMlike', drop.terms='character'), function(object, drop.terms){
    #call fit or ANOVA updated formula
})

setMethod('waldTest', signature=c(object='GLMlike', hypothesis.matrix='character'), function(object, hypothesis.matrix){


})
