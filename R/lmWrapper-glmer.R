setMethod('fit', signature=c(object='LMERlike', response='missing'), function(object, response, ...){
    ## Assume design exists:
    ## Call lmFit with response and object
    prefit <- .fit(object)
    if(!prefit) return(object)

    formC <- update.formula(object@formula, .response ~ .)
    formD <- update.formula(object@formula, .response>0 ~ .)
    dat <- cbind(.response=object@response, object@design)
    object@fitC <- lmer(formC, data=dat, ...)
    if(!all(pos)){
        object@fitD <- glmer(formD, data=dat, family=binomial(), ...)
        object@fitted['D'] <- length(object@fitD@optinfo$conv$lme)==0
    } 
    object@fitted['C'] <- length(object@fitC@optinfo$conv$lme)==0
    object
})


setMethod('vcovD', signature=c(object='LMERlike'), function(object){
    vcov(object@fitD)
})

setMethod('vcovC', signature=c(object='LMERlike'), function(object){
    vcov(object@fitC)
})

setMethod('coefC', signature=c(object='LMERlike'), function(object){
    fixef(object@fitC)
})


setMethod('coefD', signature=c(object='LMERlike'), function(object){
    fixef(object@fitD)
})

setMethod('logLik', signature=c(object='LMERlike'), function(object){
    sum(ifelse(object@fitted, c(logLik(object@fitC), logLik(object@fitD)), c(0,0)))
})

setMethod('dof', signature=c(object='LMERlike'), function(object){
    sum(ifelse(object@fitted, c(attr(logLik(object@fitC), 'df'), attr(logLik(object@fitD), 'df')), c(0,0)))

})
