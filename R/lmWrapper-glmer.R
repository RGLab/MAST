##' @include AllClasses.R
##' @include AllGenerics.R
setMethod('fit', signature=c(object='LMERlike', response='missing'), function(object, response, silent=TRUE, ...){
    ## Assume design exists:
    ## Call lmFit with response and object
    prefit <- .fit(object)
    if(!prefit){
        if(!silent) warning('No positive observations')
        return(object)
    }


    formC <- update.formula(object@formula, .response ~ .)
    formD <- update.formula(object@formula, .response>0 ~ .)
    dat <- cbind(.response=object@response, object@design)
    fitArgsC <- object@fitArgsC
    fitArgsD <- object@fitArgsD
    object@fitC <- do.call(lmer, c(list(formula=formC, data=dat[pos,], REML=FALSE), fitArgsC))
    if(!all(pos)){
        object@fitD <- do.call(glmer, c(list(formula=formD, data=dat, family=binomial()), fitArgsD))
        object@fitted['D'] <- length(object@fitD@optinfo$conv$lme)==0
    } 
    object@fitted['C'] <- length(object@fitC@optinfo$conv$lme)==0
    if(!silent & !all(object@fitted)) warning('At least one component failed to converge')
    object
})


setMethod('vcov', signature=c(object='LMERlike'), function(object, which, ...){
    stopifnot(which %in% c('C', 'D'))
    if(which=='C') vcov(object@fitC) else vcov(object@fitD)
})

## setMethod('vcovC', signature=c(object='LMERlike'), function(object){
##     vcov(object@fitC)
## })

## setMethod('coefC', signature=c(object='LMERlike'), function(object){
##     fixef(object@fitC)
## })


setMethod('coef', signature=c(object='LMERlike'), function(object, which, singular=TRUE, ...){
    stopifnot(which %in% c('C', 'D'))
    co <- if(which=='C') fixef(object@fitC) else fixef(object@fitD)
    if(!singular) co <- co[!is.na(co)]
    co
})

setMethod('logLik', signature=c(object='LMERlike'), function(object){
    setNames(ifelse(object@fitted, c(logLik(object@fitC), logLik(object@fitD)), c(0,0)), c('C', 'D'))
})

setMethod('dof', signature=c(object='LMERlike'), function(object){
    setNames(ifelse(object@fitted, c(attr(logLik(object@fitC), 'df'), attr(logLik(object@fitD), 'df')), c(0,0)), c('C', 'D'))

})
