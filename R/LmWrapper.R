## Methods for LMlike
setMethod('show',  signature=c(object='LMlike'), function(object){
    if(all(object@fitted)){
        cat('Fitted Continuous and discrete')
    } else if(object@fitted['C']){
        cat('Fitted Continuous')
    }else if(object@fitted['D']){
        cat('Fitted Discrete')
    }else {
        cat('Unfitted')
    }
    cat(':', as.character(object@formula), '\n')
    cat(class(object), ':', nrow(object@design), ' cases\n', sep='')
})

.fit <- function(object){
    frame <- sys.frame(-1)
    positive <- object@response>0
    object@fitted <- c(C=FALSE, D=FALSE)
    assign('pos', positive, pos=frame)
    assign('object', object, pos=frame)
    
    return(any(positive))
}


setMethod('fit', signature=c(object='LMlike', response='vector'), function(object, response, silent=TRUE, fitArgsC=list(), fitArgsD=list(), ...){
    object@response <- response
    object@fitArgsC <- fitArgsC
    object@fitArgsD <- fitArgsD
    validObject(object)
    fit(object, silent=silent, ...)
})

setMethod('coef', signature=c(object='LMlike'), function(object, which, singular=TRUE, ...){
    stopifnot(which %in% c('C', 'D'))
    co <- if(which=='C') coef(object@fitC) else coef(object@fitD)
    if(!singular) co <- co[!is.na(co)]
    co
})



setMethod('summary', signature=c(object='LMlike'), function(object){
    print('===========Discrete============')
    print(coef(object, 'D'))
    print('==========Continuous===========')
    print(coef(object, 'C'))
})

setMethod('update', signature=c(object='LMlike'), function(object, formula., ...){
    object@formula <- update.formula(object@formula, formula.)
    object@fitC <- object@fitD <- numeric(0)
    object@fitted <- c(C=FALSE, D=FALSE)
    object
})


makeChiSqTable <- function(lambda, df, test){
    stopifnot(all(names(lambda) == c('C', 'D')))
    stopifnot(all(names(df) == c('C', 'D')))
    lambdaC <- c(lambda, sum(lambda))
    dfC <- c(df, sum(df))
    tab <- cbind(lambda=lambdaC,
               df=dfC, 'Pr(>Chisq)'=pchisq(lambdaC, df=dfC, lower.tail=FALSE))
    row.names(tab) <- c('cont', 'disc', 'hurdle')
    structure(tab, test=test)
}


setMethod('waldTest', signature=c(object='LMlike', hypothesis.matrix='ANY'), function(object, hypothesis.matrix){
    if(object@fitted['C']){
            C <- car::linearHypothesis.default(object@fitC, hypothesis.matrix=hypothesis.matrix, test='Chisq', vcov.=vcov(object, which='C'), coef.=coef(object, which='C', singular=FALSE), singular.ok=TRUE)[2,c('Df', 'Chisq'), drop=TRUE]
    }else{
        C <- list(Chisq=0, Df=0)
    }
    if(object@fitted['D']){
        D <- car::linearHypothesis.default(object@fitD, hypothesis.matrix=hypothesis.matrix, test='Chisq', vcov.=vcov(object, which='D'), coef.=coef(object, which='D', singular=FALSE), singular.ok=TRUE)[2,c('Df', 'Chisq'), drop=TRUE]
    }else{
        D <- list(Chisq=0, Df=0)
    }
    makeChiSqTable(c(C[['Chisq']], D[['Chisq']]), c(C[['Df']],D[['Df']]),hypothesis.matrix)
})

setMethod('lrTest', signature=c(object='LMlike', drop.terms='character'), function(object, drop.terms){
    l0 <- logLik(object)
    F <- update.formula(object@formula, formula(sprintf(' ~. - %s', drop.terms)))
    U <- update(object, F)
    fitnew <- fit(U)
    l1 <- logLik(fitnew)
    bothfitted <- object@fitted & fitnew@fitted
    dl <- ifelse(bothfitted, -2*(l1-l0), c(0, 0))
    df <- ifelse(bothfitted, dof(object) - dof(fitnew), c(0, 0))
    makeChiSqTable(dl, df, drop.terms)
})

setMethod('residuals', signature=c(object='LMlike'), function(object, type='response', which, ...){
    which <- match.arg(which, c('Discrete', 'Continuous', 'Marginal'))
    RD <- residuals(object@fitD, type=type)
    RC <- residuals(object@fitC, type=type)
    if(which=='Discrete') return(RD)
    if(which=='Continuous') return(RC)
    if(which=='Marginal'){
        if(type != 'response') warning("Marginal residuals probably don't make sense unless predicting on the response scale")
        ## Zero inflated residuals
        RC <- object@response - predict(object@fitC, newx=object@design)
        return(RC*RD)
    }
})
