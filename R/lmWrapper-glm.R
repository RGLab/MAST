##' @include AllClasses.R
##' @include AllGenerics.R

setMethod('vcov', signature=c(object='GLMlike'), function(object, which, ...){
    stopifnot(which %in% c('C', 'D'))
    if(which=='C') stats:::summary.glm(object@fitC, dispersion=object@fitC$dispersion)$cov.scaled else stats:::summary.glm(object@fitD)$cov.scaled
})

## dispersion calculations for glm-like fitters
.dispersion <- function(object){
    object@fitC$dispersionMLE <- object@fitC$dispersion <- NA
    if(object@fitted['C']){
        df.total <- object@fitC$df.null+1
        df.residual <- object@fitC$df.residual
        ## Save unshrunken
        dMLEns <- object@fitC$deviance/df.total
        dns <- object@fitC$deviance/df.residual
        object@fitC$dispersionMLENoShrink <- dMLEns
        object@fitC$dispersionNoShrink <- dns

        ## Now shrink default component
        object@fitC$dispersionMLE <- (dMLEns*df.total + object@priorVar*object@priorDOF)/(df.total+object@priorDOF)
        object@fitC$dispersion <- (dns*df.residual+object@priorVar*object@priorDOF)/(df.residual+object@priorDOF)   
    }

    object
}

.residualsD <- function(object){
    if(object@fitted['D']){
        object@fitD$residuals <- (object@response>0)*1 - object@fitD$fitted
    } else{
        object@fitD$residuals <- NA
    }
    object
}

## Performance enhancement: consider adding a 'which' argument, because for LRT with contrasts, we need to refit only the continuous in principal
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
    ## needed so that residuals dispatches more correctly
    class(object@fitD) <- c('glm', class(object@fitD))
    object@fitted <- c(C=object@fitC$converged & object@fitC$df.residual>0, D=object@fitD$converged & object@fitD$df.residual>0)
    ## update dispersion, possibly shrinking by prior
    object <- .dispersion(object)
    
    if(!silent & !all(object@fitted)) warning('At least one component failed to converge')    
    object
})

setMethod('initialize', 'GLMlike', function(.Object, ...){
    .Object <- callNextMethod()
    .Object@modelMatrix <- model.matrix(.Object@formula, .Object@design)
    .Object
})


setMethod('logLik', signature=c(object='GLMlike'), function(object){
    L <- c(C=0, D=0)
    if(object@fitted['C']){
        s2 <- object@fitC$dispersionMLE
        dev <- object@fitC$deviance
        N <- (object@fitC$df.null+1)
        L['C'] <- -.5*N*(log(s2*2*pi) +1)
    }

    if(object@fitted['D']){
         dev <- object@fitD$deviance
         L['D'] <- -dev/2
    }
    return(L)
})

setMethod('dof', signature=c(object='GLMlike'), function(object){
    c(C=length(coef(object, 'C', singular=FALSE)), D=length(coef(object, 'D', singular=FALSE)))
})

setMethod('residuals', signature=c(object='GLMlike'), function(object, type='response', which, ...){
    which <- match.arg(which, c('Discrete', 'Continuous', 'Marginal'))
    if(type != 'response') stop("Only type='response' residuals implemented for GLMlike")
    PD <- object@fitD$fitted
    RD <- object@fitD$residuals <- (object@response>0)*1 -PD
    ## May contain NAs for non-estimible coefficients, set to zero
    coefC <- coef(object, which='C', singular=TRUE)
    coefC[is.na(coefC)] <- 0
    PC <- object@modelMatrix %*% coefC
    RC <- (object@response - PC)[object@response>0]
    if(which=='Discrete') return(RD)
    if(which=='Continuous') return(RC)

    if(which=='Marginal'){
        if(type != 'response') warning("Marginal residuals probably don't make sense unless predicting on the response scale")               
        return(object@response-PC*PD)
    }
})
