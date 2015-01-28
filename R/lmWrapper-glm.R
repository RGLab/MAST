##' @include AllClasses.R
##' @include AllGenerics.R

setMethod('initialize', 'GLMlike', function(.Object, ...){
    .Object <- callNextMethod()
    model.matrix(.Object) <- model.matrix(.Object@formula, .Object@design)
    .Object
})


## This is pinch point (up to 10% of computation time can be spent here)
setMethod('vcov', signature=c(object='GLMlike'), function(object, which, ...){
    stopifnot(which %in% c('C', 'D'))
    vc <- object@defaultVcov
    if(which=='C' & object@fitted['C']){
        vc2 <- stats:::summary.glm(object@fitC, dispersion=object@fitC$dispersion)$cov.scaled
    } else if(which=='D' & object@fitted['D']){
        vc2 <- stats:::summary.glm(object@fitD)$cov.scaled
    } else{
        vc2 <- numeric()
    }
    ok <- colnames(vc2)
    vc[ok,ok] <- vc2
    vc
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

setMethod('fit', signature=c(object='GLMlike', response='missing'), function(object, response, silent=TRUE, ...){
    prefit <- .fit(object)
    if(!prefit){
        if(!silent) warning('No positive observations')
        return(object)
    }

    fitArgsC <- object@fitArgsC
    fitArgsD <- object@fitArgsD
    object@fitC <- do.call(glm.fit, c(list(x=object@modelMatrix[pos,,drop=FALSE], y=object@response[pos], weights=object@weights[pos]), fitArgsC))
    object@fitD <- hushWarning(do.call(glm.fit, c(list(x=object@modelMatrix, y=object@weights, family=binomial()), fitArgsD)), fixed("non-integer #successes in a binomial glm"))
    ## needed so that residuals dispatches more correctly
    class(object@fitD) <- c('glm', class(object@fitD))
    object@fitted <- c(C=object@fitC$converged & object@fitC$df.residual>0, D=object@fitD$converged & object@fitD$df.residual>0)
    ## cheap additional test for convergence
    ## object@fitted['D'] <- object@fitted['D'] & (object@fitD$null.deviance >= object@fitD$deviance)
    ## update dispersion, possibly shrinking by prior
    object <- .dispersion(object)
    
    if(!silent & !all(object@fitted)) warning('At least one component failed to converge')    
    object
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


## make a row matrix
rowm <- function(C, D){
    x <- c(C=NA, D=NA)
    try({if(is.null(C) | missing(C))
        C <- NA
    if(is.null(D) | missing(D))
        D <- NA
    x <- c(C=C, D=D)
     }, silent=TRUE)
    ## dim(x) <- c(1, length(x))
    ## colnames(x) <- c('C', 'D')
    x
}

torowm <- function(x){
     ## dim(x) <- c(1, length(x))
     ## colnames(x) <- c('C', 'D')
     x
}

setMethod('summarize', signature=c(object='GLMlike'), function(object, ...){
    coefC <- coef(object, which='C')
    coefD <- coef(object, which='D')
    ## make sure covariance matrices are constant size
    ## if it's not fitted, then we'll throw an error here
    vcC <- vcov(object, 'C')
    vcD <- vcov(object, 'D')
    
    list(coefC=coefC, vcovC=vcC,
          deviance=rowm(C=object@fitC$deviance, D=object@fitD$deviance),
          df.null=rowm(C=object@fitC$df.null, D=object@fitD$df.null),
          df.resid=rowm(C=object@fitC$df.residual, D=object@fitD$df.residual),
          dispersion=rowm(C=object@fitC$dispersionMLE, D=object@fitD$dispersion),
          dispersionNoshrink=rowm(C=object@fitC$dispersionMLENoShrink, D=object@fitD$dispersion),
          loglik=torowm(logLik(object)),
          coefD=coefD, vcovD=vcD, converged=torowm(object@fitted))
  })
