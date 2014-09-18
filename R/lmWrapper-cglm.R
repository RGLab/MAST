## first four parameters give intercepts and scaling for components
bidx <- 5

## return log-likelihood given beta0, beta0prime, beta, s:
## Just sum of logistic and linear ll
hurdleConstrLL <- function(theta, X, y, pos){
    b0 <- theta['beta0']
    b0p <- theta['beta0prime']
    s <- theta['s']
    sigma <- theta['sigma']
    #sigma <- 1
    b <- theta[5:length(theta)]
    bGaus <- as.matrix(c(b0, b))
    bBinom <- as.matrix(c(b0p, b))*s

    gausLL <- sum(-(y[pos]- X[pos,,drop=FALSE] %*% bGaus)^2/(2*sigma^2)-log(sigma))
    binomEta <- X%*%bBinom
    binomLL <- sum(pos*(binomEta) - log(1+exp(binomEta)))
        cat(theta, gausLL, binomLL, '\n')
    -(gausLL+binomLL)
}

setMethod('initialize', 'ConstrainedGLMlike', function(.Object, ...){
    .Object <- callNextMethod()
    model.matrix(.Object) <- model.matrix(.Object@formula, .Object@design)
    .Object
})

makeGLMLike <- function(O, which){
    if(which=='C'){
        idx <- c(1, 5:length(O$par))
        s <- 1
        dispersion <- O$par['sigma']
    } else{
        idx <- c(2, 5:length(O$par))
        s <- O$par['s']
        dispersion <- 1
    }
    li <- list(coefficients=O$par[idx]*s, vcov=solve(O$hessian+diag(length(O$par))*1e-3)[idx,idx],dispersion=dispersion, H=O$hessian, s=s)
    
}

setMethod('fit', signature=c(object='ConstrainedGLMlike', response='missing'), function(object, response, silent=TRUE, ...){
    prefit <- .fit(object)
    if(!prefit){
        if(!silent) warning('No positive observations')
        return(object)
    }

    fitArgsC <- object@fitArgsC
    fitArgsD <- object@fitArgsD
    MM <- model.matrix(object)
    ## fixme: wont work with pivoting in LRT
    browser()
    theta <- setNames(rep(1, ncol(MM)+3), c('beta0', 'beta0prime', 's', 'sigma', colnames(MM)[-1]))
    #theta <- setNames(rep(1, ncol(MM)+2), c('beta0', 'beta0prime', 's', colnames(MM)[-1]))
    O <- optim(theta, hurdleConstrLL, X=MM, y=object@response, pos=pos, hessian=TRUE, method='BFGS')
    object@fitC <- makeGLMLike(O, 'C')
    object@fitD <- makeGLMLike(O, 'D')
    object@fitted <- c(C=TRUE, D=TRUE)
    object
})
