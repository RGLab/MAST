##' @include AllClasses.R
##' @include AllGenerics.R
setMethod('fit', signature=c(object='RidgeBGLMlike', response='missing'), function(object, response, silent=TRUE, ...){
    prefit <- .fit(object)
    if(!prefit){
        if(!silent) warning('No positive observations')
        return(object)
    }
    
    fitArgsC <- object@fitArgsC
    fitArgsD <- object@fitArgsD
    if(length(object@coefPrior)>0){
        fitArgsD$prior.mean <- object@coefPrior['loc', 'D',]
        fitArgsD$prior.scale <- object@coefPrior['scale', 'D',]
        fitArgsD$prior.df <- object@coefPrior['df', 'D', ]
        if(object@useContinuousBayes){
            fitArgsC$prior.mean <- object@coefPrior['loc', 'C',]
            fitArgsC$prior.scale <- object@coefPrior['scale', 'C',]
            fitArgsC$prior.df <- object@coefPrior['df', 'C', ]
        }
    }
    
    contFit <- if(object@useContinuousBayes) .bayesglm.fit else ridge.fit
    fitArgsC$lambda<-object@lambda
    object@fitC <- do.call(contFit, c(list(x=object@modelMatrix[pos,,drop=FALSE], y=object@response[pos],  weights=object@weightFun(object@response[pos])), fitArgsC))
    object@fitD <- do.call(.bayesglm.fit, c(list(x=object@modelMatrix, y=object@weightFun(object@response), family=binomial()), fitArgsD))
    
    object <- .glmDOF(object, pos)
    object <- .dispersion(object)
    
    if(!silent & !all(object@fitted)) warning('At least one component failed to converge')
    object
})

setReplaceMethod('model.matrix', 'RidgeBGLMlike', function(object, value){
    object <- callNextMethod()
    oldcols <- dimnames(object@coefPrior)[[3]]
    newcols <- colnames(model.matrix(object))
    keepcols <- intersect(oldcols, newcols)
    if(length(object@coefPrior)>0){
        newprior <- defaultPrior(newcols)
        newprior[,,keepcols] <- object@coefPrior[,,keepcols]
        object@coefPrior <- newprior
    }
    object
})


ridge.fit<-function (x, y, weights = rep(1, nobs), start = NULL, etastart = NULL, 
                     mustart = NULL, offset = rep(0, nobs), family = gaussian(), 
                     control = list(), intercept = TRUE,lambda=0.1) 
{
                                        #browser()
    control <- do.call("glm.control", control)
    G<-diag(lambda,NCOL(x))
    G[1,1]<-0
                                        #x[,-1]<-scale(x[,-1],center=FALSE,scale=TRUE)
    yorig<-y
    y<-scale(y,center=TRUE,scale=FALSE)
    x<-rbind(x,G)
    y<-c(y,rep(0,NCOL(x)))
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    ynames <- if (is.matrix(y)) 
                  rownames(y)
              else names(y)
    conv <- FALSE
    nobs <- NROW(y)
    nvars <- ncol(x)
    EMPTY <- nvars == 0
    if (is.null(weights)) 
        weights <- rep.int(1, nobs)
    else
        weights<-c(weights,rep(1,nobs-length(weights)))
    if (is.null(offset)) 
        offset <- rep.int(0, nobs)
    variance <- family$variance
    linkinv <- family$linkinv
    if (!is.function(variance) || !is.function(linkinv)) 
        stop("'family' argument seems not to be a valid family object", 
             call. = FALSE)
    dev.resids <- family$dev.resids
    aic <- family$aic
    mu.eta <- family$mu.eta
    unless.null <- function(x, if.null) if (is.null(x)) 
                                            if.null
                                        else x
    valideta <- unless.null(family$valideta, function(eta) TRUE)
    validmu <- unless.null(family$validmu, function(mu) TRUE)
    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }
    if (EMPTY) {
        eta <- rep.int(0, nobs) + offset
        if (!valideta(eta)) 
            stop("invalid linear predictor values in empty model", 
                 call. = FALSE)
        mu <- linkinv(eta)
        if (!validmu(mu)) 
            stop("invalid fitted means in empty model", call. = FALSE)
        dev <- sum(dev.resids(y, mu, weights))
        w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
        residuals <- (y - mu)/mu.eta(eta)
        good <- rep_len(TRUE, length(residuals))
        boundary <- conv <- TRUE
        coef <- numeric()
        iter <- 0L
    }
    else {
        coefold <- NULL
        eta <- if (!is.null(etastart)) 
                   etastart
               else if (!is.null(start)) 
                   if (length(start) != nvars) 
                       stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s", 
                                     nvars, paste(deparse(xnames), collapse = ", ")), 
                            domain = NA)
                   else {
                       coefold <- start
                       offset + as.vector(if (NCOL(x) == 1L) 
                                              x * start
                                          else x %*% start)
                   }
               else family$linkfun(mustart)
        mu <- linkinv(eta)
        if (!(validmu(mu) && valideta(eta))) 
            stop("cannot find valid starting values: please specify some", 
                 call. = FALSE)
        devold <- sum(dev.resids(y[1:(nobs-NCOL(x))], mu[1:(nobs-NCOL(x))], weights[1:(nobs-NCOL(x))]))
        boundary <- conv <- FALSE
        for (iter in 1L:control$maxit) {
            good <- weights > 0
            varmu <- variance(mu)[good]
            if (anyNA(varmu)) 
                stop("NAs in V(mu)")
            if (any(varmu == 0)) 
                stop("0s in V(mu)")
            mu.eta.val <- mu.eta(eta)
            if (any(is.na(mu.eta.val[good]))) 
                stop("NAs in d(mu)/d(eta)")
            good <- (weights > 0) & (mu.eta.val != 0)
            if (all(!good)) {
                conv <- FALSE
                warning(gettextf("no observations informative at iteration %d", 
                                 iter), domain = NA)
                break
            }
            z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
            w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
            ngoodobs <- as.integer(nobs - sum(!good))
            C_qr<-get("C_Cdqrls",getNamespace("stats"))
            fit <- .Call(C_qr, x[good, , drop = FALSE] * 
                               w, z * w, min(1e-07, control$epsilon/1000), check = FALSE)
            if (any(!is.finite(fit$coefficients))) {
                conv <- FALSE
                warning(gettextf("non-finite coefficients at iteration %d", 
                                 iter), domain = NA)
                break
            }
            if (nobs < fit$rank) 
                stop(sprintf(ngettext(nobs, "X matrix has rank %d, but only %d observation", 
                                      "X matrix has rank %d, but only %d observations"), 
                             fit$rank, nobs), domain = NA)
            start[fit$pivot] <- fit$coefficients
            eta <- drop(x %*% start)
            mu <- linkinv(eta <- eta + offset)
            dev <- sum(dev.resids(y[1:(nobs-NCOL(x))], mu[1:(nobs-NCOL(x))], weights[1:(nobs-NCOL(x))]))
            if (control$trace) 
                cat("Deviance = ", dev, " Iterations - ", iter, 
                    "\n", sep = "")
            boundary <- FALSE
            if (!is.finite(dev)) {
                if (is.null(coefold)) 
                    stop("no valid set of coefficients has been found: please supply starting values", 
                         call. = FALSE)
                warning("step size truncated due to divergence", 
                        call. = FALSE)
                ii <- 1
                while (!is.finite(dev)) {
                    if (ii > control$maxit) 
                        stop("inner loop 1; cannot correct step size", 
                             call. = FALSE)
                    ii <- ii + 1
                    start <- (start + coefold)/2
                    eta <- drop(x %*% start)
                    mu <- linkinv(eta <- eta + offset)
                    dev <- sum(dev.resids(y, mu, weights))
                }
                boundary <- TRUE
                if (control$trace) 
                    cat("Step halved: new deviance = ", dev, "\n", 
                        sep = "")
            }
            if (!(valideta(eta) && validmu(mu))) {
                if (is.null(coefold)) 
                    stop("no valid set of coefficients has been found: please supply starting values", 
                         call. = FALSE)
                warning("step size truncated: out of bounds", 
                        call. = FALSE)
                ii <- 1
                while (!(valideta(eta) && validmu(mu))) {
                    if (ii > control$maxit) 
                        stop("inner loop 2; cannot correct step size", 
                             call. = FALSE)
                    ii <- ii + 1
                    start <- (start + coefold)/2
                    eta <- drop(x %*% start)
                    mu <- linkinv(eta <- eta + offset)
                }
                boundary <- TRUE
                dev <- sum(dev.resids(y, mu, weights))
                if (control$trace) 
                    cat("Step halved: new deviance = ", dev, "\n", 
                        sep = "")
            }
            if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
                conv <- TRUE
                coef <- start
                break
            }
            else {
                devold <- dev
                coef <- coefold <- start
            }
        }
        if (!conv) 
            warning("glm.fit: algorithm did not converge", call. = FALSE)
        if (boundary) 
            warning("glm.fit: algorithm stopped at boundary value", 
                    call. = FALSE)
        eps <- 10 * .Machine$double.eps
        if (family$family == "binomial") {
            if (any(mu > 1 - eps) || any(mu < eps)) 
                warning("glm.fit: fitted probabilities numerically 0 or 1 occurred", 
                        call. = FALSE)
        }
        if (family$family == "poisson") {
            if (any(mu < eps)) 
                warning("glm.fit: fitted rates numerically 0 occurred", 
                        call. = FALSE)
        }
        nobs<-nobs-nvars
        weights<-weights[1:nobs]
        mu<-mu[1:nobs]
        n<-n[1:nobs]
        y<-y[1:(nobs)]
        eta<-eta[1:(nobs)]
        good<-good[1:(nobs)]
        ynames<-ynames[1:(nobs)]
        w<-w[1:nobs]
        fit$effects<-fit$effects[1:nobs]
        fit$qr<-fit$qr[1:(nobs),]
        dev <- sum(dev.resids(y, mu, weights))
        if (fit$rank < nvars) 
            coef[fit$pivot][seq.int(fit$rank + 1, nvars)] <- NA
        xxnames <- xnames[fit$pivot]
        residuals <- (y - mu)/mu.eta(eta)
        fit$qr <- as.matrix(fit$qr)
        nr <- min(sum(good), nvars)
        if (nr < nvars) {
            Rmat <- diag(nvars)
            Rmat[1L:nr, 1L:nvars] <- fit$qr[1L:nr, 1L:nvars]
        }
        else Rmat <- fit$qr[1L:nvars, 1L:nvars]
        Rmat <- as.matrix(Rmat)
        Rmat[row(Rmat) > col(Rmat)] <- 0
        names(coef) <- xnames
        colnames(fit$qr) <- xxnames
        dimnames(Rmat) <- list(xxnames, xxnames)
    }
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    wt <- rep.int(0, nobs)
    wt[good] <- (w^2)
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    if (!EMPTY) 
        names(fit$effects) <- c(xxnames[seq_len(fit$rank)], rep.int("", 
                                                                    sum(good) - fit$rank))
    wtdmu <- if (intercept) 
                 sum(weights * y)/sum(weights)
             else linkinv(offset)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
    rank <- if (EMPTY) 
                0
            else fit$rank
    resdf <- n.ok - rank
    aic.model <- aic(y, n, mu, weights, dev) + 2 * rank
                                        #coefficients.. add de-center the data
    coef[1]<-coef[1]+attr(scale(yorig),"scaled:center")
    mu<-mu+attr(scale(yorig),"scaled:center")
    list(coefficients = coef, residuals = residuals, fitted.values = mu, 
         effects = if (!EMPTY) fit$effects, R = if (!EMPTY) Rmat, 
         rank = rank, qr = if (!EMPTY) structure(fit[c("qr", "rank", 
                                                       "qraux", "pivot", "tol")], class = "qr"), family = family, 
         linear.predictors = eta, deviance = dev, aic = aic.model, 
         null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights, 
         df.residual = resdf, df.null = nulldf, y = yorig, converged = conv, 
         boundary = boundary)
}

