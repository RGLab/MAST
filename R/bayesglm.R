## from ARM version 1.7-7 (2015-04-07)
## code originally by Gelman and Su
## http://cran.r-project.org/web/packages/arm/index.html
##' @importFrom stats .checkMFClasses as.formula binomial contrasts<- cov2cor delete.response deviance
##' @importFrom stats df.residual dt ecdf
##' @importFrom stats median  na.pass napredict naresid nobs offset optim
##' @importFrom stats optimize p.adjust pchisq predict pt qnorm qt quantile relevel setNames sigma t.test
##' @importFrom stats terms terms.formula  weighted.residuals
##' @importFrom stats family fitted formula gaussian glm.control glm.fit lm.fit update.formula model.frame model.matrix.default


.bayesglm.fit <- function (x, y, weights = rep(1, nobs), start = NULL,
                           etastart = NULL, mustart = NULL, offset = rep(0, nobs),
                           family = gaussian(),
                           control = glm.control(), intercept = TRUE,
                           prior.mean = 0,
                           prior.scale = NULL,
                           prior.df = 1,
                           prior.mean.for.intercept = 0,
                           prior.scale.for.intercept = NULL,
                           prior.df.for.intercept = 1,
                           min.prior.scale=1e-12,
                           scaled = TRUE, print.unnormalized.log.posterior=FALSE, Warning=TRUE)
{
######### INITIALIZE #######
    nobs <- NROW(y)
    nvars <- NCOL(x)
    conv <- FALSE
    EMPTY <- nvars == 0

    output <- .bayesglm.fit.initialize.priors (family, prior.scale, prior.scale.for.intercept, nvars,
                                               prior.mean, prior.mean.for.intercept, intercept, prior.df, prior.df.for.intercept)
    prior.scale <- output$prior.scale
    prior.mean <- output$prior.mean
    prior.df <- output$prior.df

    prior.scale <- .bayesglm.fit.initialize.priorScale (scaled, family, prior.scale, prior.scale.for.intercept, y, nvars, x, min.prior.scale)
    output <- .bayesglm.fit.initialize.x (x, nvars, nobs, intercept, scaled)
    x <- output$x
    xnames <- output$xnames
    x.nobs <- output$x.nobs

    output <- .bayesglm.fit.initialize.other (nobs, y, weights, offset)
    ynames <- output$ynames
    weights <- output$weights
    offset <- output$offset

    family <- .bayesglm.fit.initialize.family (family)
    ## TODO: DL -- I would put this inside initialize.family, but this does some magic setting of variables.
    ##             What I can do is push an environment into the function and wrap it in.

    if (is.null(mustart)){
        eval(family$initialize)
    }
    else {
        mustart.keep <- mustart
        eval(family$initialize)
        mustart <- mustart.keep
    }
########################
    if (EMPTY) {
        eta <- rep.int(0, nobs) + offset
        if (!family$valideta(eta)){
            stop("invalid linear predictor values in empty model")
        }
        mu <- family$linkinv(eta)
        if (!family$validmu(mu)){
            stop("invalid fitted means in empty model")
        }
        devold <- sum(family$dev.resids(y, mu, weights))
        w <- ((weights * family$mu.eta(eta)^2)/family$variance(mu))^0.5
        residuals <- (y - mu)/family$mu.eta(eta)
        good <- rep(TRUE, length(residuals))
        boundary <- conv <- TRUE
        coef <- numeric(0)
        iter <- 0
    }
    else {
        ## 3.1 ##
        output <- .bayesglm.fit.loop.setup (etastart, start = start, nvars, xnames, offset, x, nobs, family,
                                            mustart, eta, weights, conv, prior.scale, y, x.nobs)
        coefold <- output$Coefold
        start <- output$Start
        eta <- output$eta
        mu <- output$mu
        devold <- output$devold
        boundary <- output$boundary
        prior.sd <- output$prior.sd
        dispersion <- output$dispersion
        dispersionold <- output$dispersionold
##########

        ## 3.2 ## check the link between 3.1 and 3.2
        output.new <- .bayesglm.fit.loop.main.ideal (control, x, y, nvars, nobs, weights, offset,
                                                     intercept, scaled,
                                                     start = start, etastart = eta, mustart = mu, dispersion = dispersion,
                                                     family,
                                                     prior.mean, prior.scale=prior.sd, prior.df, 
                                                     print.unnormalized.log.posterior,
                                                     Warning)
        fit <- output.new$fit
        good <- output.new$good
        z <- output.new$z
        w <- output.new$w
        ngoodobs <- output.new$ngoodobs
        prior.scale <- output.new$prior.scale
        prior.sd <- output.new$prior.sd
        eta <- output.new$eta
        mu <- output.new$mu
        dev <- output.new$dev
        dispersion <- output.new$dispersion
        start <- output.new$start
        coef <- output.new$coef
        conv <- output.new$conv
        iter <- output.new$iter
        boundary <- output.new$boundary
        Rmat <- output.new$Rmat
        residuals <- output.new$residuals
    }

    ## 4 ##
    output <- .bayesglm.fit.cleanup (
        ynames = ynames,
        residuals = residuals,
        mu = mu, eta = eta, nobs = nobs,
        weights= weights, w = w, good = good,
        linkinv = family$linkinv,
        dev.resids = family$dev.resids,
        y = y,
        intercept = intercept,
        fit = fit, offset = offset, EMPTY = EMPTY,
        dev =dev, aic = family$aic
    )



    residuals <- output$residuals
    mu <- output$mu
    eta <- output$eta
    wt <- output$wt
    weights <- output$weights
    y <- output$y
    wtdmu <- output$wtdmu
    nulldev <- output$nulldev
    n.ok <- output$n.ok
    nulldf <- output$nulldf
    rank <- output$rank
    resdf <- output$resdf
    aic.model <- output$aic.model
#######
    list(coefficients = coef,
         residuals = residuals,
         fitted.values = mu,
         effects = if (!EMPTY) fit$effects,
         R = if (!EMPTY) Rmat,
         rank = rank,
         qr = if (!EMPTY) structure(getQr(fit)[c("qr", "rank", "qraux", "pivot", "tol")], class = "qr"),
         family = family,
         linear.predictors = eta,
         deviance = dev,
         aic = aic.model,
         null.deviance = nulldev,
         iter = iter,
         weights = wt,
         prior.weights = weights,
         df.residual = resdf,
         df.null = nulldf,
         y = y,
         converged = conv,
         boundary = boundary,
         prior.mean = prior.mean,
         prior.scale = prior.scale,
         prior.df = prior.df,
         prior.sd = prior.sd,
         dispersion = dispersion)
}

.bayesglm.fit.initialize.priors <- function (family, prior.scale, prior.scale.for.intercept, nvars,
                                             prior.mean, prior.mean.for.intercept, intercept, prior.df, prior.df.for.intercept) {
##### 12.13 ####
    if (is.null(prior.scale)){
        prior.scale <- 2.5
        if (family$link == "probit"){
            prior.scale <- prior.scale*1.6
        }
    }

    if (is.null(prior.scale.for.intercept)){
        prior.scale.for.intercept <- 10
        if (family$link == "probit"){
            prior.scale.for.intercept <- prior.scale.for.intercept*1.6
        }
    }
################

    if (intercept) {
        nvars <- nvars - 1
    }
    
    if (length(prior.mean) == 1) {
        prior.mean <- rep(prior.mean, nvars)
    }
    else if (length(prior.mean) != nvars) {
        stop("Invalid length for prior.mean")
    }

    if (length(prior.scale)==1){
        prior.scale <- rep(prior.scale, nvars)
    }
    else if (length(prior.scale) != nvars) {
        stop("Invalid length for prior.scale")
    }

    if (length(prior.df) == 1) {
        prior.df <- rep(prior.df, nvars)
    }
    else if (length(prior.df) != nvars) {
        stop("Invalid length for prior.df")
    }

    if (intercept) {
        prior.mean <- c(prior.mean.for.intercept, prior.mean)
        prior.scale <- c(prior.scale.for.intercept, prior.scale)
        prior.df <- c(prior.df.for.intercept, prior.df)
    }

    list (prior.mean=prior.mean,
          prior.scale=prior.scale,
          prior.df=prior.df)
}

.bayesglm.fit.initialize.priorScale <- function (scaled, family, prior.scale, prior.scale.for.intercept, y, nvars, x, min.prior.scale) {
    if (scaled) {
        if (family$family == "gaussian"){
            prior.scale <- prior.scale * 2 * sd(y)
            prior.scale.for.intercept <- prior.scale.for.intercept * 2 * sd(y)
        }

        for (j in 1:nvars) {
            x.obs <- x[, j]
            x.obs <- x.obs[!is.na(x.obs)]
            num.categories <- length(unique(x.obs))
            x.scale <- 1
            if (num.categories == 2) {
                x.scale <- max(x.obs) - min(x.obs)
            }
            else if (num.categories > 2) {
                x.scale <- 2 * sd(x.obs)
            }
            prior.scale[j] <- prior.scale[j]/x.scale
            if (prior.scale[j] < min.prior.scale){
                prior.scale[j] <- min.prior.scale
                warning ("prior scale for variable ", j,
                         " set to min.prior.scale = ", min.prior.scale,"\n")
            }
        }
    }
    return (prior.scale)
}

.bayesglm.fit.initialize.x <- function (x, nvars, nobs, intercept, scaled) {
    x <- as.matrix (rbind(x, diag(nvars)))
    xnames <- dimnames(x)[[2]]
    x.nobs <- x[1:nobs, ,drop=FALSE]
                                        # TODO: **** I moved it here because I think it's right, and changed it to x.nobs
    if (intercept & scaled) {
        x[nobs+1,] <- colMeans(x.nobs)
    }
    list (x=x, xnames=xnames, x.nobs=x.nobs)
}

.bayesglm.fit.initialize.family <- function (family, mustart, enviroment) {
    if (!is.function(family$variance) || !is.function(family$linkinv)){
        stop("'family' argument seems not to be a valid family object")
    }
    if (is.null(family$valideta)){
        family$valideta <- function(eta) TRUE
    }

    if (is.null(family$validmu)){
        family$validmu <- function(mu) TRUE
    }

    return (family)
}

.bayesglm.fit.initialize.other <- function (nobs, y, weights, offset) {
    if (is.matrix(y)){
        ynames <- rownames(y)
    }
    else{
        ynames <- names(y)
    }

    if (is.null(weights)){
        weights <- rep.int(1, nobs)
    }

    if (is.null(offset)){
        offset <- rep.int(0, nobs)
    }

    list (ynames=ynames,
          weights=weights,
          offset=offset)
}

.bayesglm.fit.loop.setup <- function (etastart, start, nvars, xnames, offset, x, nobs, family,
                                      mustart, eta, weights, conv, prior.scale, y, x.nobs) {
    coefold <- NULL
    if (!is.null(etastart)) {
        eta <- etastart
    } else if (!is.null(start)) {
        if (length(start) != nvars){
            if(start==0&length(start)==1){
                start <- rep(0, nvars)
                eta <- offset + as.vector(ifelse((NCOL(x) == 1), x.nobs[,1]*start, x.nobs %*% start))
            }else{
                stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
                              nvars, paste(deparse(xnames), collapse = ", ")), domain = NA)
            }
        }else {
            coefold <- start
            eta <- offset + as.vector(ifelse((NCOL(x) == 1), x.nobs[,1]*start, x.nobs %*% start))
        }
    }
    else {
        eta <- family$linkfun(mustart)
    }

    mu <- family$linkinv(eta)
    if (!isTRUE(family$validmu(mu) && family$valideta(eta)))
        stop("cannot find valid starting values: please specify some")
    devold <- sum(family$dev.resids(y, mu, weights))
    boundary <- conv <- FALSE
    prior.sd <- prior.scale
                                        #========Andy 2008.7.8=============
    dispersion <- ifelse((family$family %in% c("poisson", "binomial")),  1, var(y)/10000)
                                        #==================================
    dispersionold <- dispersion
    list (  Start = start,
          Coefold = coefold,
          eta=eta,
          mu=mu,
          devold=devold,
          boundary=boundary,
          prior.sd=prior.sd,
          dispersion=dispersion,
          dispersionold=dispersionold)
}

.bayesglm.fit.loop.initializeState <- function (start, etastart, mustart, dispersion, offset, x.nobs, var.y, nvars, family, weights, prior.sd, y) {
    coefold <- NULL
    if (!is.null(etastart)) {
        eta <- etastart
    }
    else if (!is.null(start)) {
        coefold <- as.matrix (start)
        eta <- drop (x.nobs %*% start) + offset
        ## TODO: should be able to drop this. coefold is the original start.
                                        #coefold <- start
    }
    else {
        eta <- family$linkfun(mustart)
    }
                                        #    if (family$family %in% c("poisson", "binomial")) {
                                        #        dispersion <- 1
                                        #    } else{
                                        #        dispersion <- var.y / 10000
                                        #    }
    mu <- family$linkinv(eta)
    mu.eta.val <- family$mu.eta(eta)
    dev <- sum (family$dev.resids (y, mu, weights))
    list (  Start = start, Coefold = coefold, eta=eta,
          mu=mu,
          mu.eta.val=mu.eta.val,
          varmu=family$variance(mu),
          good=(weights > 0) & (mu.eta.val != 0),
          dispersion=dispersion,
          dev=dev,
          conv=FALSE,
          boundary=FALSE,
          prior.sd=prior.sd)
}



.bayesglm.fit.loop.validateInputs <- function(etastart, start, nvars, xnames) {
    if (is.null(etastart) & !is.null(start) & length(start) != nvars) {
        stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
                      nvars, paste(xnames, collapse = ", ")), domain = NA)
    }
}


.bayesglm.fit.loop.validateState <- function (state, family, control, iter, dispersionold, devold) {
    if (all(!state$good)) {
        warning("no observations informative at iteration ", state$iter)
        return (FALSE)
    }
    if (is.null(state$fit) == FALSE && any (!is.finite(state$fit$coefficients))) {
        warning("non-finite coefficients at iteration ", iter)
        return (FALSE)
    }
    if (any(is.na(state$varmu[state$good]))){
        stop("NAs in V(mu)")
    }
    if (any(state$varmu[state$good] == 0)){
        stop("0s in V(mu)")
    }

    if (any(is.na(state$mu.eta.val[state$good]))){
        stop("NAs in d(mu)/d(eta)")
    }

    if (is.infinite(state$dispersion)){
        stop("dispersion blows up due to divergence!")
    }

    if (iter > 1 & abs(state$dev - devold)/(0.1 + abs(state$dev)) <  control$epsilon &
        abs(state$dispersion - dispersionold)/(0.1 + abs(state$dispersion)) < control$epsilon) {
        return (FALSE)
    }
    return (TRUE)
}

.bayesglm.fit.loop.updateState <- function (state, priors, family, #fortran.call.parameters,
                                            offset, weights,
                                            y, x, x.nobs, nvars, nobs,
                                            intercept, scaled, control) {
    z <- (state$eta[state$good] - offset[state$good]) + (y[state$good] - state$mu[state$good]) / state$mu.eta.val[state$good]
    z.star <- c(z, priors$mean)
    w <- sqrt((weights[state$good] * state$mu.eta.val[state$good]^2)/state$varmu[state$good])
    w.star <- c(w, sqrt(state$dispersion)/priors$scale)
    good.star <- c(state$good, rep(TRUE, nvars))
    fit <- lm.fit(x = as.matrix(x[good.star, ])*w.star, y = z.star*w.star)
    start <- state$Start
    coefold <- state$Start

                                        #fit <- .Fortran("dqrls",
                                        #            qr = x[good.star, ] * w.star,
                                        #            n = sum (good.star),
                                        #            p = nvars,
                                        #            y = w.star * z.star,
                                        #            ny = fortran.call.parameters$ny,
                                        #            tol = fortran.call.parameters$tol,
                                        #            coefficients = fortran.call.parameters$coefficients,
                                        #            residuals = fortran.call.parameters$residuals,
                                        #            effects = fortran.call.parameters$effects,
                                        #            rank = fortran.call.parameters$rank,
                                        #            pivot = fortran.call.parameters$pivot,
                                        #            qraux = fortran.call.parameters$qraux,
                                        #            work = fortran.call.parameters$work,
                                        #            PACKAGE = "base")


                                        #state$prior.sd <- priors$scale
    if (all (priors$df==Inf) == FALSE) {
        colMeans.x <- colMeans (x.nobs)
        centered.coefs <- fit$coefficients
        if(NCOL(x.nobs)==1){
            V.coefs <- chol2inv(fit$qr$qr[1:nvars])
        }else{
            V.coefs <- chol2inv(fit$qr$qr[1:nvars, 1:nvars, drop = FALSE])
        }
        diag.V.coefs <- diag(V.coefs)
        sampling.var <- diag.V.coefs
                                        # DL: sampling.var[1] <- crossprod(colMeans.x, V.coefs) %*% colMeans.x
        if(intercept & scaled){
            centered.coefs[1] <- sum(fit$coefficients*colMeans.x)
            sampling.var[1] <- crossprod (crossprod(V.coefs, colMeans.x), colMeans.x)
        }
        sd.tmp <- ((centered.coefs - priors$mean)^2 + sampling.var * state$dispersion + priors$df * state$prior.sd^2)/(1 + priors$df)
        sd.coef <- sqrt(sd.tmp)
        state$prior.sd[priors$df != Inf] <- sd.coef[priors$df != Inf]
    }

    if(NCOL(x.nobs)==1){
        predictions <- x.nobs * fit$coefficients
    }else{
        predictions <- x.nobs %*% fit$coefficients
    }

    start[fit$qr$pivot] <- fit$coefficients

    if (!(family$family %in% c("poisson", "binomial"))) {
        if (exists ("V.coefs") == FALSE) {
            if(NCOL(x.nobs)==1){
                V.coefs <- chol2inv(fit$qr$qr[1:nvars])
            }else{
                V.coefs <- chol2inv(fit$qr$qr[1:nvars, 1:nvars, drop = FALSE])
            }
        }
                                        #mse.resid <- mean((w * (z - x.nobs %*% fit$coefficients))^2) ## LOCAL VARIABLE
                                        #mse.resid <- mean ( (fit$y[1:nobs] - w * predictions)^2)
        ## mse.uncertainty <- mean(diag(x.nobs %*% V.coefs %*% t(x.nobs))) * state$dispersion
        mse.resid <- mean ( ((z.star*w.star)[1:sum(state$good)] - w * predictions[state$good,])^2)
        mse.uncertainty <- max (0, mean(rowSums(( x.nobs %*% V.coefs ) * x.nobs)) * state$dispersion) #faster  ## LOCAL VARIABLE
        state$dispersion <- mse.resid + mse.uncertainty
    }

    state$eta <- drop(predictions) + offset
    state$mu <- family$linkinv(state$eta)
    state$mu.eta.val <- family$mu.eta(state$eta)
    dev <- sum(family$dev.resids(y, state$mu, weights))

    if (!is.finite (dev)||(!is.finite(state$dispersion)) || !isTRUE(family$valideta(state$eta) && family$validmu(state$mu))) {
        if ((!is.finite(dev))|(!is.finite(state$dispersion))) {
            if (is.null(coefold)) {
                stop("no valid set of coefficients has been found: please supply starting values", call. = FALSE)
            }
            warning("step size truncated due to divergence", call. = FALSE)
        } else if (!isTRUE(family$valideta(state$eta) && family$validmu(state$mu))) {
            warning("step size truncated: out of bounds", call. = FALSE)
        }
        ii <- 1
        while (ii <= control$maxit & !is.finite (dev)) {
            ii <- ii + 1
            start <- (start + coefold) / 2
            state$eta <- if(NCOL(x.nobs)==1){
                             drop(x.nobs * start)
                         }else{
                             drop(x.nobs %*% start)
                         }
            state$mu <- family$linkinv(state$eta + offset)
            dev <- sum (family$dev.resids (y, state$mu, weights))
        }
        if (ii > control$maxit) {
            stop("inner loop 1; cannot correct step size")
        }

        ii <- 1
        while (ii <= control$maxit &  !isTRUE(family$valideta(state$eta) && family$validmu(state$mu))) {
            ii <- ii + 1
            start <- (start + coefold) / 2
            state$eta <- if(NCOL(x.nobs)==1){
                             drop(x.nobs * start)
                         }else{
                             drop(x.nobs %*% start)
                         }
            state$mu <- family$linkinv(state$eta + offset)
        }
        if (ii > control$maxit) {
            stop("inner loop 2; cannot correct step size")
        }


        state$boundary <- TRUE
        if (control$trace){
            cat("Step halved: new deviance =", dev, "\n")
        }
    }

    list (  Start = start,
          Coefold = coefold,
          eta=state$eta,
          mu=state$mu,
          mu.eta.val=state$mu.eta.val,
          varmu=family$variance(state$mu),
          good=(weights > 0) & (state$mu.eta.val != 0),
          dispersion=state$dispersion,
          dev=dev,
          fit=fit,
          conv=FALSE,
          boundary=state$boundary,
          prior.sd=state$prior.sd,
          z=z,
          w=w)
}


.bayesglm.fit.loop.initializePriors <- function (prior.mean, prior.scale, prior.df) {
    list(mean=prior.mean,
         scale=prior.scale,
         df=prior.df)

}

.bayesglm.fit.loop.print <- function (state, priors, family, print.unnormalized.log.posterior, intercept, y) {
    if (print.unnormalized.log.posterior && family$family == "binomial") {
        logprior <- sum( dt( state$fit$coefficients, priors$df , priors$mean, log = TRUE ) )
                                        #xb <- invlogit( x.nobs %*% coefs.hat )
        xb <- invlogit (state$eta - offset)
        loglikelihood <- sum( log( c( xb[ y == 1 ], 1 - xb[ y == 0 ] ) ) )
        cat( "log prior: ", logprior, ", log likelihood: ", loglikelihood, ",
                        unnormalized log posterior: ", loglikelihood +logprior, "\n" ,sep="")
    }
}


.bayesglm.fit.loop.createAuxillaryItems <- function (state, nvars, nobs, xnames, offset) {
    if (state$fit$rank < nvars) {
        state$fit$coefficients[seq(state$fit$rank + 1, nvars)] <- NA
    }
    residuals <- rep.int(NA, nobs)
    residuals[state$good] <- state$z[state$good] - (state$eta[state$good] - offset[state$good])

    state$fit$qr$qr <- as.matrix(state$fit$qr$qr)

    nr <- min(sum(state$good), nvars)
    if (nr < nvars) {
        Rmat <- diag(x = 0, nvars)
        Rmat[1:nr, 1:nvars] <- state$fit$qr$qr[1:nr, 1:nvars, drop=FALSE]
    }
    else{
        if(NCOL(state$fit$qr$qr)==1){
            Rmat <- state$fit$qr$qr[1:nvars]
        }else{
            Rmat <- state$fit$qr$qr[1:nvars, 1:nvars, drop=FALSE]
        }
    }
    Rmat <- as.matrix(Rmat)
    Rmat[row(Rmat) > col(Rmat)] <- 0
    names (state$fit$coefficients) <- xnames
    colnames(state$fit$qr$qr) <- xnames
    dimnames(Rmat) <- list(xnames, xnames)
    list (residuals=residuals,
          Rmat=Rmat,
          state=state)
}

.bayesglm.fit.loop.printWarnings <- function (Warning, state, family) {
    if(Warning){
        if (!state$conv){
            warning("algorithm did not converge")
        }
        if (state$boundary){
            warning("algorithm stopped at boundary value")
        }
        eps <- 10 * .Machine$double.eps
        if (family$family == "binomial") {
            if (any(state$mu > 1 - eps) || any(state$mu < eps)) {
                warning("fitted probabilities numerically 0 or 1 occurred")
            }
        }
        if (family$family == "poisson") {
            if (any(state$mu < eps)){
                warning("fitted rates numerically 0 occurred")
            }
        }
    }
}


.bayesglm.fit.loop.main.ideal <- function (control, x, y, nvars, nobs, weights, offset,
                                           intercept, scaled,
                                           start, etastart, mustart, dispersion,
                                           family,
                                           prior.mean, prior.scale, prior.df, 
                                           print.unnormalized.log.posterior,
                                           Warning) {

    xnames <- dimnames(x)[[2]]
    x.nobs <- x[1:nobs, ,drop=FALSE]

    .bayesglm.fit.loop.validateInputs (etastart, start, nvars, dimnames(x)[[2]])

                                        # invalid starting point (should be moved out of here):
                                        # is.null(start)==false & length(start) != nvars

    priors <- .bayesglm.fit.loop.initializePriors (prior.mean = prior.mean, prior.scale = prior.scale, prior.df = prior.df)

    state <- .bayesglm.fit.loop.initializeState (start, etastart, mustart, dispersion, offset, x.nobs, var(y), nvars, family, weights, prior.sd = priors$scale, y)
                                        #core elements of state:
                                        # good: derived from eta
                                        # eta
                                        # x, y, nvars, nobs: invariants

                                        #fortran.call.parameters <- .bayesglm.fit.loop.initialMemoryAllocation (control$epsilon, nvars, sum (state$good))

    for (iter in 1:control$maxit) {
        dispersionold <- state$dispersion
        devold <- state$dev
        state <- .bayesglm.fit.loop.updateState (state, priors, family, #fortran.call.parameters,
                                                 offset, weights,
                                                 y, x, x.nobs, nvars, nobs,
                                                 intercept, scaled, control)
        if (control$trace){
            cat("Deviance =", state$dev, "Iterations -", iter, "\n")
        }

        if (.bayesglm.fit.loop.validateState(state, family, control, iter, dispersionold, devold) == FALSE) {
            state$conv <- TRUE
            break
        }
        .bayesglm.fit.loop.print(state, priors, family, print.unnormalized.log.posterior, intercept, y)
    }

    .bayesglm.fit.loop.printWarnings(Warning, state, family)
    output <- .bayesglm.fit.loop.createAuxillaryItems(state, nvars, nobs, xnames, offset)
    ## residuals=residuals
    ## Rmat=Rmat
    ## state

    list (fit=output$state$fit,
          good=output$state$good,
          z=output$state$z,
          ## z.star=output$state$z.star,
          w=output$state$w,
          ## w.star=output$state$w.star,
          ngoodobs=sum (output$state$good),
          prior.scale=priors$scale,
          prior.sd=output$state$prior.sd,
          eta=output$state$eta,
          mu=output$state$mu,
          dev=output$state$dev,
          dispersion=output$state$dispersion,
          dispersionold=dispersionold,
          start=output$state$fit$coefficients,
          coefold=output$state$fit$coefficients,
          devold=devold,
          conv=output$state$conv,
          iter=iter,
          boundary=output$state$boundary,
          Rmat=output$Rmat,
          residuals=output$residuals)
}


.bayesglm.fit.cleanup <- function (ynames, residuals, mu, eta, nobs, weights, w, good, linkinv, dev.resids, y, intercept, fit, offset, EMPTY, dev, aic) {
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    wt <- rep.int(0, nobs)
    wt[good] <- w^2
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    wtdmu <- if (intercept){
                 sum(weights * y)/sum(weights)
             }
             else{
                 linkinv(offset)
             }
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(w == 0)
    nulldf <- n.ok - as.integer(intercept)

    rank <- if (EMPTY){
                0
            }
            else{
                fit$rank
            }
    resdf <- n.ok - rank
    aic.model <- aic(y, n.ok, mu, weights, dev) + 2 * rank

    list (residuals=residuals,
          mu=mu,
          eta=eta,
          wt = wt,
          weights = weights,
          y=y,
          wtdmu=wtdmu,
          nulldev=nulldev,
          n.ok=n.ok,
          nulldf=nulldf,
          rank=rank,
          resdf=resdf,
          aic.model=aic.model)
}


predict.bayesglm <- function (object, newdata = NULL, type = c("link", "response",
                                                               "terms"), se.fit = FALSE, dispersion = NULL, terms = NULL,
                              na.action = na.pass, ...)
{
    type <- match.arg(type)
    na.act <- object$na.action
    object$na.action <- NULL
    if (!se.fit) {
        if (missing(newdata)) {
            pred <- switch(type, link = object$linear.predictors,
                           response = object$fitted.values, terms = predictLM(object,
                                                                              se.fit = se.fit, scale = 1, type = "terms",
                                                                              terms = terms))
            if (!is.null(na.act))
                pred <- napredict(na.act, pred)
        }
        else {
            pred <- predictLM(object, newdata, se.fit, scale = 1,
                              type = ifelse(type == "link", "response", type),
                              terms = terms, na.action = na.action)
            switch(type, response = {
                pred <- family(object)$linkinv(pred)
            }, link = , terms = )
        }
    }
    else {
        if (inherits(object, "survreg"))
            dispersion <- 1
        if (is.null(dispersion) || dispersion == 0)
            dispersion <- summary(object, dispersion = dispersion)$dispersion
        residual.scale <- as.vector(sqrt(dispersion))
        pred <- predictLM(object, newdata, se.fit, scale = residual.scale,
                          type = ifelse(type == "link", "response", type),
                          terms = terms, na.action = na.action)
        fit <- pred$fit
        se.fit <- pred$se.fit
        switch(type, response = {
            se.fit <- se.fit * abs(family(object)$mu.eta(fit))
            fit <- family(object)$linkinv(fit)
        }, link = , terms = )
        if (missing(newdata) && !is.null(na.act)) {
            fit <- napredict(na.act, fit)
            se.fit <- napredict(na.act, se.fit)
        }
        pred <- list(fit = fit, se.fit = se.fit, residual.scale = residual.scale)
    }
    pred
}

predictLM <- function (object, newdata, se.fit = FALSE, scale = NULL, df = Inf,
                       interval = c("none", "confidence", "prediction"), level = 0.95,
                       type = c("response", "terms"), terms = NULL, na.action = na.pass,
                       pred.var = res.var/weights, weights = 1, ...)
{
    tt <- terms(object)
    keep.order <- object$keep.order
    drop.baseline <- object$drop.baseline
    if (!inherits(object, "lm"))
        warning("calling predict.lm(<fake-lm-object>) ...")
    if (missing(newdata) || is.null(newdata)) {
        mm <- X <- model.matrix(object)
        mmDone <- TRUE
        offset <- object$offset
    }
    else {
        Terms <- delete.response(tt)
        m <- model.frame(Terms, newdata, na.action = na.action,
                         xlev = object$xlevels)
        if (!is.null(cl <- attr(Terms, "dataClasses")))
            .checkMFClasses(cl, m)
        X <- model.matrixBayes(Terms, m, contrasts.arg = object$contrasts, keep.order = keep.order, drop.baseline = drop.baseline)
        offset <- rep(0, nrow(X))
        if (!is.null(off.num <- attr(tt, "offset")))
            for (i in off.num) offset <- offset + eval(attr(tt,
                                                            "variables")[[i + 1]], newdata)
        if (!is.null(object$call$offset))
            offset <- offset + eval(object$call$offset, newdata)
        mmDone <- FALSE
    }
    n <- length(object$residuals)
    p <- object$rank
    p1 <- seq_len(p)
    piv <- if (p)
               getQr(object)$pivot[p1]
    if (p < ncol(X) && !(missing(newdata) || is.null(newdata)))
        warning("prediction from a rank-deficient fit may be misleading")
    beta <- object$coefficients
    predictor <- drop(X[, piv, drop = FALSE] %*% beta[piv])
    if (!is.null(offset))
        predictor <- predictor + offset
    interval <- match.arg(interval)
    if (interval == "prediction") {
        if (missing(newdata))
            warning("Predictions on current data refer to _future_ responses\n")
        if (missing(newdata) && missing(weights)) {
            w <- .weights.default(object)
            if (!is.null(w)) {
                weights <- w
                warning("Assuming prediction variance inversely proportional to weights used for fitting\n")
            }
        }
        if (!missing(newdata) && missing(weights) && !is.null(object$weights) &&
            missing(pred.var))
            warning("Assuming constant prediction variance even though model fit is weighted\n")
        if (inherits(weights, "formula")) {
            if (length(weights) != 2L)
                stop("'weights' as formula should be one-sided")
            d <- if (missing(newdata) || is.null(newdata))
                     model.frame(object)
                 else newdata
            weights <- eval(weights[[2L]], d, environment(weights))
        }
    }
    type <- match.arg(type)
    if (se.fit || interval != "none") {
        res.var <- if (is.null(scale)) {
                       r <- object$residuals
                       w <- object$weights
                       rss <- sum(if (is.null(w)) r^2 else r^2 * w)
                       df <- object$df.residual
                       rss/df
                   }
                   else scale^2
        if (type != "terms") {
            if (p > 0) {
                XRinv <- if (missing(newdata) && is.null(w))
                             qr.Q(getQr(object))[, p1, drop = FALSE]
                         else X[, piv] %*% qr.solve(qr.R(getQr(object))[p1, p1])
                ip <- drop(XRinv^2 %*% rep(res.var, p))
            }
            else ip <- rep(0, n)
        }
    }
    if (type == "terms") {
        if (!mmDone) {
            mm <- model.matrixBayes(object, keep.order = keep.order, drop.baseline = drop.baseline)
            mmDone <- TRUE
        }
        aa <- attr(mm, "assign")
        ll <- attr(tt, "term.labels")
        hasintercept <- attr(tt, "intercept") > 0L
        if (hasintercept)
            ll <- c("(Intercept)", ll)
        aaa <- factor(aa, labels = ll)
        asgn <- split(order(aa), aaa)
        if (hasintercept) {
            asgn$"(Intercept)" <- NULL
            if (!mmDone) {
                mm <- model.matrixBayes(object, keep.order = keep.order, drop.baseline = drop.baseline)
                mmDone <- TRUE
            }
            avx <- colMeans(mm)
            termsconst <- sum(avx[piv] * beta[piv])
        }
        nterms <- length(asgn)
        if (nterms > 0) {
            predictor <- matrix(ncol = nterms, nrow = NROW(X))
            dimnames(predictor) <- list(rownames(X), names(asgn))
            if (se.fit || interval != "none") {
                ip <- matrix(ncol = nterms, nrow = NROW(X))
                dimnames(ip) <- list(rownames(X), names(asgn))
                Rinv <- qr.solve(qr.R(getQr(object))[p1, p1])
            }
            if (hasintercept)
                X <- sweep(X, 2L, avx, check.margin = FALSE)
            unpiv <- rep.int(0L, NCOL(X))
            unpiv[piv] <- p1
            for (i in seq.int(1L, nterms, length.out = nterms)) {
                iipiv <- asgn[[i]]
                ii <- unpiv[iipiv]
                iipiv[ii == 0L] <- 0L
                predictor[, i] <- if (any(iipiv > 0L))
                                      X[, iipiv, drop = FALSE] %*% beta[iipiv]
                                  else 0
                if (se.fit || interval != "none")
                    ip[, i] <- if (any(iipiv > 0L))
                                   as.matrix(X[, iipiv, drop = FALSE] %*% Rinv[ii,
                                                                             , drop = FALSE])^2 %*% rep.int(res.var,
                                                                                                            p)
                else 0
            }
            if (!is.null(terms)) {
                predictor <- predictor[, terms, drop = FALSE]
                if (se.fit)
                    ip <- ip[, terms, drop = FALSE]
            }
        }
        else {
            predictor <- ip <- matrix(0, n, 0L)
        }
        attr(predictor, "constant") <- if (hasintercept)
                                           termsconst
                                       else 0
    }
    if (interval != "none") {
        tfrac <- qt((1 - level)/2, df)
        hwid <- tfrac * switch(interval, confidence = sqrt(ip),
                               prediction = sqrt(ip + pred.var))
        if (type != "terms") {
            predictor <- cbind(predictor, predictor + hwid %o%
                                          c(1, -1))
            colnames(predictor) <- c("fit", "lwr", "upr")
        }
        else {
            if (!is.null(terms))
                hwid <- hwid[, terms, drop = FALSE]
            lwr <- predictor + hwid
            upr <- predictor - hwid
        }
    }
    if (se.fit || interval != "none") {
        se <- sqrt(ip)
        if (type == "terms" && !is.null(terms) && !se.fit)
            se <- se[, terms, drop = FALSE]
    }
    if (missing(newdata) && !is.null(na.act <- object$na.action)) {
        predictor <- napredict(na.act, predictor)
        if (se.fit)
            se <- napredict(na.act, se)
    }
    if (type == "terms" && interval != "none") {
        if (missing(newdata) && !is.null(na.act)) {
            lwr <- napredict(na.act, lwr)
            upr <- napredict(na.act, upr)
        }
        list(fit = predictor, se.fit = se, lwr = lwr, upr = upr,
             df = df, residual.scale = sqrt(res.var))
    }
    else if (se.fit)
        list(fit = predictor, se.fit = se, df = df, residual.scale = sqrt(res.var))
    else predictor
}



getQr <- function(x, ...){
    if (is.null(r <- x$qr))
        stop("lm object does not have a proper 'qr' component.\n Rank zero or should not have used lm(.., qr=FALSE).")
    r
}

logit <- function (x) {
    log(x/(1-x))
}

##' Inverse of logistic transformation
##' 
##' @export
##' @param x numeric
##' @return numeric
##' @examples
##' x <- 1:5
##' invlogit(log(x/(1-x)))
invlogit <- function (x) {
    1/(1+exp(-x))
}

model.matrixBayes <- function(object, data = environment(object),
                              contrasts.arg = NULL, xlev = NULL, keep.order=FALSE, drop.baseline=FALSE,...)
{
                                        #class(object) <- c("terms", "formula")
    t <- if( missing( data ) ) { 
             terms( object ) 
         }else{ 
             terms.formula(object, data = data, keep.order=keep.order) 
         }
    attr(t, "intercept") <- attr(object, "intercept")
    if (is.null(attr(data, "terms"))){ 
        data <- model.frame(object, data, xlev=xlev) 
    }else {
        reorder <- match(sapply(attr(t,"variables"), deparse, width.cutoff=500)[-1], names(data))
        if (any(is.na(reorder))) {
            stop( "model frame and formula mismatch in model.matrix()" ) 
        }
        if(!identical(reorder, seq_len(ncol(data)))) {
            data <- data[,reorder, drop = FALSE] 
        }
    }
    int <- attr(t, "response")
    if(length(data)) {      # otherwise no rhs terms, so skip all this
        
        if (drop.baseline){
            contr.funs <- as.character(getOption("contrasts"))
        }else{
            contr.funs <- as.character(list("contr.bayes.unordered", "contr.bayes.ordered"))
        }
        
        namD <- names(data)
        ## turn any character columns into factors
        for(i in namD)
            if(is.character( data[[i]] ) ) {
                data[[i]] <- factor(data[[i]])
                warning( gettextf( "variable '%s' converted to a factor", i ), domain = NA)
            }
        isF <- vapply(data, function(x) is.factor(x) || is.logical(x), NA)        
        isF[int] <- FALSE
        isOF <- vapply(data, is.ordered, NA)
        for( nn in namD[isF] )            # drop response
            if( is.null( attr( data[[nn]], "contrasts" ) ) ) {
                contrasts( data[[nn]] ) <- contr.funs[1 + isOF[nn]]
            }
        ## it might be safer to have numerical contrasts:
        ##    get(contr.funs[1 + isOF[nn]])(nlevels(data[[nn]]))
        if ( !is.null( contrasts.arg ) && is.list( contrasts.arg ) ) {
            if ( is.null( namC <- names( contrasts.arg ) ) ) {
                stop( "invalid 'contrasts.arg' argument" )
            }
            for (nn in namC) {
                if ( is.na( ni <- match( nn, namD ) ) ) {
                    warning( gettextf( "variable '%s' is absent, its contrast will be ignored", nn ), domain = NA )
                }
                else {
                    ca <- contrasts.arg[[nn]]
                    if( is.matrix( ca ) ) {
                        contrasts( data[[ni]], ncol( ca ) ) <- ca
                    }
                    else { 
                        contrasts( data[[ni]] ) <- contrasts.arg[[nn]]
                    }
                }
            }
        }
    } else {               # internal model.matrix needs some variable
        isF  <-  FALSE
        data <- data.frame(x=rep(0, nrow(data)))
    }
                                        #ans  <- .Internal( model.matrix( t, data ) )
    ans  <- model.matrix.default(object=t, data=data)
    cons <- if(any(isF)){
                lapply( data[isF], function(x) attr( x,  "contrasts") ) 
            }else { NULL }
    attr(ans, "contrasts" ) <- cons
    ans
}

.weights.default <- function (object, ...) 
{
    wts <- object$weights
    if (is.null(wts)) 
        wts
    else napredict(object$na.action, wts)
}
