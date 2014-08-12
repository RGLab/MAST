## Likelihood functions and other helpers for shrunken dispersion estimates for zlm

## rNg: residual Ng: Ng -p, where p is the dimension of the model
## SSg: residual sum of squares
getMarginalHyperLikelihood <- function(rNg, SSg, deriv=FALSE){
    if(!deriv){
        fun <- function(theta){
            stopifnot(names(theta)==c('a0', 'b0'))
            a0 <- theta['a0']
            b0 <- theta['b0']

            Li <- -lbeta(rNg/2, a0)-rNg/2*log(b0)-log(1+SSg/(2*b0))*(rNg/2+a0)
            return(sum(Li))
        }
    } else{
        fun <- function(theta){
            stopifnot(names(theta)==c('a0', 'b0'))
            a0 <- theta['a0']
            b0 <- theta['b0']
            score_a0_i <- digamma(rNg/2+a0)-digamma(a0)-log(1+SSg/(2*b0))
            score_b0_i <- (a0*SSg-rNg*b0)/(SSg*b0+2*b0^2)
            return(c(a0=sum(score_a0_i), b0=sum(score_b0_i)))
        }
    }
    fun
}

## probably need a global optimization routine--plus there are multiple roots potentially.
## or just a good starting value
solveMoM <- function(rNg, SSg){
    rbar <- mean(SSg/rNg)
    rbarbar <- mean(SSg^2/(rNg*(rNg+2)))
    a0mom <- function(a0) (2*(a0-1)^2*rbar^2  -rbarbar^2*((a0-2)*(a0-4)))^2
    
    a0slv <- optimize(a0mom, c(0, 10))
    a0 <- a0slv$minimum
    b0 <- (a0-1)*rbar
    c(a0, b0)
}

##' @importFrom plyr aaply
getSSg_rNg <- function(sca, mm){
    aaply(exprs(sca), 2, function(y){
            SSg <- NA
            rNg <- NA
            try({
                pos <- y>0
                yp <- y[pos]
                mp <- mm[pos,]
                QR <- qr(mp)
                resid <- qr.resid(QR, yp)
                SSg <- crossprod(resid)
                rNg <- length(yp)-QR$rank
                   }, silent=TRUE)
            return(c(SSg=SSg, rNg=rNg))
        })
}

##' Estimate hyperparameters for hierarchical variance model for continuous component
##'
##' \code{ebayesControl} is a named list with (optional) components 'method' (one of 'MOM' or 'MLE') and 'model' (one of 'H0' or 'H1')
##' method MOM uses a method-of-moments estimator, while MLE using the marginal likelihood.
##' H0 model estimates the precisions using the intercept alone in each gene, while H1 fits the full model specified by \code{formula}
##' @param sca \code{SingleCellAssay}
##' @param ebayesControl list with (optional) components 'method', 'model'.  See details.
##' @param Formula a formula (using variables in \code{cData(sca)} used when \code{model='H1'}.
##' @param truncate Genes with sample precisions exceeding this value are discarded when estimating the hyper parameters
##' @return \code{numeric} of length two, giving the hyperparameters in terms of a variance (\code{v}) and prior observations (\code{df}), inside a \code{structure}, with component \code{hess}, giving the Fisher Information of the hyperparameters.
ebayes <- function(sca, ebayesControl, Formula, truncate=Inf){
     ## Empirical bayes method
    defaultCtl <- list(method='MLE', model='H0')
    if (is.null(ebayesControl)){
    ebayesControl <- list()
  }
    missingControl <- setdiff(names(ebayesControl), names(ebayesControl))
    ebayesControl[missingControl] <- defaultCtl[missingControl]
    method <- match.arg(ebayesControl[['method']], c('MOM', 'MLE'))
    model <- match.arg(ebayesControl[['model']], c('H0', 'H1'))

    ee <- exprs(sca)
    ee[ee==0] <- NA
    
    if(model == 'H0'){
        ee <- scale(ee, scale=FALSE, center=TRUE)
        ## Global variance
        rNg <- colSums(!is.na(ee), na.rm=TRUE)-1
        SSg <- colSums(ee^2, na.rm=TRUE)
        valid <- rNg>0 & rNg/SSg < truncate
        rNg <- rNg[valid]
        SSg <- SSg[valid]
    } else if(model == 'H1'){
        mm <- model.matrix(Formula, cData(sca))

        allfits <- getSSg_rNg(sca, mm)
        valid <- apply(!is.na(allfits), 1, all) & allfits[, 'rNg']/allfits[, 'SSg']<truncate
        valid[is.na(valid)] <- FALSE
        SSg <- allfits[valid,'SSg']
        rNg <- allfits[valid, 'rNg']
    }

    if(method == 'MLE'){
        fn <- getMarginalHyperLikelihood(rNg, SSg, deriv=FALSE)
        grad <- getMarginalHyperLikelihood(rNg, SSg, deriv=TRUE)
        O <- optim(c(a0=1, b0=1), fn, gr=grad, method='L-BFGS', lower=.001, upper=Inf, control=list(fnscale=-1), hessian=TRUE)
        if(O$convergence!=0) warning('Hyper parameter estimation might have failed', O$message)
        #O <- optim(c(a0=1, b0=1), fn, method='L-BFGS', lower=.001, upper=Inf, control=list(fnscale=-1))
        th <- O$par
    } else if(method == 'MOM'){
        th <- solveMoM(rNg, SSg)
        O <- list(hessian=NA)
    }

    v <- max(th['b0']/th['a0'], 0)
    df <- max(2*th['a0'], 0)
    structure(c(v=v, df=df), hess=O$hessian)
}
