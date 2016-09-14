.expit <- function(eta) exp(eta)/(1+exp(eta))
## derivative of expit w/r/t beta
.dexpit <- function(eta) exp(eta)/(1+exp(eta))^2

## dot product of contrast and Coef
safeContrastDP <- function(contrast, Coef) uncomplexify(complexifyNA(Coef) %*% t(contrast))

## quadratic form of `contrast` about `vc` to give the covariance of the transformation defined by `contrast`
safeContrastQF <- function(contrast, vc) uncomplexify(tcrossprod(contrast, complexifyNA(vc)) %*% t(contrast))

##' Calculate log-fold changes from hurdle model components
##'
##' Using the delta method, estimate the log-fold change from a state given by a vector contrast0 and the state(s) given by contrast1.
##'
##' The log-fold change is defined as follows.  For each gene, let \eqn{u(x)} be the expected value of the continuous component, given a covariate x and the estimated coefficients \code{coefC}, ie, \eqn{u(x)=} \code{crossprod(x, coefC)}.
##' Likewise, Let \eqn{v(x)=} \code{1/(1+exp(-crossprod(coefD, x)))} be the expected value of the discrete component.
##'  The log fold change from contrast0 to contrast1 is defined as
##' \deqn{u(contrast1)v(contrast1)-u(contrast0)v(contrast0).}
##' Note that for this to be a log-fold change, then the regression for u must have been fit on the log scale.  This is returned in the matrix \code{logFC}.
##' An approximation of the variance of \code{logFC} (applying the delta method to  formula defined above) is provided in \code{varLogFC}.
##' @param zlmfit ZlmFit output
##' @param contrast0 vector of coefficients giving baseline contrast, or a \code{\link{Hypothesis}}.  If missing, then the '(Intercept)' is used as baseline.
##' @param contrast1 matrix of coefficients giving comparison contrasts, or a \code{\link{Hypothesis}}.  If missing, then all non-(Intercept) coefficients are compared.
##' @seealso Hypothesis
##' @return list of matrices `logFC` and `varLogFC`, giving the log-fold-changes for each contrast (columns) and genes (rows) and the estimated sampling variance thereof
##' @examples
##' data(vbetaFA)
##' zz <- zlm.SingleCellAssay( ~ Stim.Condition+Population, vbetaFA[1:5,])
##' ##log-fold changes in terms of intercept (which is Stim(SEB) and CD154+VbetaResponsive)
##' lfcStim <- logFC(zz)
##' ##If we want to compare against unstim, we can try the following
##' coefnames <- colnames(coef(zz, 'D'))
##' contrast0 <- setNames(rep(0, length(coefnames)), coefnames)
##' contrast0[c('(Intercept)', 'Stim.ConditionUnstim')] <- 1
##' contrast1 <- diag(length(coefnames))
##' rownames(contrast1)<-colnames(contrast1)<-coefnames
##' contrast1['(Intercept)',]<-1
##' lfcUnstim <- logFC(zz, contrast0, contrast1)
##' ##log-fold change with itself is 0
##' stopifnot(all(lfcUnstim$logFC[,2]==0))
##' ##inverse of log-fold change with Stim as reference
##' stopifnot(all(lfcStim$logFC[,1]==(-lfcUnstim$logFC[,1])))
##' ##As a data.table:
##' getLogFC(zz)
##' @export
logFC <- function(zlmfit, contrast0, contrast1){
    coname <- colnames(coef(zlmfit, 'D'))
    genes <- rownames(coef(zlmfit, 'D'))
    if(missing(contrast0)){
        if(! '(Intercept)' %in% coname) stop("Assuming comparision to intercept, but I can't figure out what coefficient that corresponds to.  Provide `contrast0`.")
        contrast0 <- t(as.matrix(setNames(rep(0, length(coname)), coname)))
        contrast0[,'(Intercept)'] <- 1
        rownames(contrast0) <- '(Intercept)'
    } else if(inherits(contrast0, 'Hypothesis')){
        gh <- generateHypothesis(contrast0, coname)
        contrast0 <- gh@contrastMatrix
    } 
    if(missing(contrast1)){
        contrast1 <- cbind(0, diag(1, nrow=length(coname)-1))
        colnames(contrast1) <- coname
        rownames(contrast1) <- setdiff(coname, rownames(contrast0))
        ## add contrast0
        contrast1 <- t(contrast1)+as.numeric(contrast0)
    } else if(inherits(contrast1, 'Hypothesis')){
        gh <- generateHypothesis(contrast1, coname)
        contrast1 <- gh@contrastMatrix
    }
    contrast1 <- t(contrast1)
    Contr <- rbind(contrast0, contrast1)
    ## Get expectation and variance of contrasts, discrete and continuous
    ## expectation of contrast.  genes x contrast
    mu.cont <- safeContrastDP(Contr, coef(zlmfit, 'C'))
    eta.disc <- safeContrastDP(Contr, coef(zlmfit, 'D'))
    ## variance of contrasts.  genes x contrast x contrast
    vcont <- aaply(vcov(zlmfit, 'C'), 3, safeContrastQF, contrast=Contr)
    vcd <- vcov(zlmfit, 'D')
    ## variance of contrast _times_ jacobian due to expit transformation
    vdisc <- aaply(seq_along(genes), 1, function(i){
        ## variance-covariance of discrete beta
        vcc <- safeContrastQF(Contr, vcd[,,i])
        ## gradient of expit function. it's actually a diagonal matrix.
        jacobian <- .dexpit(eta.disc[i,])
        ## same as pre and post multiplying by diagonal jacobian
        vcc*tcrossprod(jacobian)
    })
    mu.disc <- .expit(eta.disc)

    ## expectation of all contrasts, product of discrete and continuous
    mu.prod <- mu.cont*mu.disc
    
    ## variance of the above
    ## diagonals of covariance matrices
    diagIdx <- as.matrix(expand.grid(i=seq_along(genes), j=seq_len(nrow(Contr))))
    diagIdx <- cbind(diagIdx, k=diagIdx[,'j'])
    dvcont <- matrix(vcont[diagIdx], nrow=length(genes))
    dvdisc <- matrix(vdisc[diagIdx], nrow=length(genes))
    ## variance of product is product of variance + higher order terms 
    v.prod <- dvcont*dvdisc + mu.cont^2*dvdisc + mu.disc^2*dvcont
    
    ## covariance between contrast1 and contrast0
    covc1c0Idx <- diagIdx[-(seq_along(genes)),] #drop j=k=1 indices, ie, 1,1 entry of covariances
    covc1c0Idx[,'k'] <- 1
    c1c0cont <- matrix(vcont[covc1c0Idx], nrow=length(genes))
    c1c0disc <- matrix(vdisc[covc1c0Idx], nrow=length(genes))
    ## covariance of product is product of covariance + higher order terms
    cov.prod <- c1c0cont*c1c0disc+genewiseMult(mu.cont[,1], mu.cont[,-1,drop=FALSE])*c1c0disc + genewiseMult(mu.disc[,1], mu.disc[,-1,drop=FALSE])*c1c0cont

    ## _differences_ between expectation of baseline contrast and all others
    ## baseline is in row 1
    lfc <-mu.prod[,-1,drop=FALSE]-mu.prod[,1,drop=TRUE]
    ## var of difference is sum of variances minus 2*covariance
    vlfc <- v.prod[,-1,drop=FALSE]+v.prod[,1,drop=TRUE]-2*cov.prod
    dimnames(lfc) <- dimnames(vlfc) <- list(primerid=genes, contrast=rownames(contrast1))
    list(logFC=lfc,varLogFC=vlfc)
}

## attempting to be careful with a sneaky use of recycling of rowvec
genewiseMult <- function(rowvec, rowMajorMatrix){
    res <- rowvec*rowMajorMatrix
    ## debugging assertion
    ## stopifnot(all.equal(res[3,], rowvec[3]*rowMajorMatrix[3,]))
    res
}

if(getRversion() >= "2.15.1") globalVariables(c(
                                  'primerid',
                                  'z', 
                                  'varLogFC')) #getLogFC

##' @describeIn logFC Return results as a perhaps friendlier \code{data.table}
##' @export
getLogFC <- function(zlmfit, contrast0, contrast1){
    lfc <- logFC(zlmfit, contrast0, contrast1)
    logFC <- dcast.data.table(data.table(melt(lfc)), primerid + contrast ~ L1, fun.aggregate=mean, na.rm=TRUE)
    logFC[,primerid:=as.character(primerid)]
    logFC[,z:=logFC/sqrt(varLogFC)]
    setkey(logFC,primerid)
    logFC <- logFC[dimnames(lfc[[1]])$primerid,]
    logFC
}
