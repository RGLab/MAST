.expit <- function(eta) exp(eta)/(1+exp(eta))
## derivative of expit w/r/t beta
.dexpit <- function(eta) exp(eta)/(1+exp(eta))^2

## dot product of contrast and Coef
safeContrastDP <- function(contrast, Coef) uncomplexify(complexifyNA(Coef) %*% t(contrast))

## quadratic form of `contrast` about `vc` to give the covariance of the transformation defined by `contrast`
safeContrastQF <- function(contrast, vc) uncomplexify(tcrossprod(contrast, complexifyNA(vc)) %*% t(contrast))

##' Calculate log-fold changes from hurdle model components
##'
##' Using the delta method, estimate the log-fold change from a state given by a vector contrast0 and the state(s) given by contrast1
##' @param zlmfit ZlmFit output
##' @param contrast0 vector of coefficients giving baseline contrast.  If missing, then the '(Intercept)' is used as baseline.
##' @param contrast1 matrix of coefficients giving comparison contrasts.  If missing, then all non-(Intercept) coefficients are compared.
##' @return list of matrices `logFC` and `varLogFC`, giving the log-fold-changes for each contrast (columns) and genes (rows) and the estimated sampling variance thereof
##' @export
logFC <- function(zlmfit, contrast0, contrast1){
    coname <- colnames(coef(zlmfit, 'D'))
    genes <- rownames(coef(zlmfit, 'D'))
    if(missing(contrast0)){
        if(! '(Intercept)' %in% coname) stop("Assuming comparision to intercept, but I can't figure out what coefficient that corresponds to.  Provide `contrast0`.")
        contrast0 <- t(as.matrix(setNames(rep(0, length(coname)), coname)))
        contrast0[,'(Intercept)'] <- 1
        rownames(contrast0) <- '(Intercept)'
    }
    if(missing(contrast1)){
        contrast1 <- cbind(0, diag(1, nrow=length(coname)-1))
        colnames(contrast1) <- coname
        rownames(contrast1) <- setdiff(coname, rownames(contrast0))
        ## add contrast0
        contrast1 <- t(t(contrast1)+as.numeric(contrast0))
    } else{
        contrast1 <- t(contrast1)
    }
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
    ## debugging assertion
    ## stopifnot(all.equal(dvcont[1,], diag(vcont[1,,]), check.attributes=FALSE))
    ## stopifnot(all.equal(dvdisc[10,], diag(vdisc[10,,]), check.attributes=FALSE))
    ## variance of product is product of variance + higher order terms 
    v.prod <- dvcont*dvdisc + mu.cont^2*dvdisc + mu.disc^2*dvcont
    
    ## covariance between contrast1 and contrast0
    covc1c0Idx <- diagIdx[-(seq_along(genes)),] #drop j=k=1 indices, ie, 1,1 entry of covariances
    covc1c0Idx[,'k'] <- 1
    c1c0cont <- matrix(vcont[covc1c0Idx], nrow=length(genes))
    c1c0disc <- matrix(vdisc[covc1c0Idx], nrow=length(genes))
    ## debugging assertion
    ## stopifnot(all.equal(c1c0cont[1,], vcont[1,2:nrow(Contr),1], check.attributes=FALSE))
    ## stopifnot(all.equal(c1c0disc[10,], vdisc[10,2:nrow(Contr),1], check.attributes=FALSE))
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

##' @import data.table
##' @importFrom reshape melt
##' @describeIn logFC Return results as a perhaps friendlier \code{data.table}
##' @export
getLogFC <- function(zlmfit, contrast0, contrast1){
    lfc <- logFC(zlmfit, contrast0, contrast1)
    logFC <- dcast.data.table(data.table(melt(lfc)), primerid + contrast ~ L1)
    logFC[,primerid:=as.character(primerid)]
    logFC[,z:=logFC/sqrt(varLogFC)]
    setkey(logFC,primerid)
    logFC=logFC[dimnames(lfc[[1]])$primerid,]
    logFC
}
