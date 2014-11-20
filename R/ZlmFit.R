## static constructor of summary components
initSummaries <- function(genes, coefNames){
    p <- length(coefNames)
    ng <- length(genes)
    coefD <- coefC <- matrix(NA, nrow=ng, ncol=p, dimnames=list(primerid=genes, coef=coefNames))
    vcovC <- vcovD <- array(NA, dim=c(ng, p, p), dimnames=list(primerid=genes, coefNames, coefNames))
converged <- dispersion <- df.resid <- df.null <- deviance <- matrix(rep(NA, ng*2), nrow=ng, ncol=2, dimnames=list(genes, c('C', 'D')))
    as.environment(list(coefC=coefC, vcovC=vcovC, df.resid=df.resid, df.null=df.null, deviance=deviance, dispersion=dispersion, coefD=coefD, vcovD=vcov, converged=converged))
}

Glue <- function(...) abind(..., rev.along=0)

collectSummaries <- function(listOfSummaries){
    summaries <- list()
    summaries[['coefC']] <- do.call(rbind, lapply(listOfSummaries, '[[', 'coefC'))
    summaries[['vcovC']] <- do.call(Glue, lapply(listOfSummaries, '[[', 'vcovC'))
    summaries[['df.resid']] <- do.call(rbind, lapply(listOfSummaries, '[[', 'df.resid'))
    summaries[['df.null']] <- do.call(rbind, lapply(listOfSummaries, '[[', 'df.null'))
    summaries[['deviance']] <- do.call(rbind, lapply(listOfSummaries, '[[', 'deviance'))
    summaries[['dispersion']] <- do.call(rbind, lapply(listOfSummaries, '[[', 'dispersion'))
summaries[['dispersionNoshrink']] <- do.call(rbind, lapply(listOfSummaries, '[[', 'dispersionNoshrink'))
    summaries[['converged']] <- do.call(rbind, lapply(listOfSummaries, '[[', 'converged'))
    summaries[['loglik']] <- do.call(rbind, lapply(listOfSummaries, '[[', 'loglik'))
    summaries[['coefD']] <- do.call(rbind, lapply(listOfSummaries, '[[', 'coefD'))
    summaries[['vcovD']] <- do.call(Glue, lapply(listOfSummaries, '[[', 'vcovD'))

    summaries
}

## Refit using new MM and calculate stats
## hString is human-readable hypothesis
.lrtZlmFit <- function(zlmfit, newMM, hString){
    o1 <- zlmfit
    LMlike <- o1@LMlike
    model.matrix(LMlike) <- newMM
    message('Refitting on reduced model...')
    o0 <- zlm.SingleCellAssay(sca=o1@sca, LMlike=LMlike)
    lambda <- -2*(o0@loglik-o1@loglik)
    lambda <- ifelse(o0@converged & o1@converged, lambda, 0)
    df <- o0@df.resid-o1@df.resid
    df <- ifelse(o0@converged & o1@converged, df, 0)
    cst <- makeChiSqTable(as.data.frame(lambda), as.data.frame(df), hString)
    dimnames(cst) <- list(primerid=fData(zlmfit@sca)$primerid, test.type=dimnames(cst)[[2]], metric=dimnames(cst)[[3]])
    cst
}

##' Likelihood ratio test
##'
##' A 3D array with first dimension being the genes,
##' next dimension giving information about the test
##' (the degrees of freedom, Chisq statistic, and P value), and final dimension
##' being the value of these quantities on the
##' discrete, continuous and hurdle (combined) levels.
##' @param object ZlmFit
##' @param hypothesis See Details
##' @return 3D array
setMethod('lrTest',  signature=c(object='ZlmFit', hypothesis='character'), function(object, hypothesis){
    o1 <- object
    LMlike <- o1@LMlike
    F <- update.formula(LMlike@formula, formula(sprintf(' ~. - %s', hypothesis)))
    if(F==LMlike@formula) stop('Removing term ', sQuote(hypothesis), " doesn't actually alter the model, maybe due to marginality? Try specifying individual coefficents as a `CoefficientHypothesis`.")
    LMlike <- update(LMlike, F)
    .lrtZlmFit(o1, LMlike@modelMatrix, hypothesis)
})

setMethod('lrTest', signature=c(object='ZlmFit', hypothesis='CoefficientHypothesis'), function(object, hypothesis){
    h <- generateHypothesis(hypothesis, colnames(object@coefD))
    testIdx <- h@transformed
    newMM <- model.matrix(object@LMlike)[,-testIdx, drop=FALSE]
    .lrtZlmFit(object, newMM, hypothesis@.Data)
})

setMethod('lrTest', signature=c(object='ZlmFit', hypothesis='Hypothesis'), function(object, hypothesis){
    ## original fit
    h <- generateHypothesis(hypothesis, colnames(object@coefD))
    ## call using coefficient matrix
    lrTest(object, h@transformed)
})

## contrast matrices
setMethod('lrTest', signature=c(object='ZlmFit', hypothesis='matrix'), function(object, hypothesis){
    ## original fit
    LMlike <- object@LMlike
    MM <- .rotateMM(LMlike, hypothesis)
    testIdx <- attr(MM, 'testIdx')
    .lrtZlmFit(object, MM[,-testIdx, drop=FALSE], 'Contrast Matrix')
})

##' Wald test
##'
##' A 3D array with first dimension being the genes,
##' next dimension giving information about the test
##' (the degrees of freedom, Chisq statistic, and P value), and final dimension
##' being the value of these quantities on the
##' discrete, continuous and hurdle (combined) levels.
##' @param object ZlmFit
##' @param hypothesis See Details
##' @return 3D array
setMethod('waldTest',  signature=c(object='ZlmFit', hypothesis='matrix'), function(object, hypothesis){
    coefC <- coef(object, 'C')
    coefD <- coef(object, 'D')
    vcovC <- vcov(object, 'C')
    vcovD <- vcov(object, 'D')
    converged <- object@converged
    genes <- rownames(coefC)
    tests <- aaply(seq_along(genes), 1, function(i){
         .waldTest(coefC[i,],
                coefD[i,],
                vcovC[,,i],
                vcovD[,,i],
                hypothesis, converged[i,])
    }, .drop=FALSE)
    dimnames(tests)[[1]] <- genes
    names(dimnames(tests)) <- c('primerid', 'test.type', 'metric')
    tests
})

setMethod('waldTest',  signature=c(object='ZlmFit', hypothesis='CoefficientHypothesis'), function(object, hypothesis){
    cm <- .makeContrastMatrixFromCoefficientHypothesis(hypothesis,colnames(object@coefC))
    waldTest(object, cm)
})

setMethod('waldTest',  signature=c(object='ZlmFit', hypothesis='Hypothesis'), function(object, hypothesis){
    h <- generateHypothesis(hypothesis, colnames(object@coefD))
    waldTest(object, h@transformed)
})


setMethod('show', signature=c(object='ZlmFit'), function(object){
    cat('Fitted zlm on ', ncol(object@sca), ' genes and ', nrow(object@sca), ' cells.\n Using ', class(object@LMlike), ' to fit.\n')
})

setMethod('coef', signature=c(object='ZlmFit'), function(object, which, ...){
    which <- match.arg(which, c('C', 'D'))
    if(which=='C') object@coefC else object@coefD
})

setMethod('vcov', signature=c(object='ZlmFit'), function(object, which, ...){
    which <- match.arg(which, c('C', 'D'))
    if(which=='C') object@vcovC else object@vcovD
})

##' Standard error on coefficients
##'
##' @param object ZlmFit
##' @param which  character vector, one of "C" (continuous) or "D" (discrete) specifying which component should be returned
##' @return matrix of standard errors
##' @importMethodsFrom arm se.coef
##' @export
setMethod('se.coef', signature=c(object='ZlmFit'), function(object, which, ...){
    which <- match.arg(which, c('C', 'D'))
    vc <- if(which=='C') object@vcovC else object@vcovD
    se <- sqrt(aaply(vc, 3, diag))
    rownames(se) <- fData(object@sca)$primerid
    se
})

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
    stopifnot(all.equal(dvcont[1,], diag(vcont[1,,]), check.attributes=FALSE))
    stopifnot(all.equal(dvdisc[10,], diag(vdisc[10,,]), check.attributes=FALSE))
    ## variance of product is product of variance + higher order terms 
    v.prod <- dvcont*dvdisc + mu.cont^2*dvdisc + mu.disc^2*dvcont
    
    ## covariance between contrast1 and contrast0
    covc1c0Idx <- diagIdx[-(seq_along(genes)),] #drop j=k=1 indices, ie, 1,1 entry of covariances
    covc1c0Idx[,'k'] <- 1
    c1c0cont <- matrix(vcont[covc1c0Idx], nrow=length(genes))
    c1c0disc <- matrix(vdisc[covc1c0Idx], nrow=length(genes))
    ## debugging assertion
    stopifnot(all.equal(c1c0cont[1,], vcont[1,2:nrow(Contr),1], check.attributes=FALSE))
    stopifnot(all.equal(c1c0disc[10,], vdisc[10,2:nrow(Contr),1], check.attributes=FALSE))
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
    stopifnot(all.equal(res[3,], rowvec[3]*rowMajorMatrix[3,]))
    res
}
