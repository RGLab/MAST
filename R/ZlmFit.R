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

##'  @describeIn ZlmFit Returns an array with likelihood-ratio tests on contrasts defined using \code{CoefficientHypothesis()}.
setMethod('lrTest', signature=c(object='ZlmFit', hypothesis='CoefficientHypothesis'), function(object, hypothesis){
    h <- generateHypothesis(hypothesis, colnames(object@coefD))
    testIdx <- h@transformed
    newMM <- model.matrix(object@LMlike)[,-testIdx, drop=FALSE]
    .lrtZlmFit(object, newMM, hypothesis@.Data)
})

##' @describeIn ZlmFit Returns an array with likelihood-ratio tests specified by \code{Hypothesis}, which can be a call to \link{Hypothesis} or \link{CoefficientHypothesis} or a contrast matrix.
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
##' 
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

##' @describeIn ZlmFit Returns an array with Wald Tests on contrasts defined using \code{CoefficientHypothesis()}.
setMethod('waldTest',  signature=c(object='ZlmFit', hypothesis='CoefficientHypothesis'), function(object, hypothesis){
    cm <- .makeContrastMatrixFromCoefficientHypothesis(hypothesis,colnames(object@coefC))
    waldTest(object, cm)
})

##' @param hypothesis call to \link{Hypothesis} or \link{CoefficientHypothesis} or a matrix giving such contrasts.
##' @describeIn ZlmFit Returns an array with Wald Tests on contrasts defined in \code{Hypothesis()}
setMethod('waldTest',  signature=c(object='ZlmFit', hypothesis='Hypothesis'), function(object, hypothesis){
    h <- generateHypothesis(hypothesis, colnames(object@coefD))
    waldTest(object, h@transformed)
})

setMethod('show', signature=c(object='ZlmFit'), function(object){
    cat('Fitted zlm on ', ncol(object@sca), ' genes and ', nrow(object@sca), ' cells.\n Using ', class(object@LMlike), ' to fit.\n')
})

##' @describeIn ZlmFit Returns the matrix of coefficients for component \code{which}.
setMethod('coef', signature=c(object='ZlmFit'), function(object, which, ...){
    which <- match.arg(which, c('C', 'D'))
    if(which=='C') object@coefC else object@coefD
})

##' @describeIn ZlmFit Returns an array of variance/covariance matrices for component \code{which}.
##' @export
setMethod('vcov', signature=c(object='ZlmFit'), function(object, which, ...){
    which <- match.arg(which, c('C', 'D'))
    if(which=='C') object@vcovC else object@vcovD
})


##' @param which  character vector, one of "C" (continuous) or "D" (discrete) specifying which component should be returned
##' @importMethodsFrom arm se.coef
##' @describeIn ZlmFit Returns a matrix of standard error estimates for coefficients on component \code{which}.
##' @export
setMethod('se.coef', signature=c(object='ZlmFit'), function(object, which, ...){
    which <- match.arg(which, c('C', 'D'))
    vc <- if(which=='C') object@vcovC else object@vcovD
    se <- sqrt(aaply(vc, 3, diag))
    rownames(se) <- fData(object@sca)$primerid
    se
})

##' @param object ZlmFit
##' @param logFC If TRUE, calcualte log-fold changes, or output from a call to \code{getLogFC}.
##' @param ... \code{waldTests} or \code{lrTests} on \code{ZlmFit} to be combined in the output
##' @export
##' @describeIn ZlmFit  Returns a \code{data.table} summary of fit (invisibly).
setMethod('summary', signature=c(object='ZlmFit'), function(object, logFC=FALSE,  ...){
    message('Combining coefficients and standard errors')
    coefAndCI <- aaply( c(C='C', D='D'), 1, function(component){
        ## coefficients for each gene
        coefs <- coef(object, which=component)
        ## standard errors for each gene
        se2 <- se.coef(object, which=component)*2
        names(dimnames(se2)) <- names(dimnames(coefs)) <- c('primerid', 'Coefficient')
        ci.lo <- coefs-se2
        ci.hi <- coefs+se2
        z <- coefs/(se2/2)
        abind(coef=coefs, z=z, ci.lo=ci.lo, ci.hi=ci.hi, rev.along=0, hier.names=TRUE)
    })
    names(dimnames(coefAndCI)) <- c('component', 'primerid', 'contrast', 'metric')
    dt <- dcast.data.table(data.table(melt(coefAndCI)), primerid + component + contrast ~ metric)
    setkey(dt, primerid, contrast)
 
    if(is.logical(logFC) && logFC){
        message("Calculating log-fold changes")
        logFC <- getLogFC(zlmfit=object)
    }
    
    if(!is.logical(logFC)){
        lfc <- logFC
        lfc[,se2:=2*sqrt(varLogFC)]
        setnames(lfc, 'logFC', 'coef')
        lfc <- lfc[,c('component', 'ci.lo', 'ci.hi', 'varLogFC', 'se2'):=list('logFC', coef-se2, coef+se2, NULL, NULL)]
        dt <- rbind(dt, lfc, fill=TRUE)
    }
    
    vargs <- list(...)
    if(length(vargs)>0){
        browser()
        MTests <- do.call(melt, vargs)
    }
    return(dt)
})
