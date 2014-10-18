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

##' Bootstrap a zlmfit
##'
##' Sample cells with replacement to find bootstrapped distribution of coefficients
##' @param zlmfit class \code{ZlmFit}
##' @param R number of bootstrap replicates
##' @return array of bootstrapped coefficients
##' @importFrom plyr raply
##' @export
bootVcov1 <- function(zlmfit, R=999){
    sca <- zlmfit@sca
    N <- nrow(sca)
    manyvc <- raply(R, {
        s <- sample(N, replace=TRUE)
        newsca <- sca[s,]
        z <- zlm.SingleCellAssay(sca=newsca, LMlike=zlmfit@LMlike)
        abind(C=coef(z, 'C'), D=coef(z, 'D'), rev.along=0)
    })

   manyvc
    
}

## use dfbetas? Incomplete...
bootVcov2 <- function(zlmfit, R=999){
    sca <- zlmfit@sca
    N <- nrow(sca)
    z <- zlm.SingleCellAssay(sca=newsca, LMlike=zlmfit@LMlike, hook=function(lml){
        cont <- lml@fitC
        disc <- lml@fitD
        class(disc) <-class(cont) <- 'glm'
        dfc <- lm.influence(cont)$coefficients
        dfd <- lm.influence(disc)$coefficients
        list(cont=dfc, disc=dfd)
        })
    ## do stuff to this
}


