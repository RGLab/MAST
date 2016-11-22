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
    testable <- o0@converged & o1@converged
    testable[,'C'] <-  testable[,'C'] & (o0@df.resid > 1 & o1@df.resid > 1)[,'C']
    lambda <- ifelse(testable, lambda, 0)
    df <- o0@df.resid-o1@df.resid
    df <- ifelse(testable, df, 0)
    cst <- makeChiSqTable(as.data.frame(lambda), as.data.frame(df), hString)
    dimnames(cst) <- list(primerid=mcols(zlmfit@sca)$primerid, test.type=dimnames(cst)[[2]], metric=dimnames(cst)[[3]])
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
    oldMM <- model.matrix(LMlike)
    newF <- update.formula(LMlike@formula, formula(sprintf(' ~. - %s', hypothesis)))
    LMlike <- update(LMlike, newF)
    newMM <- model.matrix(LMlike)

    ##test if design matrix has different column space
    ## Same if Q * Q^T X = X
    ## (old MM lies in the projection onto the new MM)
    ## There must be a better way to do this??
    newq <- qr.Q(qr(newMM))
    diff <- any(abs(newq %*% crossprod(newq, oldMM) - oldMM)>1e-6)
    if(!diff) stop('Removing term ', sQuote(hypothesis), " doesn't actually alter the model, maybe due to marginality? Try specifying individual coefficents as a `CoefficientHypothesis`.")
    .lrtZlmFit(o1, LMlike@modelMatrix, hypothesis)
})

##' @describeIn ZlmFit Returns an array with likelihood-ratio tests on contrasts defined using \code{CoefficientHypothesis()}.
##' @return see "Methods (by generic)"
setMethod('lrTest', signature=c(object='ZlmFit', hypothesis='CoefficientHypothesis'), function(object, hypothesis){
    h <- generateHypothesis(hypothesis, colnames(object@coefD))
    testIdx <- h@index
    newMM <- model.matrix(object@LMlike)[,-testIdx, drop=FALSE]
    .lrtZlmFit(object, newMM, hypothesis@.Data)
})

##' @describeIn ZlmFit Returns an array with likelihood-ratio tests specified by \code{Hypothesis}, which is a \link{Hypothesis}.
setMethod('lrTest', signature=c(object='ZlmFit', hypothesis='Hypothesis'), function(object, hypothesis){
    ## original fit
    h <- generateHypothesis(hypothesis, colnames(object@coefD))
    ## call using coefficient matrix
    lrTest(object, h@contrastMatrix)
})

## contrast matrices
##' @describeIn ZlmFit Returns an array with likelihood-ratio tests specified by \code{Hypothesis}, which is a contrast \code{matrix}.
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
    h <- generateHypothesis(hypothesis,  colnames(object@coefD))
    waldTest(object, h@contrastMatrix)
})

##' @param hypothesis call to \link{Hypothesis} or \link{CoefficientHypothesis} or a matrix giving such contrasts.
##' @describeIn ZlmFit Returns an array with Wald Tests on contrasts defined in \code{Hypothesis()}
setMethod('waldTest',  signature=c(object='ZlmFit', hypothesis='Hypothesis'), function(object, hypothesis){
    h <- generateHypothesis(hypothesis, colnames(object@coefD))
    waldTest(object, h@contrastMatrix)
})

##' @describeIn show print info on ZlmFit
setMethod('show', signature=c(object='ZlmFit'), function(object){
    cat('Fitted zlm on', nrow(object@sca), 'genes and', ncol(object@sca), 'cells.\n Using', class(object@LMlike), as.character(object@LMlike@formula), '\n')
})

##' @describeIn ZlmFit Returns the matrix of coefficients for component \code{which}.
##' @param ... ignored
##' @param object \code{ZlmFit}
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
##' @describeIn ZlmFit Returns a matrix of standard error estimates for coefficients on component \code{which}.
##' @export
setMethod('se.coef', signature=c(object='ZlmFit'), function(object, which, ...){
    which <- match.arg(which, c('C', 'D'))
    vc <- if(which=='C') object@vcovC else object@vcovD
    se <- sqrt(aaply(vc, 3, diag))
    rownames(se) <- mcols(object@sca)$primerid
    se
})

normalci <- function(center, se, level){
    zstar <- qnorm(1-(1-level)/2)
    list(coef=center, ci.lo=center-se*zstar, ci.hi=center+se*zstar)
}

##' Summarize model features from a \code{ZlmFit} object
##'
##' Returns a \code{data.table} with a special print method that shows the top 2 most significant genes by contrast.
##' This \code{data.table} contains columns: 
##' \describe{
#'   \item{primerid}{the gene}
#'   \item{component}{C=continuous, D=discrete, logFC=log fold change, S=combined using Stouffer's method, H=combined using hurdle method}
#' \item{contrast}{the coefficient/contrast of interest}
#' \item{ci.hi}{upper bound of confidence interval}
#' \item{ci.lo}{lower bound of confidence interval}
#' \item{coef}{point estimate}
#' \item{z}{z score (coefficient divided by standard error of coefficient)}
#' \item{Pr(>Chisq)}{likelihood ratio test p-value (only if \code{doLRT=TRUE})}
#' }
#' Some of these columns will contain NAs if they are not applicable for a particular component or contrast.
##' @param object A \code{ZlmFit} object
##' @param logFC If TRUE, calculate log-fold changes, or output from a call to \code{getLogFC}.
##' @param doLRT if TRUE, calculate lrTests on each coefficient, or a character vector of such coefficients to consider.
##' @param level what level of confidence coefficient to return.  Defaults to 95 percent. 
##' @param ... ignored
##' @return \code{data.table}
##' @seealso print.summaryZlmFit
##' @examples
##' data(vbetaFA)
##' z <- zlm(~Stim.Condition, vbetaFA[1:5,])
##' zs <- summary(z)
##' names(zs)
##' print(zs)
##' ##Select `datatable` copmonent to get normal print method
##' zs$datatable
##' @export
setMethod('summary', signature=c(object='ZlmFit'), function(object, logFC=TRUE,  doLRT=FALSE, level=.95, ...){
    message('Combining coefficients and standard errors')

    
    coefAndCI <- aaply( c(C='C', D='D'), 1, function(component){
        ## coefficients for each gene
        coefs <- coef(object, which=component)
        ## standard errors for each gene
        se <- se.coef(object, which=component)
        names(dimnames(se)) <- names(dimnames(coefs)) <- c('primerid', 'Coefficient')
        ci <- normalci(coefs, se, level)
        z <- coefs/se
        abind(coef=coefs, z=z, ci.lo=ci[['ci.lo']], ci.hi=ci[['ci.hi']], rev.along=0, hier.names=TRUE)
    })
    names(dimnames(coefAndCI)) <- c('component', 'primerid', 'contrast', 'metric')
    dt <- dcast.data.table(data.table(melt(coefAndCI, as.is=TRUE)), primerid + component + contrast ~ metric)
    setkey(dt, primerid, contrast)
    stouffer <- dt[,list(z=sum(z)/sqrt(sum(!is.na(z))), component='S'), keyby=list(primerid, contrast)]
    dt <- rbind(dt, stouffer, fill=TRUE)
    
    if(is.logical(logFC) && logFC){
        message("Calculating log-fold changes")
        logFC <- getLogFC(zlmfit=object)
    }
    
    if(!is.logical(logFC)){
        lfc <- logFC
        ci <- normalci(lfc$logFC, sqrt(lfc$varLogFC), level=level)
        setnames(lfc, 'logFC', 'coef')
        lfc <- lfc[,c('component', 'ci.lo', 'ci.hi', 'varLogFC'):=list('logFC', ci$ci.lo, ci$ci.hi, NULL)]
        dt <- rbind(dt, lfc, fill=TRUE)
    }
    setkey(dt, contrast, z)

    if(is.logical(doLRT) && doLRT){
        doLRT <- setdiff(colnames(coef(object, 'D')), '(Intercept)')
    }
    if(!is.logical(doLRT)){
        message('Calculating likelihood ratio tests')
        llrt <- lapply(doLRT, function(x) lrTest(object, CoefficientHypothesis(x))[,,'Pr(>Chisq)'])
        names(llrt) <-  doLRT
        llrt <- data.table(melt(llrt))
        setnames(llrt, c('test.type', 'L1', 'value'), c('component', 'contrast', 'Pr(>Chisq)'))
        llrt[,':='(component=c(cont='C', disc='D', hurdle='H')[component],
                   primerid=as.character(primerid))]
        setkey(llrt, primerid, component, contrast)
        setkey(dt, primerid, component, contrast)
        dt <- merge(llrt, dt, all.x=TRUE, all.y=TRUE)
        dt[is.na(component), component :='H']
    }
    out <- list(datatable=dt)
    class(out) <- c('summaryZlmFit', 'list')
    out
})

if(getRversion() >= "2.15.1") globalVariables(c(
                                  'contrast',
                                  'metric', 
                                  '.',
                                  'value'))

##' Print summary of a ZlmFit
##'
##' Shows the top `n` genes by z score on `by`
##' @param x output from summary(ZlmFit)
##' @param n number of genes to show
##' @param by one of 'C' , 'D' or 'logFC' for continuous, discrete and log fold change z-scores for each contrast
##' @param ... ignored
##' @return prints a pretty table and invisibly returns a \code{data.table} representing the table.
##' @seealso summary,ZlmFit-method
##' @export
print.summaryZlmFit <- function(x, n=2, by='logFC', ...){
    x <- x$datatable
    ns <- seq_len(n)
    dt <- x[contrast!='(Intercept)',]
    by <- match.arg(by, c('logFC', 'D', 'C'))
    if(length(intersect(dt$component, 'logFC'))==0 & by=='logFC'){
        warning('Log fold change not calculated, selecting discrete', call.=FALSE)
        by <- 'D'
    }
    dt[,metric:=ifelse(is.na(z), 0, -abs(z))]
    describe <- switch(by, logFC='log fold change Z-score', D='Wald Z-scores on discrete', C='Wald Z-scores tests on continuous')
    
    ## if('Pr(>Chisq)' %in% names(dt)){
    ##     keyv <- 'contrast'
    ##     dt[,metric:=`Pr(>Chisq)`]
    ##     jj <- quote(.(primerid=primerid[ns], value=`Pr(>Chisq)`[ns]))
    ##     sel <- quote(component=='H')
    ## } else
    dt <- dt[component==by,]
    setkey(dt, contrast, metric)
    pid <- dt[,list(primerid=primerid[ns]),by=contrast]
    setkey(dt, primerid)
    setkey(pid, primerid, contrast)
    dts <- dt[unique(pid[,list(primerid)]), list(contrast=contrast, primerid=primerid, z=z)]
    setkey(dts, primerid, contrast)
    dts[pid,value:=sprintf('%6.1f*',  z)]
    dts[!pid,value:=sprintf('%6.1f ',  z)]
    cdts <- as.data.frame(dcast.data.table(dts, primerid ~ contrast))
    cat('Fitted zlm with top', n, 'genes per contrast:\n')
    cat('(', describe, ')\n')
    print(cdts, right=FALSE, row.names=FALSE)
    invisible(cdts)
}
