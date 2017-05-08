## This is a horrible hack, and should be rewritten to use the lmer internals.
## But in the meantime: we make the fixed effects model matrix, then call the formula method for lmer/glmer
## This is to allow us to do arbitrary LRT tests and add/drop columns of the design

## Details:
## Invariants:
## 1. model.matrix contains fixed effects--so whenever we set the model.matrix, we'll delete the random effects
## 2. formula contains full model (fixed and random), so we can update it normally
## 3. Random portion of the model will be parsed off

## Construction:
## 1 & 2
## Fitting:
## establish pseudodesign and mutilate the formula

getREvars <- function(Formula){
    termNames <- labels(terms(Formula))
    hasRE <- str_detect(termNames, fixed('|'))
    ## collapse all variables into something that can be used for model.frame
    REvar <- str_replace_all(paste(termNames[hasRE], collapse='+', sep='+'), '[|]+', '+')
    ## save portion of formula that contained random effects
    REform <- paste(sprintf('(%s)', termNames[hasRE]), collapse='+')   
    FEform<- paste(sprintf('%s', termNames[!hasRE]), collapse='+')
    if(str_trim(FEform)=='') FEform <- '1'

    ## REvar: Random effects variables concatenated with +
    ## REform: the actual formula specifying the random effects
    ## FEform: the actual formula specifying the fixed effects
    ## All are character vectors of length 1.
    list(vars=REvar, REform=REform, FEform=FEform)
}

toAdditiveString <- function(string){
    if(length(string)>1)
        string <- paste(string, collapse='+')
    paste0('~', string)
}

toAdditiveFormula <- function(string){
    string <- as.formula(toAdditiveString(string))
}

##' @export
##' @describeIn LMERlike update the formula or design matrix
##' @param formula. \code{formula}
##' @param design  something coercible to a \code{data.frame}
setMethod('update', signature=c(object='LMERlike'), function(object, formula., design, ...){
    if(!missing(formula.)){
        object@formula <- update.formula(object@formula, formula.)
    }
    reComponents <- getREvars(object@formula)
    if(!missing(design)){
        object@design <- as(design, 'data.frame')
    }
    model.matrix(object) <- model.matrix(as.formula(paste0('~', reComponents$FEform)), object@design, ...)
    object@fitC <- object@fitD <- numeric(0)
    object@fitted <- c(C=FALSE, D=FALSE)
    object
})


setMethod('initialize', 'LMERlike', function(.Object, ...){
    .Object <- callNextMethod()
    reComponents <- getREvars(.Object@formula)
    model.matrix(.Object) <- model.matrix(as.formula(paste0('~', reComponents$FEform)), .Object@design)
    .Object
})


setReplaceMethod('model.matrix', signature=c(object='LMERlike'), function(object, value){
    reComponents <- getREvars(object@formula)
    object <- callNextMethod()
    object@pseudoMM <- as.data.frame(cbind(model.matrix(object),
                                           model.frame(toAdditiveFormula(reComponents$vars), object@design)))
    object
})

## lmerMM <- function (formula, data = NULL, REML = TRUE, control = lmerControl(), 
##     start = NULL, verbose = 0L, subset, weights, na.action, offset, 
##     contrasts = NULL, devFunOnly = FALSE, modelMatrix, ...) 
## {
##     mc <- mcout <- match.call()
##     missCtrl <- missing(control)
##     if (!missCtrl && !inherits(control, "lmerControl")) {
##         if (!is.list(control)) 
##             stop("'control' is not a list; use lmerControl()")
##         warning("passing control as list is deprecated: please use lmerControl() instead", 
##             immediate. = TRUE)
##         control <- do.call(lmerControl, control)
##     }
##     if (!is.null(list(...)[["family"]])) {
##         warning("calling lmer with 'family' is deprecated; please use glmer() instead")
##         mc[[1]] <- quote(lme4::glmer)
##         if (missCtrl) 
##             mc$control <- glmerControl()
##         return(eval(mc, parent.frame(1L)))
##     }
##     mc$control <- control
##     mc[[1]] <- quote(lme4::lFormula)
##     lmod <- eval(mc, parent.frame(1L))
##     lmod$X <- modelMatrix
##     mcout$formula <- lmod$formula
##     lmod$formula <- NULL
##     devfun <- do.call(mkLmerDevfun, c(lmod, list(start = start, 
##         verbose = verbose, control = control)))
##     if (devFunOnly) 
##         return(devfun)
##     opt <- optimizeLmer(devfun, optimizer = control$optimizer, 
##         restart_edge = control$restart_edge, boundary.tol = control$boundary.tol, 
##         control = control$optCtrl, verbose = verbose, start = start, 
##         calc.derivs = control$calc.derivs, use.last.params = control$use.last.params)
##     cc <- checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, 
##         lbound = environment(devfun)$lower)
##     mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr, 
##         mcout, lme4conv = cc)
## }

## glmerMM <- function (formula, data = NULL, family = gaussian, control = glmerControl(), 
##     start = NULL, verbose = 0L, nAGQ = 1L, subset, weights, na.action, 
##     offset, contrasts = NULL, mustart, etastart, devFunOnly = FALSE, modelMatrix,
##     ...) 
## {
##     if (!inherits(control, "glmerControl")) {
##         if (!is.list(control)) 
##             stop("'control' is not a list; use glmerControl()")
##         msg <- "Use control=glmerControl(..) instead of passing a list"
##         if (length(cl <- class(control))) 
##             msg <- paste(msg, "of class", dQuote(cl[1]))
##         warning(msg, immediate. = TRUE)
##         control <- do.call(glmerControl, control)
##     }
##     mc <- mcout <- match.call()
##     if (is.character(family)) 
##         family <- get(family, mode = "function", envir = parent.frame(2))
##     if (is.function(family)) 
##         family <- family()
##     if (isTRUE(all.equal(family, gaussian()))) {
##         warning("calling glmer() with family=gaussian (identity link) as a shortcut to lmer() is deprecated;", 
##             " please call lmer() directly")
##         mc[[1]] <- quote(lme4::lmer)
##         mc["family"] <- NULL
##         return(eval(mc, parent.frame()))
##     }
##     mc[[1]] <- quote(lme4::glFormula)
##     glmod <- eval(mc, parent.frame(1L))
##     glmod$X <- modelMatrix
##     mcout$formula <- glmod$formula
##     glmod$formula <- NULL
##     devfun <- do.call(mkGlmerDevfun, c(glmod, list(verbose = verbose, 
##         control = control, nAGQ = 0)))
##     if (nAGQ == 0 && devFunOnly) 
##         return(devfun)
##     if (is.list(start) && !is.null(start$fixef)) 
##         if (nAGQ == 0) 
##             stop("should not specify both start$fixef and nAGQ==0")
##     opt <- optimizeGlmer(devfun, optimizer = control$optimizer[[1]], 
##         restart_edge = if (nAGQ == 0) 
##             control$restart_edge
##         else FALSE, boundary.tol = if (nAGQ == 0) 
##             control$boundary.tol
##         else 0, control = control$optCtrl, start = start, nAGQ = 0, 
##         verbose = verbose, calc.derivs = FALSE)
##     if (nAGQ > 0L) {
##         start <- updateStart(start, theta = opt$par)
##         devfun <- updateGlmerDevfun(devfun, glmod$reTrms, nAGQ = nAGQ)
##         if (devFunOnly) 
##             return(devfun)
##         opt <- optimizeGlmer(devfun, optimizer = control$optimizer[[2]], 
##             restart_edge = control$restart_edge, boundary.tol = control$boundary.tol, 
##             control = control$optCtrl, start = start, nAGQ = nAGQ, 
##             verbose = verbose, stage = 2, calc.derivs = control$calc.derivs, 
##             use.last.params = control$use.last.params)
##     }
##     cc <- if (!control$calc.derivs) 
##         NULL
##     else {
##         if (verbose > 10) 
##             cat("checking convergence\n")
##         checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, 
##             lbound = environment(devfun)$lower)
##     }
##     mcout <- call('LMERlike')
##     mkMerMod(environment(devfun), opt, glmod$reTrms, fr = glmod$fr, 
##         mcout, lme4conv = cc)
## }


##' @include AllClasses.R
##' @include AllGenerics.R
##' @param silent mute some warnings emitted from the underlying modeling functions
##' @rdname fit
setMethod('fit', signature=c(object='LMERlike', response='missing'), function(object, response, silent=TRUE, ...){
    prefit <- .fit(object)
    if(!prefit){
        if(!silent) warning('No positive observations')
        return(object)
    }

    fitArgsC <- object@fitArgsC
    fitArgsD <- object@fitArgsD

    ## Mutilate the formula and replace it with the colnames of the fixed effects
    ## 
    reComp <- getREvars(object@formula)
    protoForm <- sprintf('~ 0 + %s + %s',
                         paste(escapeSymbols(colnames(model.matrix(object))), collapse='+'),
                         reComp$REform)
    
    formC <- as.formula(paste0('response ', protoForm))
    formD <- as.formula(paste0('response>0', protoForm))
    dat <- cbind(response=object@response, object@pseudoMM)
    
    if(inherits(object, 'bLMERlike')){
        cfun <- blme::blmer
        dfun <- blme::bglmer
    } else{
        cfun <- lme4::lmer
        dfun <- lme4::glmer
    }
    
    if(any(pos)){
        datpos <- dat[pos,]
        object@fitC <- do.call(cfun, c(list(formula=formC, data=quote(datpos), REML=FALSE), fitArgsC))
        ok <- length(object@fitC@optinfo$conv$lme4)==0
        object@fitted['C'] <- TRUE
        if(!ok){
            object@optimMsg['C'] <- object@fitC@optinfo$conv$lme4$messages[1]
            object@fitted['C'] <- !object@strictConvergence
        }
    }
    if(!all(pos)){
        object@fitD <- do.call(dfun, c(list(formula=formD, data=quote(dat), family=binomial()), fitArgsD))
        object@fitted['D'] <- length(object@fitD@optinfo$conv$lme)==0
        ok <- length(object@fitD@optinfo$conv$lme4)==0
        object@fitted['D'] <- TRUE
        if(!ok){
            object@optimMsg['D'] <- object@fitD@optinfo$conv$lme4$messages[[1]]
            object@fitted['D'] <- !object@strictConvergence
        }
    } 
    if(!silent & !all(object@fitted)) warning('At least one component failed to converge')
    object    

})

#' @describeIn LMERlike return the variance/covariance of component \code{which}
#' @param object \code{LMERlike}
#' @param which \code{character}, one of 'C', 'D'.
#' @param ... In the case of \code{vcov}, ignored.  In the case of \code{update}, passed to \code{model.matrix}.
#' @return see the section "Methods (by generic)"
setMethod('vcov', signature=c(object='LMERlike'), function(object, which, ...){
    stopifnot(which %in% c('C', 'D'))
    vc <- object@defaultVcov

    if(which=='C' & object@fitted['C']){
        V <- vcov(object@fitC)
    } else if(which=='D' & object@fitted['D']){
        V <- vcov(object@fitD)
    } else{
        V <- matrix(nrow=0, ncol=0)
    }
    nm <- str_replace_all(colnames(V), fixed('`'), '')
    dimnames(V) <- list(nm, nm)
    ok <- colnames(V)
    vc[ok,ok] <- as.numeric(V)
    vc
})

if(getRversion() >= "2.15.1") globalVariables(c('fixef', 'lmer', 'glmer'))
#' @describeIn LMERlike return the coefficients.  The horrendous hack is attempted to be undone.
#' @param singular \code{logical}. Should NA coefficients be returned?
setMethod('coef', signature=c(object='LMERlike'), function(object, which, singular=TRUE, ...){
    stopifnot(which %in% c('C', 'D'))
    co <- setNames(rep(NA, ncol(model.matrix(object))), colnames(model.matrix(object)))
    if(which=='C' & object@fitted['C']){
        co <- fixef(object@fitC)}
    else if(object@fitted['D']){
        co <- fixef(object@fitD)
    }
    if(!singular) co <- co[!is.na(co)]
    conm <- names(co)
    ## because of backtick shenanigans
    names(co) <- str_replace_all(conm, fixed('`'), '')
    co
})

##' @describeIn LMERlike return the log-likelihood
setMethod('logLik', signature=c(object='LMERlike'), function(object){
    L <- c(C=0, D=0)
    if(object@fitted['C']) L['C'] <- logLik(object@fitC)
    if(object@fitted['D']) L['D'] <- logLik(object@fitD)
    L
})

setMethod('dof', signature=c(object='LMERlike'), function(object){
    setNames(ifelse(object@fitted, c(attr(logLik(object@fitC), 'df'), attr(logLik(object@fitD), 'df')), c(0,0)), c('C', 'D'))

})

setMethod('summarize', signature=c(object='LMERlike'), function(object, ...){

    li <- list(coefC=coef(object, which='C'), vcovC=vcov(object, 'C'),
               deviance=rowm(deviance(object@fitC), deviance(object@fitD)),
               df.null=rowm(nobs(object@fitC),nobs(object@fitD)),
               dispersion=rowm(sigma(object@fitC), NA),
               coefD=coef(object, which='D'), vcovD=vcov(object, 'D'),
               loglik=torowm(logLik(object)),
               converged=torowm(object@fitted))
    
    li[['df.resid']] <- li[['df.null']]-c(sum(!is.na(li[['coefC']])), sum(!is.na(li[['coefD']])))
    li[['dispersionNoshrink']] <- li[['dispersion']]
    li
})
