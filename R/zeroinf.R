methodDict <- data.table(keyword=c('glm', 'glmer', 'lmer', 'bayesglm','ridge'),
                         lmMethod=c('GLMlike', 'LMERlike','LMERlike', 'BayesGLMlike','RidgeBGLMlike'),
                         implementsEbayes=c(TRUE, FALSE, FALSE, TRUE, TRUE))


if(getRversion() >= "2.15.1") globalVariables(c(
                  'keyword',
                 'lmMethod', 
                  'implementsEbayes')) #zlm, zlm.SingleCellAssay


residualsHook <- function(fit){
    residuals(fit, which='Marginal')
}

revealHook <- function(zlm){
    return(attr(zlm, 'hookOut'))
}

##' @importFrom plyr laply
collectResiduals <- function(zlm, sca, newLayerName='Residuals'){
    if(any(newLayerBool <- assayNames(sca) %in% newLayerName)){
        warning('Overwriting layer', newLayerName)
        i <- which(newLayerBool)
    } else{
        i <- length(assays(sca))+1
    }
    mat <- laply(revealHook(zlm), function(x) x)
    assay(sca, i) <- mat
    assayNames(sca, i) <- newLayerName
    sca
}

## for each column in exprMat, fit obj and return an environment containing the init summaries
## exprMat should have column names containing the names of the genes
## hook is an (optional) function called on each fitted object
                      

##' Convenience function for running a zero-inflated regression
##'
##' Fits a hurdle model on zero-inflated continuous data in which the zero process
##' is modeled as a logistic regression
##' and (conditional on the the response being >0), the continuous process is Gaussian, ie, a linear regression.
##' @param formula model formula
##' @param data a data.frame, list, environment or SingleCellAssay in which formula is evaluated
##' @param method one of 'glm', 'glmer' or 'bayesglm'.  See MAST:::methodDict for other possibilities.
##' @param silent if TRUE suppress common errors from fitting continuous part
##' @param ... passed to \code{fit}, and eventually to the linear model fitting function
##' @return list with "disc"rete part and "cont"inuous part 
##' @export
##' @examples
##' data<- data.frame(x=rnorm(500), z=rbinom(500, 1, .3))
##' logit.y <- with(data, x*2 + z*2); mu.y <- with(data, 10+10*x+10*z + rnorm(500))
##' y <- (runif(500)<exp(logit.y)/(1+exp(logit.y)))*1
##' y[y>0] <- mu.y[y>0]
##' data$y <- y
##' fit <- zlm(y ~ x+z, data)
##' summary.glm(fit$disc)
##' summary.glm(fit$cont)
##'
##' @seealso GLMlike, LMERlike, BayesGLMlike
##' @import stringr
zlm <- function(formula, data, method='bayesglm',silent=TRUE, ...){
    ## perhaps we should be generic, but since we are dispatching on second argument, which might be an S3 class, let's just do this instead.
    if(inherits(data, 'SingleCellAssay')){
        return(zlm.SingleCellAssay(formula, sca=data, method=method, silent=silent, ...))
    }
    
    if(!inherits(data, 'data.frame')) stop("'data' must be data.frame, not matrix or array")
    if(!is(formula, 'formula')) stop("'formula' must be class 'formula'")

    ## get response
    resp <- eval(formula[[2]], data)
    RHS <- removeResponse(formula, warn=FALSE)
    
    obj <- new(methodDict[keyword==method, lmMethod], formula=RHS, design=data, response=resp)
    obj <- fit(obj)
    list(cont=obj@fitC, disc=obj@fitD)
}

summary.zlm <- function(out){
  summary(out$cont)
  summary(out$disc)
}
 
##' Zero-inflated regression for SingleCellAssay 
##'
##' For each gene in sca, fits the hurdle model in \code{formula} (linear for et>0), logistic for et==0 vs et>0.
##' Return an object of class \code{ZlmFit} containing slots giving the coefficients, variance-covariance matrices, etc.
##' After each gene, optionally run the function on the fit named by 'hook'
##'
##' @section Empirical Bayes variance regularization:
##' The empirical bayes regularization of the gene variance assumes that the precision (1/variance) is drawn from a
##' gamma distribution with unknown parameters.
##' These parameters are estimated by considering the distribution of sample variances over all genes.
##' The procedure used for this is determined from
##' \code{ebayesControl}, a named list with components 'method' (one of 'MOM' or 'MLE') and 'model' (one of 'H0' or 'H1')
##' method MOM uses a method-of-moments estimator, while MLE using the marginal likelihood.
##' H0 model estimates the precisions using the intercept alone in each gene, while H1 fits the full model specified by \code{formula}
##'
##' @param formula a formula with the measurement variable on the LHS and predictors present in colData on the RHS
##' @param sca SingleCellAssay object
##' @param method character vector, either 'glm', 'glmer' or 'bayesglm'
##' @param silent Silence common problems with fitting some genes
##' @param ebayes if TRUE, regularize variance using empirical bayes method
##' @param ebayesControl list with parameters for empirical bayes procedure.  See \link{ebayes}.
##' @param force Should we continue testing genes even after many errors have occurred?
##' @param hook a function called on the \code{fit} after each gene.
##' @param parallel If TRUE and \code{option(mc.cores)>1} then multiple cores will be used in fitting.
##' @param LMlike if provided, then the model defined in this object will be used, rather than following the formulas.  This is intended for internal use.
##' @param onlyCoef If TRUE then only an array of model coefficients will be returned (probably only useful for bootstrapping).
##' @param ... arguments passed to the S4 model object upon construction.  For example, \code{fitArgsC} and \code{fitArgsD}, or \code{coefPrior}.
##' @return a object of class \code{ZlmFit} with methods to extract coefficients, etc.
##' @export
##' @seealso ebayes, glmlike-class, ZlmFit-class, BayesGLMlike-class
##' @examples
##' \dontrun{
##' data(vbetaFA)
##' zlmVbeta <- zlm.SingleCellAssay(~ Stim.Condition, subset(vbetaFA, ncells==1))
##' slotNames(zlmVbeta)
##' coef(zlmVbeta, 'D')['CD27',]
##' 
##' vcov(zlmVbeta, 'D')[,,'CD27']
##' waldTest(zlmVbeta, CoefficientHypothesis('Stim.ConditionUnstim'))
##' }
zlm.SingleCellAssay <- function(formula, sca, method='bayesglm', silent=TRUE, ebayes=TRUE, ebayesControl=NULL, force=FALSE, hook=NULL, parallel=TRUE, LMlike, onlyCoef=FALSE, ...){
    ## Default call
    if(missing(LMlike)){
        ## Which class are we using for the fits...look it up by keyword
        method <- match.arg(method, methodDict[,keyword])
        method <- methodDict[keyword==method,lmMethod]
        
        if(!is(sca, 'SingleCellAssay')) stop("'sca' must be (or inherit) 'SingleCellAssay'")
        if(!is(formula, 'formula')) stop("'formula' must be class 'formula'")
        Formula <- removeResponse(formula)

        ## Empirical bayes method
        priorVar <- 1
        priorDOF <- 0
        if(ebayes){
            if(!methodDict[lmMethod==method,implementsEbayes]) stop('Method', method, ' does not implement empirical bayes variance shrinkage.')
            ebparm <- ebayes(sca, ebayesControl, Formula)
            priorVar <- ebparm[['v']]
            priorDOF <- ebparm[['df']]
            stopifnot(all(!is.na(ebparm)))
        }
        ## initial value of priorVar, priorDOF default to no shrinkage
        obj <- new(method, design=colData(sca), formula=Formula, priorVar=priorVar, priorDOF=priorDOF, ...)
        ## End Default Call
    } else{
        ## Refitting
        if(!missing(formula)) warning("Ignoring formula and using model defined in 'objLMLike'")
        if(!inherits(LMlike, 'LMlike')) stop("'LMlike' must inherit from class 'LMlike'")
        ## update design matrix with possibly new/permuted colData
                                        ##obj <- update(LMlike, design=colData(sca))
        obj <- LMlike
       }
    
    ## avoiding repeated calls to the S4 object speeds calls on large sca
    ## due to overzealous copying semantics on R's part
    ee <- exprs(sca)
    genes <- colnames(ee)
    ng <- length(genes)
    MM <- model.matrix(obj)
    coefNames <- colnames(MM)
    ## to facilitate our call to mclapply
    listEE <- setNames(seq_len(ng), genes)
    ## in hopes of finding a typical gene
    upperQgene <- which(rank(freq(sca), ties.method='random')==floor(.75*ng))
    obj <- fit(obj, ee[,upperQgene], silent=silent)

    ## called internally to do fitting, but want to get local variables in scope of function
    nerror <- 0
    .fitGeneSet <- function(idx){
        ## initialize outputs
        hookOut <- NULL
        tt <- try({
            obj <- fit(obj, response=ee[,idx], silent=silent, quick=TRUE)
            if(!is.null(hook)) hookOut <- hook(obj)
            nerror <- 0
            if((idx %% 20)==0) message('.', appendLF=FALSE)
        })

        if(is(tt, 'try-error')){
            obj@fitC <- obj@fitD <- NULL
            obj@fitted <- c(C=FALSE, D=FALSE)
            message('!', appendLF=FALSE)
            nerror <- nerror+1
            if(nerror>5 & !force) {
                stop("We seem to be having a lot of problems here...are your tests specified correctly?  \n If you're sure, set force=TRUE.", tt)                
            }
        }
        if(onlyCoef) return(cbind(C=coef(obj, 'C'), D=coef(obj, 'D')))
        summaries <- summarize(obj)
        structure(summaries, hookOut=hookOut)
    }


    if(!parallel || getOption('mc.cores', 1L)==1){
        listOfSummaries <- lapply(listEE, .fitGeneSet)
    } else{
    listOfSummaries <- parallel::mclapply(listEE, .fitGeneSet, mc.preschedule=TRUE, mc.silent=silent)
}

    if(onlyCoef){
        out <- do.call(abind, c(listOfSummaries, rev.along=0))
        return(aperm(out, c(3,1,2)))
    }
        
    ## test for try-errors
    cls <- sapply(listOfSummaries, function(x) class(x))
    complain <- if(force) warning else stop
    if(mean(cls=='try-error')>.5) complain('Lots of errors here..something is amiss.')

    ## gethooks
    hookOut <- NULL
    if(!is.null(hook)) hookOut <- lapply(listOfSummaries, attr, which='hookOut')

    
    
    message('\nDone!')
    summaries <- collectSummaries(listOfSummaries)

    ## add rest of slots, plus class name
    summaries[['LMlike']] <- obj
    summaries[['sca']] <- sca
    summaries[['priorVar']] <- obj@priorVar
    summaries[['priorDOF']] <- obj@priorDOF
    summaries[['hookOut']] <- hookOut
    summaries[['Class']] <- 'ZlmFit'
    ## everything we need to call new
    zfit <- do.call(new, as.list(summaries))
    ## tests, summarized objects, example fit, hooks
    zfit
}
