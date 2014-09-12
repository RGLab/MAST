methodDict <- data.table(keyword=c('glm', 'glmer', 'lmer', 'bayesglm'),
                         lmMethod=c('GLMlike', 'LMERlike','LMERlike', 'BayesGLMlike'),
                         lrtHypoType=c('SimpleHypothesis', 'TermHypothesis','TermHypothesis',  'SimpleHypothesis'),
                         implementsEbayes=c(TRUE, TRUE, FALSE, TRUE))


residualsHook <- function(fit){
    residuals(fit, which='Marginal')
}

revealHook <- function(zlm){
    return(attr(zlm, 'hookOut'))
}

##' @importFrom plyr laply
collectResiduals <- function(zlm, sca, newLayerName='Residuals'){
    if(newLayerName %in% dimnames(sca)[[3]]) warning('Overwriting layer', newLayerName) else     sca <- addlayer(sca, newLayerName)
    layer(sca) <- newLayerName
    mat <- t(laply(revealHook(zlm), function(x) x))
    exprs(sca) <- mat
    sca
}

##' Convenience function for running a zero-inflated regression
##'
##' Fits a hurdle model on zero-inflated continuous data in which the zero process
##' is modeled as a logistic regression
##' and (conditional on the the response being >0), the continuous process is Gaussian, ie, a linear regression.
##' @param formula model formula
##' @param data a data.frame, list or environment in which formula is evaluated
##' @param method one of 'glm' or 'glmer'.  See SingleCellAssay:::methodDict for other possibilities.
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
##' summary(fit$disc)
##' summary(fit$cont)
##'
##' @seealso GLMlike, LMERlike
zlm <- function(formula, data, method='glm',silent=TRUE, ...){
  if(!inherits(data, 'data.frame')) stop("'data' must be data.frame, not matrix or array")
  if(!is(formula, 'formula')) stop("'formula' must be class 'formula'")

  ## lm initially just to get response vector
  ## Turn glmer grouping "|" into "+" to get correct model frame
  resp <- eval(formula[[2]], data)
  obj <- new(methodDict[keyword==method, lmMethod], formula=formula, design=data, response=resp)
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
##' Conducts tests specified in "hypothesis".
##' After each gene, optionally run the function on the fit named by 'hook'
##'
##' When keep.zlm is FALSE, a 3D array with first dimension being the genes,
##' next dimension giving information about the test
##' (the degrees of freedom, Chisq statistic, and P value), and final dimension
##' being the value of these quantities on the
##' discrete, continuous and hurdle (combined) levels.
##'
##' When keep.zlm is TRUE, a list of length two is returned.
##' Component "tests" gives the above 3-D array.
##' Component "models" is a list giving the model fit for each gene.
##'
##' When \code{hypothesis} is a list, then each test specified in the list will be run, and the returned object will also be a list of 3D arrays (or component "tests" will be a list if keep.zlm is TRUE).
##'
##' The empirical bayes regularization of the gene variance assumes that the precision (1/variance) is drawn from a
##' gamma distribution with unknown parameters.
##' These parameters are estimated by considering the distribution of sample variances over all genes.
##' The procedure used for this is determined from
##' \code{ebayesControl}, a named list with components 'method' (one of 'MOM' or 'MLE') and 'model' (one of 'H0' or 'H1')
##' method MOM uses a method-of-moments estimator, while MLE using the marginal likelihood.
##' H0 model estimates the precisions using the intercept alone in each gene, while H1 fits the full model specified by \code{formula}
##'
##' @param formula a formula with the measurement variable on the LHS and predictors present in cData on the RHS
##' @param sca SingleCellAssay object
##' @param method character vector, either 'glm' or 'glmer'
##' @param hypothesis character vector, list of character vectors passed to \code{lrTest} or \code{waldTest}.  See details.
##' @param type type of test to run, one of 'Wald' or 'LRT'
##' @param onlyReturnCoefs if TRUE, don't actually test, only return a gene giving example coefficients
##' @param keep.zlm should the model objects be returned?  May be memory intensive.
##' @param .parallel currently ignored
##' @param silent Silence common problems with fitting some genes
##' @param ebayes if TRUE, regularize variance using empirical bayes method
##' @param ebayesControl list with parameters for empirical bayes procedure.  See \link{ebayes}.
##' @param hook a function called on the \code{fit} after each gene.
##' @param force Should we continue testing genes even after many errors have occurred?
##' @param ... arguments passed fit method.  For example, \code{fitArgsC} and \cpde{fitArgsD}.  These are a list of arguments passed to the underlying modeling functions.
##' @return either an array of tests (one per primer), a list of such arrays (one per hypothesis),  or a list with components "models" and "fits".
##' @export
##' @seealso ebayes, glmlike-class
##' @importFrom stringr str_split_fixed
##' @importFrom stringr fixed
##' @examples
##' \dontrun{
##' data(vbetaFA)
##' testsByGene <- zlm.SingleCellAssay(~ Stim.Condition, vbetaFA, hypothesis='Stim.ConditionUnstim', method='glm', type='Wald')
##' # genes X metric X test type
##' dimnames(testsByGene)
##'
##' modelsAndTestsByGene <- zlm.SingleCellAssay(~ Stim.Condition, vbeta.sc, hypothesis='Stim.ConditionUnstim', keep.zlm=TRUE)
##' names(modelsAndTestsByGene$models)
##' summary(modelsAndTestsByGene$models[['IL13']]$disc)
##' summary(modelsAndTestsByGene$models[['IL13']]$cont)
##'
##' ## Separate tests that Stim.Condition=0 and the Intercept=0
##' ## (The second test doesn't make sense scientifically)
##' twoTests <- zlm.SingleCellAssay(~ Stim.Condition, vbeta.sc, hypothesis=list('Stim.ConditionUnstim', '(Intercept)'))
##' length(twoTests)
##' dimnames(twoTests[[1]])
##' }
zlm.SingleCellAssay <- function(formula, sca, method='glm', hypothesis, type='Wald', onlyReturnCoefs=FALSE, keep.zlm='false', .parallel=FALSE, silent=TRUE, ebayes=FALSE, ebayesControl=NULL, force=FALSE, hook=NULL, ...){
    ## Which class are we using for the fits...look it up by keyword
    method <- match.arg(method, methodDict[,keyword])
    method <- methodDict[keyword==method,lmMethod]
    type <- match.arg(type, c('LRT', 'Wald'))
    ## What sort of test will we be doing?
    if(type=='LRT'){
        test <- lrTest        
    }else{
        test <- waldTest
    }
        if(!is(sca, 'SingleCellAssay')) stop("'sca' must be (or inherit) 'SingleCellAssay'")
    if(!is(formula, 'formula')) stop("'formula' must be class 'formula'")
    fsplit <- str_split_fixed(deparse(formula), fixed('~'), 2)
    if(nchar(fsplit[1,1])>0) message("Ignoring LHS of formula (", fsplit[1,1], ') and using exprs(sca)')
    Formula <- as.formula(paste0('~', fsplit[1,2]))

    ## Empirical bayes method
    priorVar <- 1
    priorDOF <- 0
    if(ebayes){
        if(!methodDict[method,'implementsEbayes']) stop('Method', method, ' does not implement empirical bayes variance shrinkage.')
        ebparm <- ebayes(sca, ebayesControl, Formula)
        priorVar <- ebparm['v']
        priorDOF <- ebparm['df']
    }
    ## initial value of priorVar, priorDOF default to no shrinkage
    obj <- new(method, design=cData(sca), formula=Formula, priorVar=priorVar, priorDOF=priorDOF)
    

    ## always set hypothesis to be enclosed in a list
    if(!inherits(hypothesis, 'list'))
        hypothesis <- list(hypothesis)
    nhypo <- length(hypothesis)
    ## avoiding repeated calls to the S4 object speeds calls on large sca
    ## due to overzealous copying semantics on R's part
    ee <- exprs(sca)    
    genes <- colnames(ee)
    ng <- length(genes)
    ## in hopes of finding a typical gene to get coefficients
    upperQgene <- which(rank(freq(sca), ties='random')==floor(.75*ng))
    obj <- fit(obj, ee[,upperQgene], silent=silent, ...)
    if(onlyReturnCoefs){
        print(show(obj))
        return(invisible(obj))
        }

    ## gets hypothesis in suitable form for testing by now comparing it to the model names
    ## this will throw an error if there are some columns misnamed
    MM <- model.matrix(obj)
    if(listType(hypothesis) %in% c('CoefficientHypothesis','Hypothesis')){
        hypothesis <- lapply(hypothesis, generateHypothesis, terms=colnames(MM))
    } else if(listType(hypothesis) != 'character'){
        stop("hypothesis must be 'character', 'CoefficientHypothesis' or 'Hypothesis'.")
    }

    ## if(qr(MM)$rank<ncol(MM)) warning('Rank deficient design may make some methods unhappy')

    ## initialize junk
    testNames <- makeChiSqTable(c(0, 0), c(1, 1), '')
    vcovNames <- coefNames <- colnames(MM)

    ltests <- setNames(vector(mode='list', length=nhypo), names(hypothesis))
    hookOut <- if(!is.null(hook)) setNames(vector(mode='list', length=ng), genes) else NULL
    for(h in seq_len(nhypo)){
        ltests[[h]] <- array(0, dim=c(ng, nrow(testNames), ncol(testNames)), dimnames=list(primerid=genes, test.type=row.names(testNames), metric=colnames(testNames)))
        ltests[[h]][,,'Pr(>Chisq)'] <- 1
}

    ## Main loop.
    ## Todo: coefs, vcov, etc
    ## error counter--stop if exceeds 5 in a row
     nerror <- 0
 innerCatch <- ''
    for(i in seq_len(ng)){
        outerCatch <- try({
            obj <- fit(obj, response=ee[,i], silent=silent, ...)
            if(!is.null(hook)) hookOut[[i]] <- hook(obj)
            for(h in seq_len(nhypo)){
                 innerCatch <- try({ltests[[h]][i,,] <- test(obj, hypothesis[[h]])}, silent=silent)
            }  
        }, silent=silent)
        if(is(outerCatch, 'try-error') || is(innerCatch, 'try-error')){
            message('!', appendLF=FALSE)
            nerror <- nerror+1
            if(nerror>5 & !force) {
                msg <- c(outerCatch, innerCatch)[c(is(outerCatch, 'try-error'), is(innerCatch, 'try-error'))]
                stop("We seem to be having a lot of problems here...are your tests specified correctly?  \n If you're sure, set force=TRUE.", msg)
                }
            next
        }     ## Made it through, reset error counter
        nerror <- 0
        message('.', appendLF=FALSE)
    }
    message('\nDone!')
    if(length(ltests)==1) ltests <- ltests[[1]]
    structure(ltests, obj=obj, hookOut=hookOut)
    
}
