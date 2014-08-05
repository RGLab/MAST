methodDict <- c('glm'='GLMlike', 'glmer'='LMERlike', 'lmer'='LMERlike', 'bayesglm'='BayesGLMlike', 'shrunkglm'='ShrunkenGLMlike')

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
  obj <- new(methodDict[method], formula=formula, design=data, response=resp)
  obj <- fit(obj)
  list(cont=obj@fitC, disc=obj@fitD)
}

summary.zlm <- function(out){
  summary(out$cont)
  summary(out$disc)
}


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
        cat('Start:', date(), '\n')
        mm <- model.matrix(Formula, cData(sca))

        allfits <- getSSg_rNg(sca, mm)
        valid <- apply(!is.na(allfits), 1, all) & allfits[, 'rNg']/allfits[, 'SSg']<truncate
        valid[is.na(valid)] <- FALSE
        SSg <- allfits[valid,'SSg']
        rNg <- allfits[valid, 'rNg']
        cat('End:', date(), '\n')
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

 

##' zero-inflated regression for SingleCellAssay 
##'
##' For each gene in sca, fits the hurdle model in \code{formula} (linear for et>0), logistic for et==0 vs et>0.
##' Conducts tests specified in "hypothesis".
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
##' @param hypothesis character vector or list of character vectors passed to \code{lrTest} or \code{waldTest}.  See details.
##' @param type type of test to run, one of 'Wald' or 'LRT'
##' @param keep.zlm should the model objects be returned?  May be memory intensive.
##' @param .parallel currently ignored
##' @param silent Silence common problems with fitting some genes
##' @param ebayes if TRUE, regularize variance using empirical bayes method
##' @param ebayesControl list with parameters for empirical bayes procedure.  See details.
##' @param ... 
##' @param force Should we continue testing genes even after many errors have occurred?
##' @param onlyReturnCoefs if TRUE, don't actually test, only return a gene giving example coefficients
##' @param fitArgs list of arguments passed to glm/glmer
##' @return either an array of tests (one per primer), a list of such arrays (one per hypothesis),  or a list with components "models" and "fits".
##' @export
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
zlm.SingleCellAssay <- function(formula, sca, method='glm', hypothesis, type='Wald', onlyReturnCoefs=FALSE, keep.zlm='false', .parallel=FALSE, silent=TRUE, ebayes=FALSE, ebayesControl=NULL, force=FALSE, ...){
    method <- match.arg(method, names(methodDict))
    method <- methodDict[method]
    type <- match.arg(type, c('LRT', 'Wald'))
    test <- if(type=='LRT') lrTest else waldTest
        if(!is(sca, 'SingleCellAssay')) stop("'sca' must be (or inherit) 'SingleCellAssay'")
    if(!is(formula, 'formula')) stop("'formula' must be class 'formula'")
    fsplit <- str_split_fixed(deparse(formula), fixed('~'), 2)
    if(nchar(fsplit[1,1])>0) message("Ignoring LHS of formula (", fsplit[1,1], ') and using exprs(sca)')
    Formula <- as.formula(paste0('~', fsplit[1,2]))

    if(ebayes){
        ## Empirical bayes method
        if(method != 'ShrunkenGLMlike') warning('Selecting method "ShrunkenGLMlike" since ebayes=TRUE.')
        ebparm <- ebayes(sca, ebayesControl, Formula)
        obj <- new('ShrunkenGLMlike', design=cData(sca), formula=Formula, priorVar=ebparm['v'], priorDOF=ebparm['df'])
    } else{
        obj <- new(method, design=cData(sca), formula=Formula)
    }

    if(is.character(hypothesis)){
        hypothesis <- list(hypothesis)
    }
    nhypo <- length(hypothesis)
    
    genes <- colnames(exprs(sca))
    ng <- length(genes)
    upperQgene <- which(rank(freq(sca), ties='random')==floor(.75*ng))
    obj <- fit(obj, exprs(sca)[,upperQgene], silent=silent, ...)
    
    if(onlyReturnCoefs){
        print(summary(obj))
        invisible(obj)
        }

    testNames <- makeChiSqTable(c(0, 0), c(1, 1), '')
    coefNames <- names(coef(obj, 'C'))
    vcovNames <- colnames(vcov(obj, 'C'))
    ltests <- setNames(vector(mode='list', length=nhypo), names(hypothesis))
    for(h in seq_len(nhypo)){
        ltests[[h]] <- array(0, dim=c(ng, nrow(testNames), ncol(testNames)), dimnames=list(primerid=genes, test.type=row.names(testNames), metric=colnames(testNames)))
        ltests[[h]][,,'Pr(>Chisq)'] <- 1
}

    ## Main loop.  Not a very R-like expression, but prevents allocating a huge amount of memory in case we have many genes.
    ## Todo: coefs, vcov, etc
    ## error counter--stop if exceeds 5 in a row
     nerror <- 0
 innerCatch <- ''
    for(i in seq_len(ng)){
        outerCatch <- try({
            obj <- fit(obj, response=exprs(sca)[,i], silent=silent, ...)
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
    structure(ltests, obj=obj)
    
}
