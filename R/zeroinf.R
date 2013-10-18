##' Run a zero-inflated regression
##'
##' Fits a hurdle model on zero-inflated continuous data in which the zero process
##' is modeled as a logistic regression
##' and (conditional on the the response being >0), the continuous process is Gaussian, ie, a linear regression.
##' @param formula model formula
##' @param data a data.frame, list or environment in which formula is evaluated
##' @param lm.fun a function that takes a formula and data the arguments family='binomial' and family='gaussian', eg, \code{glm} or \code{glmer}.
##' @param silent if TRUE suppress common errors from fitting continuous part
##' @param subset ignored
##' @param ... passed to lm.fun
##' @return list of class 'zlm' with "disc"rete part and "cont"inuous part 
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
##' #Accepts arguments like drop1 and compares likelihood ratios
##' test.zlm(fit, type='LRT', hypothesis.matrix='x')
##' #Accepts arguments like car::linearHypothesisTest
##' test.zlm(fit, type='Wald', hypothesis.matrix=c('x=2', 'z=2'))
##' #Little evidence for diffence in discrete, big evidence in continuous
##' @importFrom stringr str_detect
zlm <- function(formula, data,lm.fun=glm,silent=TRUE, subset, ...){
  #if(!inherits(data, 'data.frame')) stop("'data' must be data.frame, not matrix or array")
  if(!missing(subset)) warning('subset ignored')
  if(!inherits(formula, 'formula')) stop("'formula' must be class 'formula'")

  ## lm initially just to get pull response vector
  ## Turn glmer grouping "|" into "+" to get correct model frame
  sanitize.formula <- as.formula(gsub('[|]', '+', deparse(formula)))
  ## Throw error on NA, because otherwise the next line will fail mysteriously
  init <- tryCatch(model.frame(sanitize.formula, data, na.action=na.fail), error=function(e) if(str_detect(as.character(e), 'missing')) stop('NAs in response or predictors not allowed; please remove before fitting') else stop(e) )
  
  data[,'pos'] <-   model.response(init)>0
  cont <- try(lm.fun(formula, data, subset=pos, family='gaussian', ...), silent=silent)
  if(inherits(cont, 'try-error')){
    if(!silent) warning('Some factors were not present among the positive part')
    cont <- lm(0~0)
  }                                     
  
  formula.split <- strsplit(deparse(formula), '~')[[1]]
  lhs <- formula.split[1]
  if(str_detect(lhs, '[()]')) stop("Left hand side of formula must be unadorned variable name from 'data'")
  lhs <- 'pos'
  rhs <- formula.split[2]
  disc.formula <- paste(lhs, '~', rhs)
  disc <- lm.fun(disc.formula, data, family='binomial', ...)
  
  out <- list(cont=cont, disc=disc)
  class(out) <- 'zlm'
  out
}

is.empty.fit <- function(fit) return(length(coef(fit))==0 || summary(fit)$df.residual==0)

summary.zlm <- function(out){
  summary(out$cont)
  summary(out$disc)
}

##' Likelihood ratio test for hurdle model
##'
##' Do LR test separately on continuous and discrete portions
##' Combine for testing hypothesis.matrix
##'
##' This just internally calls lht from package car on the discrete and continuous models.
##' It tests the provided hypothesis.matrix using a Chi-Squared 
##' @param model output from zlm
##' @param hypothesis.matrix argument passed to lht, or string naming a single variable to be dropped from the model, or a formula that is a subset of 'model'
##' @param type Test using Wald test or Likelihood Ratio test
##' @param silent Silence common errors in testing
##' @return array containing the discrete, continuous and combined tests
##' @importFrom car linearHypothesis lht matchCoef
##' @seealso zlm
##' @export
test.zlm <- function(model, hypothesis.matrix, type='Wald', silent=TRUE){
    type <- match.arg(type, c('Wald', 'LRT'))
    if(length(type)!= 1 || (type != 'Wald'&& type != 'LRT')) stop("'type' must equal 'Wald' or 'LRT'")
    if(type=='LRT'){
        if(length(hypothesis.matrix) > 1 && !is.formula(hypothesis.matrix)) stop("Currently only support testing single factors when 'type'='LRT' and length of 'hypothesis.matrix' > 1")
        if(!inherits(model$disc, 'lm')) stop('Currently only support type=LRT with glm fits')
    }
    
  mer.variant <- any('chisq' %in% eval(formals(getS3method('linearHypothesis', class(model$disc)))$test)) #don't ask
    ## Get names to agree from output of all the different variants
  chisq <- 'Chisq'
  pchisq <- 'Pr(>Chisq)'
  names.drop1.cont <- c('Df', 'scaled dev.', 'Pr(>Chi)')
  names.drop1.disc <- c('Df', 'LRT', 'Pr(>Chi)')
  rename.drop1.cont <- c('scaled dev.'='Chisq', 'Pr(>Chi)'='Pr(>Chisq)')
  rename.drop1.disc <- c('LRT'='Chisq', 'Pr(>Chi)'='Pr(>Chisq)')
  if(mer.variant) {
    chisq <- 'chisq'
    pchisq <- 'Pr(> Chisq)'
  }
  if(type=='Wald'){
  tt <- try({
    cont <- lht(model$cont, hypothesis.matrix, test=chisq, singular.ok=TRUE)
  }, silent=silent)
  disc <- lht(model$disc, hypothesis.matrix, test=chisq, singular.ok=TRUE)
} else if(type=='LRT'){
    tt <- try({
    if(summary(model$cont)$df.residual==0) stop('No degrees of freedom left') #otherwise drop1 throws an obscure error
    cont <- rename(
        cbind(Res.df=NA, drop1(model$cont, hypothesis.matrix, test='LRT')[, names.drop1.cont]),
        rename.drop1.cont)
}, silent=silent)
        disc <- rename(
        cbind(Res.df=NA, drop1(model$disc, hypothesis.matrix, test='LRT')[, names.drop1.disc]),
        rename.drop1.disc)
    } else{
 stop('ruhroh')
}

    
  if(inherits(tt, 'try-error') || !all(dim(cont) == dim(disc))){
    cont <- rep(0, length(as.matrix(disc)))
    dim(cont) <- dim(disc)
    dimnames(cont) <- dimnames(disc)
    cont[,pchisq] <- NA 
  }

  res <- abind(disc, cont, disc+cont, rev.along=0)
  dimnames(res)[[3]] <- c('disc', 'cont', 'hurdle')
  res[,pchisq,3] <- sapply(seq_len(nrow(disc)), function(i)
                                 pchisq(res[i,'Chisq',3], df=res[i, 'Df', 3], lower.tail=FALSE))
      dm <- dimnames(res)
      names(dm) <- c('', 'metric', 'test.type')
  if(mer.variant) {
    dm[[2]][dm[[2]]==pchisq] <- 'Pr(>Chisq)'
  }
      dimnames(res) <- dm
  res
}

##' zero-inflated regression for SingleCellAssay 
##'
##' For each gene in sca, fits the hurdle model in \code{formula} (linear for et>0), logistic for et==0 vs et>0.
##' Conducts tests using hypothesis.matrix.
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
##' @title zlm.SingleCellAssay
##' @param formula a formula with the measurement variable on the LHS and predictors present in cData on the RHS
##' @param sca SingleCellAssay object
##' @param lm.fun a function accepting lm-style arguments and a family argument
##' @param hypothesis.matrix names of coefficients to test in lht form if type='Wald', otherwise a character vector if type ='LRT'.
##' @param type type of test to run, one of 'Wald' or 'LRT'
##' @param hypo.fun a function taking a model as input and returning output suitable for hypothesis.matrix
##' @param keep.zlm should the model objects be kept
##' @param .parallel run fits using parallel processing.  must have doParallel
##' @param .drop see ldply
##' @param .inform see ldply
##' @param silent Silence common problems with fitting some genes
##' @param ... passed to lm.fun
##' @return either an array of tests (one per primer) or a list
##' @export
##' @importFrom car lht
##' @importFrom plyr laply
##' @importFrom plyr llply
##' @importFrom plyr dlply
##' @seealso zlm
##' @seealso test.zlm
##' @examples
##' \dontrun{
##' data(vbeta)
##' vbeta <- computeEtFromCt(vbeta)
##' vbeta.sc <- FluidigmAssay(vbeta, idvars='Sample.ID', primerid='Gene', measurement='Et', ncells='Number.of.Cells', cellvars='Stim.Condition')
##' testsByGene <- zlm.SingleCellAssay(Et ~ Stim.Condition, vbeta.sc, hypothesis.matrix='Stim.ConditionUnstim')
##' # genes X metric X test type
##' dimnames(testsByGene)
##'
##' modelsAndTestsByGene <- zlm.SingleCellAssay(Et ~ Stim.Condition, vbeta.sc, hypothesis.matrix='Stim.ConditionUnstim', keep.zlm=TRUE)
##' names(modelsAndTestsByGene$models)
##' summary(modelsAndTestsByGene$models[['IL13']]$disc)
##' summary(modelsAndTestsByGene$models[['IL13']]$cont)
##' }
zlm.SingleCellAssay <- function(formula, sca, lm.fun=glm, hypothesis.matrix, type='Wald', hypo.fun=NULL, keep.zlm=FALSE, .parallel=FALSE, .drop=TRUE, .inform=FALSE, silent=TRUE, ...){

  
    m <- SingleCellAssay:::melt(sca)

    if(.drop) m <- droplevels(m)

    models <- dlply(m, ~primerid, function(set){
        zlm(formula, set, lm.fun, silent, ...)
        }, .drop=.drop, .inform=.inform, .parallel=.parallel)

    fit.primerid <- function(model){
            if(!is.null(hypo.fun) && inherits(hypo.fun, 'function')){
              hypothesis.matrix <- hypo.fun(model)
            }
            test.zlm(model, hypothesis.matrix, type=type, silent=silent)
    }
    
    test.models <- llply(models, fit.primerid, .inform=.inform)

    tests <- laply(test.models, function(x){
      x[2,-1,]
    }, .inform=.inform)

    if(keep.zlm){
      return(list(tests=tests, models=models))
    }

    return(tests)
}
