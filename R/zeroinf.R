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
##' #Compares likelihood ratios, currently only can drop an entire term
##' test.zlm(fit, type='LRT', hypothesis.matrix='x')
##' #Accepts arguments like car::linearHypothesisTest
##' test.zlm(fit, type='Wald', hypothesis.matrix=c('x=2', 'z=2'))
##' #Little evidence for difference in discrete, big evidence in continuous
##' @importFrom stringr str_detect
zlm <- function(formula, data,lm.fun=glm,silent=TRUE, subset, ...){
  #if(!inherits(data, 'data.frame')) stop("'data' must be data.frame, not matrix or array")
  if(!missing(subset)) warning('subset ignored')
  if(!inherits(formula, 'formula')) stop("'formula' must be class 'formula'")

  ## lm initially just to get response vector
  ## Turn glmer grouping "|" into "+" to get correct model frame
  sanitize.formula <- as.formula(gsub('[|]', '+', deparse(formula)))
  ## Throw error on NA, because otherwise the next line will fail mysteriously
  init <- tryCatch(model.frame(sanitize.formula, data, na.action=na.fail), error=function(e) if(str_detect(as.character(e), 'missing')) stop('NAs in response or predictors not allowed; please remove before fitting') else stop(e) )
  
  data[,'pos'] <-   model.response(init)>0
  dp <- data[data$pos,]

  ## new glmer doesn't like being called with a function alias
  useGlmer <-  exists('glmer') && identical(lm.fun, glmer)
  
  cont <- try({
      if(useGlmer){
      lmer(formula, dp, ...)
      } else{
      lm.fun(formula, data, subset=pos, family='gaussian', ...)
      }
  }, silent=silent)
  if(inherits(cont, 'try-error')){
    if(!silent) warning('Some factors were not present among the positive part')
    cont <- lm(0~0)
  }                                     

  converged <- TRUE
  formula.split <- strsplit(deparse(formula), '~')[[1]]
  lhs <- formula.split[1]
  if(str_detect(lhs, '[()]')) stop("Left hand side of formula must be unadorned variable name from 'data'")
  lhs <- 'pos'
  rhs <- formula.split[2]
  disc.formula <- paste(lhs, '~', rhs)
  if(useGlmer){
      disc <- try({
       glmer(disc.formula, data, family='binomial', ...)
  }, silent=silent)
      } else{
          disc <- try({
              lm.fun(disc.formula, data, family='binomial', ...)
  }, silent=silent)
      }
  if(inherits(disc, 'try-error')){
      disc <- lm(0~0)
      converged <- FALSE
  }
  
  
  out <- list(cont=cont, disc=disc, converged=converged)
  class(out) <- 'zlm'
  out
}

is.empty.fit <- function(fit) return(length(coef(fit))==0 || (inherits(fit, 'lm') && summary(fit)$df.residual==0))

summary.zlm <- function(out){
  summary(out$cont)
  summary(out$disc)
}

## Because lmer reports the posterior mode estimates for each 
coefFun <- function(model){
    if(inherits(model, 'glmerMod')){
        return(fixef(model))
    }else{
        return(coef(model))
    }
    
}





pretest.lrt <- function(model, hypothesis.matrix){
    Terms <- labels(terms(model))
    whichTerm <- attr(model.matrix(model), 'assign')
    names(whichTerm) <- names(coefFun(model))
    ## terms that were tested
    testTerm <- whichTerm[hypothesis.matrix]
    ## How many of each that were fit were tested?
    allTerm <- ifelse(table(testTerm)>0, table(testTerm)==table(testTerm), TRUE)
    if(any(!allTerm)) stop('Currently must test an entire term at once when type="LRT"')
     drop.terms <- as.formula(paste('~', paste(Terms[unique(testTerm)], collapse='+')))
}

test.lrt <- function(model, drop.terms, part, ...){
    if(!is.formula(drop.terms) || length(labels(terms(drop.terms)))  > 1) stop("Currently only support testing single factors when 'type'='LRT'")
    newlme4 <- inherits(model, 'glmerMod') | inherits(model, 'lmerMod')
    ## No longer true with lme4 > 1.0
    if(!(inherits(model, 'lm') || newlme4)) stop('Currently only support type=LRT with glm fits') 
    if(part=='cont' && inherits(model, 'lm') && summary(model)$df.residual==0) stop('No degrees of freedom left') #otherwise drop1 throws an obscure error


    
    if(newlme4){
        names.drop1.disc <- c('Df', 'LRT', 'Pr(Chi)')
        rename.drop1.disc <- c('LRT'='Chisq', 'Pr(Chi)'='Pr(>Chisq)')
        names.drop1.cont <- c('Df', 'LRT', 'Pr(>Chi)')
    rename.drop1.cont <- c('LRT'='Chisq', 'Pr(>Chi)'='Pr(>Chisq)')

    } else{
        names.drop1.disc <- c('Df', 'LRT', 'Pr(>Chi)')
        rename.drop1.disc <- c('LRT'='Chisq', 'Pr(>Chi)'='Pr(>Chisq)')
          names.drop1.cont <- c('Df', 'scaled dev.', 'Pr(>Chi)')
    rename.drop1.cont <- c('scaled dev.'='Chisq', 'Pr(>Chi)'='Pr(>Chisq)')

    }

    if(part=='cont'){
        return(rename(
            cbind(Res.df=NA, drop1(model, drop.terms, test='Chisq', ...)[, names.drop1.cont]),
            rename.drop1.cont))
    } else{
        return(rename(
            cbind(Res.df=NA, drop1(model, drop.terms, test='Chisq', ...)[, names.drop1.disc]),
            rename.drop1.disc))
    }
}

test.wald <- function(model, hypothesis.matrix, part){
    cls<-class(model)
    if(any(cls%in%c("glmerMod","lmerMod")))
      cls<-"mer"
    
      mer.variant <- any('chisq' %in% (eval(formals(getS3method('linearHypothesis', cls))$test))) #don't ask
    ## Get names to agree from output of all the different variants
  chisq <- 'Chisq'
  pchisq <- 'Pr(>Chisq)'
        if(mer.variant) {
    chisq <- 'chisq'
    pchisq <- 'Pr(> Chisq)'
  } 
    lh.out <- lht(model, hypothesis.matrix, test=chisq, singular.ok=TRUE)
           if(mer.variant) {
            lh.out <- rename(lh.out, c('chisq'='Chisq', 'Pr(> Chisq)'='Pr(>Chisq)'))
  }
      lh.out


}

##' Likelihood ratio test for hurdle model
##'
##' Do LR test separately on continuous and discrete portions
##' Combine for testing hypothesis.matrix
##'
##' It tests the provided hypothesis.matrix using a Chi-Squared Wald or LRT test
##' @param model output from zlm
##' @param hypothesis.matrix argument passed to lht, or string naming a single variable to be dropped from the model, or a formula that is a subset of 'model'
##' @param type Test using Wald test or Likelihood Ratio test
##' @param silent Silence common errors in testing
##' @return array containing the discrete, continuous and combined tests
##' @importFrom car linearHypothesis lht matchCoefs
##' @seealso zlm
##' @export
test.zlm <- function(model, hypothesis.matrix, type='Wald', silent=TRUE, ...){
    type <- match.arg(type, c('Wald', 'LRT'))
    if(length(type)!= 1 || (type != 'Wald'&& type != 'LRT')) stop("'type' must equal 'Wald' or 'LRT'")
    
  if(type=='Wald'){
  tt <- try({
    cont <- test.wald(model$cont, hypothesis.matrix, part='cont')
  }, silent=silent)
  disc <- test.wald(model$disc, hypothesis.matrix, part='disc')
} else if(type=='LRT'){
    drop.terms <- pretest.lrt(model$disc, hypothesis.matrix)
    tt <- try({
    cont <- test.lrt(model$cont, drop.terms, part='cont', ...)
}, silent=silent)
    disc <- test.lrt(model$disc, drop.terms, part='disc', ...)
    } else{
 stop('ruhroh')
}
    
  if(inherits(tt, 'try-error') || !all(dim(cont) == dim(disc))){
    cont <- rep(0, length(as.matrix(disc)))
    dim(cont) <- dim(disc)
    dimnames(cont) <- dimnames(disc)
    cont[,'Pr(>Chisq)'] <- NA 
  }

  res <- abind(disc, cont, disc+cont, rev.along=0)
  dimnames(res)[[3]] <- c('disc', 'cont', 'hurdle')
  res[,'Pr(>Chisq)',3] <- sapply(seq_len(nrow(disc)), function(i)
                                 pchisq(res[i,'Chisq',3], df=res[i, 'Df', 3], lower.tail=FALSE))
      dm <- dimnames(res)
      names(dm) <- c('', 'metric', 'test.type')
      dimnames(res) <- dm
  res
}

##' try to guess form of LHT
##'
##' @param hypo.terms 
##' @param model 
##' @return character vector of terms to test
##' @importFrom stringr fixed
guessContrast <- function(hypo.terms, model){
    stopifnot(is.character(hypo.terms))
    coef.names <- names(coefFun(model$disc))
            detected <- laply(hypo.terms, function(term) str_detect(coef.names, fixed(term)), .drop=FALSE)
    anyDetected <- apply(detected, 2, any)
            message(paste("Testing coefficients" , paste(hypo.terms, collapse=', ')))
    if(any(str_detect(hypo.terms, '[+=-]'))) warning("'+', '-' or '=' found in variable names, contrasts may not work as intended")
    if(all(!anyDetected)) stop(sprintf("terms %s did not match any coefficients in %s", paste(hypo.terms, collapse=' ,'), paste(coef.names, collapse=' ,')))
    coef.names[anyDetected]
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
##' @param hypothesis.matrix Deprecated
##' @param hypo.terms character vector giving terms to drop from model
##' @param hypo.contrasts specific contrasts to test in form expected by lht
##' @param type type of test to run, one of 'Wald' or 'LRT'
##' @param keep.zlm should the model objects be kept?
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
zlm.SingleCellAssay <- function(formula, sca, lm.fun=glm, hypothesis.matrix, hypo.terms, hypo.contrasts, type='Wald', keep.zlm=FALSE, .parallel=FALSE, .drop=TRUE, .inform=FALSE, silent=TRUE, ...){
    nmiss <- missing(hypo.contrasts)*1 +missing(hypo.terms)*1 + missing(hypothesis.matrix)*1
    if(nmiss != 2) stop("Specify one and only one of 'hypothesis.matrix', 'hypo.terms', or 'hypo.contrasts'")

     hypo.contrasts.missing <- missing(hypo.contrasts)
    if(!missing(hypothesis.matrix)){
        message("'hypothesis.matrix' is deprecated, use 'hypo.terms' or 'hypo.contrasts'")
        if(type=='LRT') hypo.terms <- hypothesis.matrix
        if(type=='Wald') {
            hypo.contrasts <- hypothesis.matrix
            hypo.contrasts.missing <- FALSE
        }
   }
    
    m <- SingleCellAssay:::melt(sca)
    if(.drop) m <- droplevels(m)

    ##keeping legacy support for logicals
    keep.zlm <- casefold(as.character(keep.zlm))  
    keep.zlm <- match.arg(keep.zlm, c('true', 'false', 'coefficients'))

    
    
    fit.primerid <- function(set, ...){
        model <- zlm(formula, set, lm.fun, silent, ...)

        if(hypo.contrasts.missing) hypo.contrasts <- guessContrast(hypo.terms, model)
        if(model$converged){
        test <- test.zlm(model, hypo.contrasts, type=type, silent=silent, ...)[2,-1,]
    } else{
        ## hack, but this will be a pain to do correctly
        test <- rep(NA, 9)
        dim(test) <- c(3, 3)
        dimnames(test) <- list(metric=c('Df', 'Chisq', 'Pr(>Chisq)'), test.type=c('disc', 'cont', 'hurdle'))
    }
        switch(keep.zlm,
               true=list(model=model, test=test),
               false=list(test=test),              #todo: write coefs function
               coefficients=list(coef=coefs(model), test=test)
               )
    }

    geneTests <- dlply(m, 'primerid', fit.primerid, .parallel=.parallel, .drop=.drop, .inform=.inform, ...)
    tests <- laply(geneTests, function(x) x$test)
    if(keep.zlm=='true') return(list(models=llply(geneTests, function(x) x$model), tests=tests))
    if(keep.zlm=='coefficients') return(list(coefs=llply(geneTests, function(x) x$coef), tests=tests))
    return(tests)
}
