##' @import Biobase
##' @importFrom plyr rbind.fill
##' @import methods
##' @include AllGenerics.R
NULL

##' MAST: Model-based Analysis of Single- cell Transcriptomics
##'
##' Methods for analysing single cell assay data using hurdle models.
##'
##' This packages provides data structures and functions for statistical analysis of single-cell assay data such as Fluidigm single cell gene expression assays.
##' @references Finak, et al.  MAST: a flexible statistical framework for assessing transcriptional changes and characterizing heterogeneity in single-cell RNA sequencing data.  Genome Biology (2015).
"_PACKAGE"

##' Vbeta Data Set
##' @docType data
##' @name vbeta
##' @rdname vbeta-dataset
##' @format a data frame with 11 columns.
##' Column \code{Ct} contains the cycle threshold, with NA denoting that the threshold was never crossed.  So it is inversely proportional to the log2 mRNA, and should be negated (and NAs set to zero) if it is used as a expression measurement for a \code{FluidigmAssay}.
NULL


##' Vbeta Data Set, FluidigmAssay
##' @docType data
##' @name vbetaFA
##' @rdname vbetaFA-dataset
##' @format a \code{FluidigmAssay} of the vbeta data set.
##' @seealso \code{\link{vbeta}}, \code{\link{FromFlatDF}}
NULL


##' MAITs data set, RNASeq
##' @docType data
##' @name maits
##' @rdname maits-dataset
##' @format a \code{list} containing an expression matrix (\code{expressionmat}), cell \code{cdat} and feature \code{fdat}.
##' @seealso \code{\link{FromMatrix}}
NULL


Mandatory_Featurevars <- character()
Mandatory_Cellvars <- character()

##' @import SummarizedExperiment
##' @import S4Vectors
##' @importMethodsFrom S4Vectors mcols
##' @importMethodsFrom SummarizedExperiment colData assays assay
setClass('SingleCellAssay', contains='SummarizedExperiment0',
         slots=list(cmap='character', fmap='character'),
         prototype=list(cmap=Mandatory_Cellvars,
                        fmap=Mandatory_Featurevars))

Fluidigm_Cellvars <- c(Mandatory_Cellvars, ncells='ncells')
setClass('FluidigmAssay', contains='SingleCellAssay', prototype=list(cmap=Fluidigm_Cellvars))



## Classes
##' Linear Model-like Class
##'
##' Wrapper around modeling function to make them behave enough alike that Wald tests and Likelihood ratio are easy to do.
##' To implement a new type of zero-inflated model, extend this class.
##' Depending on how different the method is, you will definitely need to override the \code{fit} method, and possibly the \code{model.matrix}, \code{model.matrix<-}, \code{update}, \code{coef}, \code{vcov}, and \code{logLik} methods.
##'
##' @section Slots:
##' \describe{
##' \item{design}{a data.frame from which variables are taken for the right hand side of the regression}
##' \item{fitC}{The continuous fit}
##' \item{fitD}{The discrete fit}
##' \item{response}{The left hand side of the regression}
##' \item{fitted}{A \code{logical} with components "C" and "D", TRUE if the respective component has converged}
##' \item{formula}{A \code{formula} for the regression}
##' \item{fitArgsC}{}
##' \item{fitArgsD}{Both \code{list}s giving arguments that will be passed to the fitter (such as convergence criteria or case weights)}
##' }
##' @seealso coef
##' @seealso lrTest
##' @seealso waldTest
##' @seealso vcov
##' @seealso logLik
setClass('LMlike',
         slots=c(design='ANY', modelMatrix='matrix', fitC='ANY', fitD='ANY', response='ANY', fitted='logical', formula='formula', fitArgsD='list', fitArgsC='list', priorVar='numeric', priorDOF='numeric',
                 ## this speeds construction of coef and vcov, which is a pinch point in zlm
                 defaultCoef='numeric',
                 defaultVcov='matrix'),
         prototype=list(fitted =c(C=FALSE, D=FALSE), formula=formula(0~0),modelMatrix=matrix(nrow=0, ncol=0), priorVar=0, priorDOF=0), validity=function(object){
             stopifnot( all(c("C", "D") %in% names(object@fitted)))
             if(length(object@response)>0){
                 if(any(is.na(object@response))) stop('NAs not permitted in response')
                 if(length(object@response)!=nrow(object@design)) stop('Response length differs from design length')
                 ##if(nrow(object@design) != nrow(object@modelMatrix)) stop('Design length differs from model.matrix length')
             }
         })

##' Wrapper for regular glm/lm
##'
##' @slot weightFun function to map expression values to probabilities of expression.  Currently unused.
setClass('GLMlike', contains='LMlike', slots=c(weightFun='function'), prototype=list(weightFun=function(x){
    ifelse(x>0, 1, 0)
}))

##' Initialize a prior to be used a prior for BayeGLMlike/BayesGLMlike2
##'
##' @param names character vector of coefficients.  The `(Intercept)` will be ignored.
##' @return 3d array, with leading dimension giving the prior 'loc'ation, 'scale' and degrees of freedom (df),
##' second dimension giving the component ('C'ontinuous or 'D'iscrete)
##' and trailing dimension giving the coefficient to which the prior applies.
##' The location is initialized to be 0, the scale to 2, and degrees of freedom of 1, following the default of bayesglm.
##' @export
##' @examples
##' dp <- defaultPrior('Stim.ConditionUnstim')
##' \dontrun{
##' data(vbetaFA)
##' zlmVbeta <- zlm.SingleCellAssay(~ Stim.Condition, vbeta.sc, method='bayesglm', coefPrior=dp)
##' }
defaultPrior <- function(names){
    names <- setdiff(names, '(Intercept)')
    p <- length(names)
    ar <- array(rep(c(0, 2.5, 1), times=2*p), dim=c(3, 2,p), dimnames=list(metric=c('loc', 'scale', 'df'), comp=c('C', 'D'), names))
                                        #if(p>0)     ar['scale',,names=='(Intercept)'] <- 10
    ar
}


##' Wrapper for bayesian GLM
##'
##' @slot prior \code{numeric} optional 3d array used to specify prior for coefficients
##' @slot useContinuousBayes \code{logical} should \code{bayesglm} be used to fit the continuous component as well?
setClass('BayesGLMlike', contains='GLMlike', slots=c(coefPrior='array', useContinuousBayes='logical'),
         prototype=list(prior=defaultPrior(character(0)), useContinuousBayes=FALSE),
         validity=function(object){
             ## if(length(object@coefPrior>0))
             ##     if(dim(object@coefPrior)[3] != sum(colnames(model.matrix(object))!='(Intercept)')) stop('prior must have same number of components as model.matrix')
             TRUE
         })
setClass('BayesGLMlikeWeight', contains='BayesGLMlike')


##' Wrapper for lmer/glmer
##'
##' A horrendous hack is employed in order to do arbitrary likelihood ratio tests: the model matrix is built, the names possibly mangled, then fed in as a symbolic formula to glmer/lmer.
##' This is necessary because there is no (easy) way to specify an arbitrary fixed-effect model matrix in glmer.
##' @slot pseudoMM part of this horrendous hack.
##' @slot strictConvergence \code{logical} return results even when the optimizer or *lmer complains about convergence
##' @slot optimMsg \code{character} record warnings from lme.  \code{NA_character_} means no warnings.
setClass('LMERlike', contains='LMlike', slots=c(pseudoMM='data.frame',
                                                optimMsg='character',
                                                strictConvergence='logical'),
         validity=function(object){
    if(length(object@response)>0 & nrow(object@pseudoMM)>0){
        stopifnot(nrow(object@pseudoMM)==length(object@response))
    }
    if(object@priorDOF!=0) stop('Empirical bayes shrinkage not implemented for lmer/glmer.')
    },
    prototype=list(strictConvergence=TRUE, optimMsg=c(C=NA_character_, D=NA_character_))
    )

setClass('bLMERlike', contains='LMERlike')

setClass('ConstrainedGLMlike', contains='LMlike')
setClass('RidgeBGLMlike',contains="BayesGLMlike",slots=c(lambda='numeric'),prototype = list(lambda=0.1) )

## Ways to specify hypothesis
setClass('Hypothesis', contains='character', slots=list(contrastMatrix='matrix'))
setClass('CoefficientHypothesis', contains='Hypothesis', slots=list(index='numeric'))

##' An S4 class to hold the output of a call to zlm
##'
##' This holds output from a call to zlm.SingleCellAssay.  Many methods are defined to operate on it.  See below.
##' @slot coefC matrix of continuous coefficients
##' @slot coefD matrix of discrete coefficients
##' @slot vcovC array of variance/covariance matrices for coefficients
##' @slot vcovD array of variance/covariance matrices for coefficients
##' @slot LMlike the LmWrapper object used
##' @slot sca the \code{SingleCellAssay} object used
##' @slot deviance matrix of deviances
##' @slot loglik matrix of loglikelihoods
##' @slot df.null matrix of null (intercept only) degrees of freedom
##' @slot df.resid matrix of residual DOF
##' @slot dispersion matrix of dispersions (after shrinkage)
##' @slot dispersionNoShrink matrix of dispersion (before shrinkage)
##' @slot priorDOF shrinkage weight in terms of number of psuedo-obs
##' @slot priorVar shrinkage target
##' @slot converged output that may optionally be set by the underlying modeling function
##' @slot hookOut a list of length ngenes containing output from a hook function, if \code{zlm} was called with one
##' @seealso zlm.SingleCellAssay summary,ZlmFit-method
##' @examples
##' data(vbetaFA)
##' zlmVbeta <- zlm.SingleCellAssay(~ Stim.Condition+Population, subset(vbetaFA, ncells==1)[1:10,])
##' #Coefficients and standard errors
##' coef(zlmVbeta, 'D')
##' coef(zlmVbeta, 'C')
##' se.coef(zlmVbeta, 'C')
##' #Test for a Population effect by dropping the whole term (a 5 degree of freedom test)
##' lrTest(zlmVbeta, 'Population')
##' #Test only if the VbetaResponsive cells differ from the baseline group
##' lrTest(zlmVbeta, CoefficientHypothesis('PopulationVbetaResponsive'))
##' # Test if there is a difference between CD154+/Unresponsive and CD154-/Unresponsive.
##' # Note that because we parse the expression
##' # the columns must be enclosed in backquotes
##' # to protect the \quote{+} and \quote{-} characters.
##' lrTest(zlmVbeta, Hypothesis('`PopulationCD154+VbetaUnresponsive` - `PopulationCD154-VbetaUnresponsive`'))
##' waldTest(zlmVbeta, Hypothesis('`PopulationCD154+VbetaUnresponsive` - `PopulationCD154-VbetaUnresponsive`'))
setClass('ZlmFit', slots=list(coefC='matrix', coefD='matrix', vcovC='array', vcovD='array', LMlike='LMlike', sca='SummarizedExperiment0', deviance='matrix', loglik='matrix', df.null='matrix', df.resid='matrix', dispersion='matrix', dispersionNoshrink='matrix', priorDOF='numeric', priorVar='numeric', converged='matrix', hookOut='ANY'))

##' An S4 class for Gene Set Enrichment output
##'
##' This holds output from a call to gseaAfterBoot.
##' It primarily provides a summary method.
##' @slot tests array: gene sets X {discrete,continuous} X {stat, variance, degrees of freedom, avg correlation} X {test, null}
##' @slot bootR number of bootstrap replicates
##' @seealso gseaAfterBoot
##' @seealso calcZ
##' @seealso summary,GSEATests-method
setClass('GSEATests', slots=list(tests='array', bootR='numeric'))
