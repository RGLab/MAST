##' @import Biobase
##' @importFrom plyr rbind.fill
##' @import methods
##' @include AllGenerics.R
NULL



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
##' @seealso \code{\link{vbeta}}, \code{\link{FluidigmAssay}}
NULL

Mandatory_Featurevars <- NULL#c('primerid')
Mandatory_Cellvars <- NULL#c('wellKey')

##' @import SummarizedExperiment
##' @importFrom GenomicRanges colData
##' @import S4Vectors
##' @importFrom S4Vectors mcols
setClass('SingleCellAssay', contains='SummarizedExperiment0')

setClass('Mapping', contains='list')
setMethod('initialize', 'Mapping', function(.Object, keys=NULL, values=NULL, ...){
  .Object <- callNextMethod(.Object, ...)
  if(!is.null(keys)){
    if(!is.null(values)) values <- rep(NA, length(keys))
    if(!is.character(keys)) stop('keys must be character')
    .Object@.Data <- vector(mode='list', length=length(keys))
    names(.Object@.Data) <- keys
    for(i in seq_along(.Object@.Data)) .Object@.Data[[i]] <- values[[i]]
  }
  
  .Object
})

##' @describeIn show
setMethod('show', 'Mapping', function(object){
  cat(class(object), ' containing : ', names(object), '\n')
})


##'SCASet is a set of SingleCellAssay objects or objects of its subclasses (i.e. FluidigmAssay)
##'The constructor \code{SCASet} should be used to make objects of this class.
##' @slot set: A \code{list} of \code{SingleCellAssays} or its subclasses
##' @seealso SCASet
setClass("SCASet",
         representation=list(set="list"),validity=function(object){
           if(all(names(object@set)!=unlist(lapply(object@set,function(x) x@id),use.names=FALSE))){
             warning("Names of the SCASet don't match the SingleCellAssay id's. Plese use the SingleCellAssay() constructor.")
             return(FALSE)
           }
           return(TRUE)
         })


## Classes
##' Linear Model-like Class
##'
##' Wrapper around modeling function to make them behave enough alike that Wald tests and Likelihood ratio are easy to do.
##' To implement a new type of zero-inflated model, extend this class.
##'
##' @section Slots:
##' \describe{
##' \item{design}{a data.frame from which variables are taken for the right hand side of the regression}
##' \item{fitC}{The continuous fit}
##' \item{fitD}{The discrete fit}
##' \item{response}{The left hand side of the regression}
##' \item{fitted}{A \code{logical} with components "C" and "D", TRUE if the respective component has converge}
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
         slots=c(design='DataFrame', modelMatrix='matrix', fitC='ANY', fitD='ANY', response='ANY', fitted='logical', formula='formula', fitArgsD='list', fitArgsC='list', priorVar='numeric', priorDOF='numeric',
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
##' @slot weightFun function to map expression values to probabilities of expression
setClass('GLMlike', contains='LMlike', slots=c(weightFun='function'), prototype=list(weightFun=function(x){
    ifelse(x>0, 1, 0)
}))

##' Initialize a prior to be used a prior for BayeGLMlike/BayesGLMlike2
##'
##' @param names character vector of coefficients.  The `(Intercept)` will be ignored.
##' @return 3d array, with leading dimension giving the prior 'loc'ation, 'scale' and degrees of freedom (df),
##' second dimension giving the component ('C'ontinuous or 'D'iscrete)
##' and trailing dimension giving the coefficient to which the prior applies.
##' The location is initialized to be 0, the scale to 2, and degrees of freedom of 1, following the default of \link{bayesglm}
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
setClass('LMERlike', contains='LMlike', slots=c(pseudoMM='data.frame'), validity=function(object){
    if(length(object@response)>0 & nrow(object@pseudoMM)>0){
        stopifnot(nrow(object@pseudoMM)==length(object@response))
    }
    if(object@priorDOF!=0) stop('Empirical bayes shrinkage not implemented for lmer/glmer.')
    })

setClass('ConstrainedGLMlike', contains='LMlike')
setClass('RidgeBGLMlike',contains="BayesGLMlike",slots=c(lambda='numeric'),prototype = list(lambda=0.1) )

## Ways to specify hypothesis
setClass('Hypothesis', contains='character', slots=list(contrastMatrix='matrix'))
setClass('CoefficientHypothesis', contains='Hypothesis', slots=list(index='numeric'))

##' An S4 class to hold the output of a call to zlm
##'
##' This holds output from a call to zlm.SingleCellAssay.  Many methods are defined to operate on it.  See below.
##' @slot coefC
##' @slot coefD matrices of coefficients
##' @slot vcovC
##' @slot vcovD array of variance/covariance matrices for coefficients
##' @slot LMlike the LmWrapper object used
##' @slot sca the \code{SingleCellAssay} object used
##' @slot deviance
##' @slot loglik
##' @slot df.null
##' @slot df.resid
##' @slot dispersion
##' @slot dispersionNoShrink
##' @slot priorDOF
##' @slot priorVar
##' @slot converged output that may optionally be set by the underlying modeling function
##' @slot hookOut a list of length ngenes containing output from a hook function, if \code{zlm} was called with one
##' @seealso zlm.SingleCellAssay summary,ZlmFit-method
setClass('ZlmFit', slots=list(coefC='matrix', coefD='matrix', vcovC='array', vcovD='array', LMlike='LMlike', sca='SummarizedExperiment0', deviance='matrix', loglik='matrix', df.null='matrix', df.resid='matrix', dispersion='matrix', dispersionNoshrink='matrix', priorDOF='numeric', priorVar='numeric', converged='matrix', hookOut='ANY'))


##' Constructor for a FluidigmAssay


##' Constructs a SCASet
##'
##' An SCASet is a list of SingleCellAssays or objects inheriting from SingleCellAssay. The type of constructor called is determined by the value of contentClass, which should be the class of the SCA inheriting object contained in this SCASet. Both the class and the constructor should exist and have the same name. The code dynamically looks to see if the a function with the same name exists, and ASSUMES it is the constructor for the class.
##' @param dataframe flat data.frame ala SingleCellAssay
##' @param splitby either a character vector naming columns or a factor or a list of factors used to split dataframe into SingleCellAssays
##' @param idvars character vector naming columns that uniquely identify a cell
##' @param primerid character vector of length 1 that names the column that containing what feature was measured
##' @param measurement character vector of length 1 that names the column containing the measurement
##' @param contentClass a character, the name of the class being constructed within this SCASet. Defaults to SingleCellAssay. Other methods may pass in other classes, i.e. FluidigmAssay.
##' @param ... passed up to SingleCellAssay or other dynamically called constructor.
##' @return SCASet
##' @note The dynamic lookup of the constructor could be made more robust. 
##' @export 
SCASet<-function(dataframe,splitby,idvars=NULL,primerid=NULL,measurement=NULL,contentClass="SingleCellAssay",...){
  if(is.character(splitby) && all(splitby %in% names(dataframe))){
  spl<-split(dataframe,dataframe[, splitby])
} else if(is.factor(splitby) || is.list(splitby) || is.character(splitby)){
  spl <- split(dataframe, splitby)
} else{
  stop("Invalid 'splitby' specification")
}
 
  set<-vector("list",length(spl))
  names(set)<-names(spl)
  for(i in seq_along(set)){
    ##construct a call using contentClass
    F <- try(getFunction(contentClass),silent=TRUE)
    if(is(F,"try-error"))
      message("Can't construct a class of type ",contentClass[[1]],". Constructor of this name doesn't exist")
      cl<-as.call(list(as.name(contentClass[[1]]),dataframe=spl[[i]],idvars=idvars,primerid=primerid,id=names(spl)[[i]], measurement=measurement,...))
    set[[i]]<-eval(cl)
  }
  new("SCASet",set=set)
}
