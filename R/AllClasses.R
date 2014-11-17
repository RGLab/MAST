##' @import Biobase
##' @import BiocGenerics
##' @importFrom plyr rbind.fill
##' @import methods
##' @include AllGenerics.R
NULL

##' DataLayer class
##' 
##' DataLayer is a 3-D array, wrapped to make it look like a matrix.
##' It is used to hold matrix-like expression data, for which we might want to keep several representations (transformations) around.
##' The number of matrix "layers" is given by the trailing dimension.
##' Dimensions 1 and 2 correspond to the "rows" and "columns" of the matrix.
##' The layer that is active can be set, and additional layers created (concatenated).
##' }
##' \section{Slots}{
##' DataLayer extends array, and has the following additional slots
##' \describe{
##'   \item{.Data}{the underlying array}
##'   \item{valid}{a \code{logical} that may optionally indicate the freshness of derived layers (if the underlying data changes).  Not currently used.}
##'   \item{layer}{which 'slice' of the array is being used}
##' }}
##' \section{Methods}{
##' \describe{
##' \item{addlayer}{Concatentate another slice onto the object}
##' \item{layername}{Return the name of the current slice}
##' \item{layer}{Return the active layer}
##' \item{layer<-}{Set the active layer}
##' \item{exprs}{Return the matrix representation of the active layer}
##' \item{exprs<-}{Replace the matrix on the current layer.}
##' }
##' @examples
##' ar <- array(1:10, dim=c(2, 5, 1))
##' dl <- new('DataLayer', .Data=ar)
##' nrow(dl) #2
##' ncol(dl) #5
##' layer(dl)
##' dl <- addlayer(dl, 'negative')
##' ex <- exprs(dl)
##' layer(dl) <- 'negative' #or could use 2
##' exprs(dl)<- -ex
##' exprs(dl)
##' @name DataLayer-class
##' @docType class 
##' @aliases DataLayer
##' @seealso \code{\link{SingleCellAssay}}, \code{\link{SingleCellAssay-class}}
setClass('DataLayer', contains='array', representation=representation(layer='numeric', valid='logical'), prototype=prototype(array(NA, dim=c(0, 0, 1)), layer=1L, valid=TRUE), validity=function(object){
  #cat('DL dim ', dim(object@.Data), '\n')
  length(dim(object@.Data))==3
  })

setClass('Mapping', contains='list')
setMethod('initialize', 'Mapping', function(.Object, keys=NULL, values=NULL, ...){
  .Object <- callNextMethod()
  if(!is.null(keys)){
    if(is.null(values)) values <- rep(NA, length(keys))
    if(!is.character(keys)) stop('keys must be character')
    .Object@.Data <- vector(mode='list', length=length(keys))
    names(.Object@.Data) <- keys
    for(i in seq_along(.Object@.Data)) .Object@.Data[[i]] <- values[[i]]
  }
  
  .Object
})

setMethod('show', 'Mapping', function(object){
  cat(class(object), ' containing : ', names(object), '\n')
})


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



SingleCellAssayValidity <- function(object){
  ##message('SingleCellAssayValidity') #DEBUG
  if(nrow(cData(object))==0 || nrow(fData(object) == 0)) return(TRUE)
  if(nrow(object)!=nrow(cData(object))){
    message('dimension mismatch between cData and nrows')
    return(FALSE)
    }
  if(ncol(object)!=nrow(fData(object))){
    message('dimension mismatch between fData and ncols')
    return(FALSE)
  }

  if(!all(fData(object)$primerid == colnames(object))){
    message("'DataLayer' column names mismatch featureData 'primerid' field")
    return(FALSE)
  }
  
  if(!all(cData(object)$wellKey == row.names(object))){
    message("'DataLayer' row names mismatch cellData 'wellKey' field")
    return(FALSE)
  }
  
  if(!all(names(object@cmap) %in% names(cData(object)))){
    message('some expected fields in cData are missing')
    return(FALSE)
  }
  if(!all(names(object@fmap) %in% names(fData(object)))){
    message('some expected fields in fData are missing')
    return(FALSE)
  }
  TRUE                                  #this stuff might not belong in the validity, it's getting called too early when subclasses of SingleCellAssay are constructed
}
                                          

##' SingleCellAssay class
##' 
##' SingleCellAssay represents an arbitrary single cell assay
##' It is meant to be flexible and is subclassed to represent specific assay
##' types like Fluidigm and NanoString. It should be constructed using the \code{SingleCellAssay}, \code{SCASet} or subclass constructors.
##' mapNames for the SingleCellAssay class are in the object \code{SingleCellAssay:::Mandatory_Cellvars}
##' mapNames for the FluidigmAssay class are in the object \code{SingleCellAssay:::FluidigmMapNames}
##' }
##' \section{Slots}{
##' SingleCellAssay extends class \code{\link{DataLayer}}, so inherits its slots and methods.  It also contains the following additional slots:
##' \describe{
##'   \item{featureData}{an \code{AnnotatedDataFrame} that describes feature-level metadata (i.e. genes)}
##'   \item{phenoData}{an \code{AnnotatedDataFrame} that describes the phenotype-level metadata (i.e. subject or experimental unit)} (not yet implemented)
##'   \item{cellData}{an \code{AnnotatedDataFrame} that describes the cell-level metadata (i.e. per individual cell)}
##'   \item{description}{a \code{data.frame}}
##'   \item{id}{a vector of type \code{character} that identifies the set of columns acting as a primary key to uniquely identify a single-cell or single-well across all wells / cells / assays / subjects / conditions in the data set.}
##' }
##' @name SingleCellAssay-class
##' @docType class 
##' @aliases SingleCellAssay-class
##' @aliases FluidigmAssay-class
##' @aliases NanoStringAssay-class
##' @aliases show,SingleCellAssay-method
##' @rdname SingleCellAssay-class
##' @seealso \code{\link{SingleCellAssay}}, \code{\link{NanoStringAssay}}, \code{\link{FluidigmAssay}}, \code{\link{DataLayer}}
setClass("SingleCellAssay",contains="DataLayer",
         representation=representation(featureData="AnnotatedDataFrame",
           phenoData="AnnotatedDataFrame",
           cellData="AnnotatedDataFrame",
           description='data.frame',
           id="ANY",
           cmap='Mapping', fmap='Mapping',
           keep.names='logical'),
         prototype=prototype(phenoData=new("AnnotatedDataFrame"),
           featureData=new("AnnotatedDataFrame"),
           cellData=new("AnnotatedDataFrame"),
           description=data.frame(),
           id=numeric(0),
           cmap=new('Mapping', keys=Mandatory_Cellvars),
           fmap=new('Mapping', keys=Mandatory_Featurevars),
           keep.names=TRUE),
         validity=SingleCellAssayValidity)


## Same as SingleCellAssay, but with additional mapNames
FluidigmMapNames <- c(Mandatory_Cellvars, 'ncells')

setClass('FluidigmAssay', contains='SingleCellAssay', prototype=prototype(cmap=new('Mapping', keys=FluidigmMapNames)),validity=SingleCellAssayValidity)

#Could write a constructor that takes a post-processing function...
setClass('NanoStringAssay', contains='FluidigmAssay',validity=SingleCellAssayValidity)


##'Holds output and diagnostics from thresholdNanoString
##'Not intended to be called by the user.
##' 
##' @section Slots:
##' \describe{
##' \item{melted}{A \code{data.frame} containing a melted version of \code{nsa}, plus the columns 'ps', giving the probability that a measurement belongs to the signal cluster, 'clusterID' the inferred cluster}
##' \item{nsa}{The thresholded \code{NanoStringAssay} with the thresholded expression in layer \code{et}}
##' \item{densities}{A \code{list} of length \code{ncol(nsa)} of marginal (mixture model) densities of each gene.}
##' \item{means}{A \code{matrix} dimension \code{ncol(nsa)} \eqn{\times} 2 given the posterior mean of each cluster.}
##' \item{props}{A \code{matrix} dimension \code{ncol(nsa)} \eqn{\times} 2 given the posterior probability of each cluster.}
##' \item{startLayer}{A \code{character} giving the initial layer that was used to generate the thresholding}
##' }
##' @seealso thresholdNanoString
##' @docType class
setClass('ThresholdedNanoString', representation=representation(melted='data.frame', nsa='NanoStringAssay', densities='list', means='matrix', props='matrix', startLayer='character'))



##'RNASeqAssay class. Doesn't require ncells
##'@exportClass RNASeqAssay
setClass('RNASeqAssay',contains='SingleCellAssay', prototype=prototype(cmap=new('Mapping',keys=Mandatory_Cellvars)),validity=SingleCellAssayValidity)

##'SCASet is a set of SingleCellAssay objects or objects of its subclasses (i.e. FluidigmAssay)
##'The constructor \code{SCASet} should be used to make objects of this class.
##' }
##' \section{Slots}{
##' \describe{
##' \item{set}{A \code{list} of \code{SingleCellAssays} or its subclasses}
##' }
##' 
##' @rdname SCASet-class
##' @docType class
##' @name SCASet-class
##' @exportClass SCASet
##' @aliases SCASet-class
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
##' @seealso fit
##' @seealso coef
##' @seealso lrTest
##' @seealso waldTest
##' @seealso vcov
##' @seealso dof
##' @seealso logLik
##' @name LMlike-class
##' @docType class
setClass('LMlike',
         slots=c(design='data.frame', modelMatrix='matrix', fitC='ANY', fitD='ANY', response='ANY', fitted='logical', formula='formula', fitArgsD='list', fitArgsC='list', priorVar='numeric', priorDOF='numeric',
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

setClass('GLMlike', contains='LMlike')

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

setClass('BayesGLMlike', contains='GLMlike', slots=c(coefPrior='array'),
         prototype=list(prior=defaultPrior(character(0))),
         validity=function(object){
             ## if(length(object@coefPrior>0))
             ##     if(dim(object@coefPrior)[3] != sum(colnames(model.matrix(object))!='(Intercept)')) stop('prior must have same number of components as model.matrix')
             TRUE
         })
setClass('BayesGLMlike2', contains='BayesGLMlike')
setClass('BayesGLMlikeWeight', contains='BayesGLMlike')

setClass('LMERlike', contains='LMlike', slots=c(pseudoMM='data.frame'), validity=function(object){
    if(length(object@response)>0 & nrow(object@pseudoMM)>0){
        stopifnot(nrow(object@pseudoMM)==length(object@response))
    }
    if(object@priorDOF!=0) stop('Empirical bayes shrinkage not implemented for lmer/glmer.')
    })

setClass('ConstrainedGLMlike', contains='LMlike')


## Ways to specify hypothesis
setClass('Hypothesis', contains='character', slots=list(transformed='matrix'))
setClass('CoefficientHypothesis', contains='character', slots=list(transformed='numeric'))

setClass('ZlmFit', slots=list(coefC='matrix', coefD='matrix', vcovC='array', vcovD='array', LMlike='LMlike', sca='SingleCellAssay', deviance='matrix', loglik='matrix', df.null='matrix', df.resid='matrix', dispersion='matrix', dispersionNoshrink='matrix', priorDOF='numeric', priorVar='numeric', converged='matrix', hookOut='ANY'))

##' SingleCellAssay: A constructor for an object of type SingleCellAssay.
##'
##' This is the constructor for the class. This class intends to ease the analysis of single cell assays, in which multiple, exchangeable, cells from an experimental unit (patient, or organism) are assayed along several (or many) dimensions, such as genes. A few examples of this might be Fluidigm gene expression chips, or single cell sequencing experiments.  The chief functionality is to make it easy to keep cellular-level metadata linked to the measurements through \code{cellData} and \code{phenoData}.  There are also subsetting and splitting measures to coerce between a SingleCellAssay, and a \link{SCASet}.
##' @param dataframe A 'flattened' \code{data.frame} or \code{data.table} containing columns giving cell and feature identifiers and  a measurement column
##' @param idvars character vector naming columns that uniquely identify a cell
##' @param primerid character vector of length 1 that names the column that identifies what feature (i.e. gene) was measured
##' @param measurement character vector of length 1 that names the column containing the measurement 
##' @param id An identifier (eg, experiment name) for the resulting object
##' @param cellvars Character vector naming columns containing additional cellular metadata
##' @param featurevars Character vector naming columns containing additional feature metadata
##' @param phenovars Character vector naming columns containing additional phenotype metadata
##' @param ... additional arguments are ignored
##' @export SingleCellAssay
##' @aliases SingleCellAssay
##' @name SingleCellAssay
##' @seealso \code{\link{FluidigmAssay}}
##' @docType methods
##' @examples
##' ## See FluidigmAssay for examples
##' \dontrun{example(FluidigmAssay)}
##' @return SingleCellAssay object
SingleCellAssay<-function(dataframe=NULL,idvars=NULL,primerid=NULL,measurement=NULL,id=numeric(0), cellvars=NULL, featurevars=NULL, phenovars=NULL, ...){
  new('SingleCellAssay', dataframe=dataframe, idvars=idvars, primerid=primerid, measurement=measurement, id=id, cellvars=cellvars, featurevars=featurevars, phenovars=phenovars)
}

checkArrayNames <- function(exprsArray, fData, cData){
    if(!is.numeric(exprsArray)) stop('`exprsArray` must be numeric')
    if(length(dim(exprsArray))<2) stop('`exprsArray` must be matrix or array')
    if(length(dim(exprsArray))<3) dim(exprsArray) <- c(dim(exprsArray), 1)
    dn <- dimnames(exprsArray)[1:2]
    dl <- new('DataLayer', .Data=exprsArray)

    pidDefault <- sprintf('p%0*d', ceiling(log10(ncol(dl)+1)), seq_len(ncol(dl)))
    wkDefault <- sprintf('wk%0*d', ceiling(log10(nrow(dl)+1)), seq_len(nrow(dl)))
    
    if(missing(fData)) fData <- data.frame(primerid=pidDefault, stringsAsFactors=FALSE)
    if(missing(cData)) cData <- data.frame(wellKey=wkDefault,  stringsAsFactors=FALSE)
    
    
    if(nrow(dl) != nrow(cData)) stop('`cData` must contain as many rows as `exprsArray`')
    if(ncol(dl) != nrow(fData)) stop('`fData` must contain as many columns as `exprsArray`')
    
    if(!('primerid' %in% names(fData))){
        warning("`fData` has no primerid.  I'll make something up.")
        fData$primerid <- pidDefault
    }
    row.names(fData) <- fData$primerid

    if(!('wellKey' %in% names(cData))){
        warning("`cData` has no wellKey.  I'll make something up.")
        cData$wellKey <- wkDefault
    }
    row.names(cData) <- wkDefault
    

    if(is.null(dn) || is.null(dn[[1]]) || is.null(dn[[2]])){
        message('No dimnames in `exprsArray`, assuming `fData` and `cData` are sorted according to `exprsArray`')
        dn <- list(wellkey=row.names(cData), primerid=row.names(fData), measure='et')            
    }
    if(!isTRUE(all.equal(dn[[1]], cData$wellKey))) stop('Order of `exprsArray` and `cData` doesn\'t match')
    if(!isTRUE(all.equal(dn[[2]], fData$primerid))) stop('Order of `exprsArray` and `fData` doesn\'t match')
    dimnames(dl) <- dn
    fData <- as(fData, 'AnnotatedDataFrame')
    cData <- as(cData, 'AnnotatedDataFrame')
    list(exprsArray=dl, fData=fData, cData=cData)
}

##' Construct a SingleCellAssay from a matrix or array of expression
##'
##' If 
##' @param class What class of object are we constructing?
##' @param exprsArray matrix or array, rows are cells, columns are genes
##' @param cData cellData data.frame or AnnotatedDataFrame
##' @param fData featureData data.frame or AnnotatedDataFrame
##' @return an object of class \code{class}
##' @examples
##' ncells <- 10
##' ngenes <- 5
##' fData <- data.frame(primerid=LETTERS[1:ngenes])
##' cData <- data.frame(wellKey=seq_len(ngenes))
##' mat <- matrix(rnorm(ncells*ngenes), nrow=ngenes)
##' sca <- FromMatrix('SingleCellAssay', mat, cData, fData)
##' stopifnot(inherits(sca, 'SingleCellAssay'))
FromMatrix <- function(class, exprsArray, cData, fData){
    can <- checkArrayNames(exprsArray, fData, cData)
    dl <- can$exprsArray
    fData <- can$fData
    cData <- can$cData
    new(class, .Data=dl, cellData=cData, featureData=fData, sort=FALSE)
}

##' Constructor for a FluidigmAssay
##'
##' Constructs a FluidigmAssay object. Differs little from the SingleCellAssay constructor. Only the \code{ncells} parameter is additionally required.
##' @inheritParams SingleCellAssay
##' @param ncells A \code{character} specifying the column which gives the number of cells per well
##' @param geneid An optional \code{character} alternate id for primers.
##' @return A FluidigmAssay object
##' @author Andrew McDavid and Greg Finak
##' @examples
##' data(vbeta)
##' colnames(vbeta)
##' vbeta <- computeEtFromCt(vbeta)
##' vbeta.fa <- FluidigmAssay(vbeta, idvars=c("Subject.ID", "Chip.Number", "Well"), primerid='Gene', measurement='Et', ncells='Number.of.Cells', geneid="Gene",  cellvars=c('Number.of.Cells', 'Population'), phenovars=c('Stim.Condition','Time'), id='vbeta all')
##' show(vbeta.fa)
##' nrow(vbeta.fa)
##' ncol(vbeta.fa)
##' head(fData(vbeta.fa)$primerid)
##' table(cData(vbeta.fa)$Subject.ID)
##' vbeta.sub <- subset(vbeta.fa, Subject.ID=='Sub01')
##' show(vbeta.sub)
##' @export
FluidigmAssay<-function(dataframe=NULL,idvars,primerid,measurement, ncells, geneid=NULL,id=numeric(0), cellvars=NULL, featurevars=NULL, phenovars=NULL, ...){
  cmap <- new('Mapping', .Data=list('ncells'=ncells))
    new('FluidigmAssay', dataframe=dataframe, idvars=idvars, primerid=primerid, measurement=measurement, id=id, cellvars=cellvars, featurevars=featurevars, phenovars=phenovars, cmap=cmap)
}

##' Constructs a SCASet
##'
##' An SCASet is a list of SingleCellAssays or objects inheriting from SingleCellAssay. The type of constructor called is determined by the value of contentClass, which should be the class of the SCA inheriting object contained in this SCASet. Both the class and the constructor should exist and have the same name. The code dynamically looks to see if the a function with the same name exists, and ASSUMES it is the constructor for the class.
##' ##' ##' TODO SCASet constructor should perhaps take a SingleCellAssay class or FluidigmClass rather than a dataframe. Then we can learn the class type for construction.
##' @title SCASet constructor
##' @param dataframe flat data.frame ala SingleCellAssay
##' @param splitby either a character vector naming columns or a factor or a list of factors used to split dataframe into SingleCellAssays
##' @param idvars character vector naming columns that uniquely identify a cell
##' @param primerid character vector of length 1 that names the column that containing what feature was measured
##' @param measurement character vector of length 1 that names the column containing the measurement
##' @param contentClass a character, the name of the class being constructed within this SCASet. Defaults to SingleCellAssay. Other methods may pass in other classes, i.e. FluidigmAssay.
##' @param ... passed up to SingleCellAssay or other dynamically called constructor.
##' @return SCASet
##' @note The dynamic lookup of the constructor could be made more robust. 
##' @aliases SCASet
##' @rdname SCAset-methods
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
