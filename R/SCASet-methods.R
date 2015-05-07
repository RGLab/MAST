setAs("SCASet","list",function(from)from@set)

setAs('list', 'SCASet', function(from){
    for(it in from){
        if(!is(it, 'SingleCellAssay')) stop("All members of 'x' must inherit from 'SingleCellAssay'")
    }
    new('SCASet', set=from)
})

##' Apply a function to each member of a SCASet
##'
##' @param X SCASet
##' @param FUN a function
##' @param ... passed to FUN
##' @importMethodsFrom BiocGenerics lapply
##' @export lapply
##' @aliases lapply,SCASet-method
setMethod("lapply",c("SCASet"),function(X,FUN,...){
  FUN<-match.fun(FUN)
  X<-as(X,"list")
  .Internal(lapply(X,FUN))
})



##' Index into an SCASet object
##'
##' Return a \code{SingleCellAssay} at index i
##' @export
##' @param x object to be subscripted
##' @param i index
##' @param j Ignored
##' @param ... Ignored
##' @return subscripted SingleCellAssay or derived class
##' @aliases [[,SCASet,ANY-method
##' @aliases [[,SCASet-method
##' @exportMethod [[
# @details \code{signature(x="SCASet", i="ANY")}: \code{x[[i]]}, where \code{i} is length-1 integer or character matching sampleNames.  Returns the SingleCellAssay at position or with sampleName \code{i}.
setMethod("[[",signature("SCASet","ANY"),function(x,i,j,...){
  if(length(i)!=1)
    stop("subscript out of bounds (index must have length 1)")
  sca<-x@set[[if(is.numeric(i))
              sampleNames(x)[[i]]
  else i]]
  return(sca)
})

##' Subset an SCASet to a smaller SCASet
##' @export
##' @aliases [,SCASet,ANY,ANY-method
##' @aliases [,SCASet-method
##' @param x \code{SCASet} to be subscripted
##' @param i index
##' @param j ignored
##' @param ... ignored
##' @param drop ignored
##' @return An \code{SCASet}
##' @details \code{signature(x="SCASet",i="ANY")}: \code{x[i]}, where \code{i} is the vector of integers or characters matching sampleNames. Returns an SCASet of \code{length(i)} of elements at positions in \code{i}.
setMethod("[",signature("SCASet","ANY"),function(x,i,j,...drop=FALSE){
    ret<-x
    ret@set <- x@set[if(is.numeric(i))
                     sampleNames(x)[i]] 
    return(ret)
})

##'get the sample names in an SCASet
##' @export
##' @aliases sampleNames,SCASet-method
##' @param object The \code{SCASet} object
##' @importMethodsFrom Biobase sampleNames
##' @importMethodsFrom Biobase "sampleNames<-"
##' @return a \code{character} vector of the names of the samples in the \code{SCASet}
setMethod("sampleNames","SCASet",function(object){
  unlist(lapply(object@set,function(x)x@id),use.names=FALSE)
})

##' @describeIn show
setMethod("show","SCASet",function(object){
  cat("SCASet of size ",length(object@set),"\n")
  cat("Samples ",paste(unlist(lapply(object,function(x)x@id),use.names=T),collapse=", "),"\n")
  invisible(NULL)
})

##' Get the length of an SCASet
##'
##' @param x SCASet
##' @return numeric
##' @aliases length,SCASet-method
setMethod('length', 'SCASet', function(x){
  length(x@set)
})
